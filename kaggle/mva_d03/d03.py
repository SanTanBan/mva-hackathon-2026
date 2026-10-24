"""Private Kaggle kernel santanubanerjee9/mva-d03: low-pass depth for mosaic aneuploidy (package D03).

PATIENT DATA (FASTQ). Runs only as a private mva-* kernel on santanubanerjee9 (CLAUDE.md, approvals of
2026-10-14 and 2026-10-18: FASTQ steps run only inside private Kaggle/HF compute, never on the laptop;
streaming the head of one FASTQ pair here approved by Santanu in a live chat on 2026-10-24).

Steps:
  1. Read the Kaggle secret HF_TOKEN (Santanu's own fine-grained, read-only token). It is passed to curl
     through the environment only; it is never printed or written anywhere.
  2. Stream the first N_PAIRS read pairs of lane 1 (R1 and R2) from the gated dataset into /tmp.
  3. Align with minimap2 -ax sr to the GRCh38 analysis set from mva-refs-build, and count primary reads
     with MAPQ >= 20 in 100-kb bins (pipeline/lowpass_depth.py count).
  4. GC per bin from the reference, then GC-corrected per-chromosome depth ratios (pipeline/lowpass_depth.py).
Kept in /kaggle/working/d03: the per-chromosome ratios for the 22 autosomes only (X and Y are left out,
because their expected depth depends on sex, which this pipeline never infers) and summary.json.
Reads, alignments and bin counts stay in /tmp and vanish with the session.
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile

REPO = "SageBio/mva-hackathon-2026-data"
N_PAIRS = int(os.environ.get("N_PAIRS", "4000000"))
W, TMPD, ENV = "/kaggle/working/d03", "/tmp/d03", "/tmp/env"
os.makedirs(W, exist_ok=True)
os.makedirs(TMPD, exist_ok=True)
t0 = time.time()
summary = {"started": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()), "n_pairs_requested": N_PAIRS}


def save():
    summary["minutes"] = round((time.time() - t0) / 60, 1)
    json.dump(summary, open(f"{W}/summary.json", "w"), indent=1)


def sh(cmd, check=True, env=None, pipefail=True):
    print("+", cmd[:200], flush=True)
    r = subprocess.run(["bash", "-c", ("set -o pipefail; " if pipefail else "") + cmd], text=True,
                       capture_output=True, env=env)
    if r.stdout.strip():
        print(r.stdout[-2000:], flush=True)
    if r.returncode:
        print(r.stderr[-2000:], flush=True)
        if check:
            summary["failed_at"] = cmd[:120]
            save()
            raise SystemExit(f"failed ({r.returncode})")
    return r


# 1. Secret (retried: the secrets service can refuse the first call).
token, errors = None, []
for attempt in range(4):
    try:
        from kaggle_secrets import UserSecretsClient
        token = UserSecretsClient().get_secret("HF_TOKEN")
        break
    except Exception as e:
        errors.append(f"{type(e).__name__}: {str(e)[:160]}")
        time.sleep(10)
summary["secret"] = {"attached": bool(token), "errors": errors,
                     "secrets_env_present": bool(os.environ.get("KAGGLE_USER_SECRETS_TOKEN"))}
save()
if not token:
    raise SystemExit("HF_TOKEN secret not available to this run")
henv = dict(os.environ, HF_TOKEN=token)

# File list from the Hub (paths only; the token goes in a header).
req = urllib.request.Request(f"https://huggingface.co/api/datasets/{REPO}/tree/main?recursive=true",
                             headers={"Authorization": f"Bearer {token}"})
with urllib.request.urlopen(req, timeout=60) as r:
    tree = [x for x in json.loads(r.read().decode()) if x.get("type") == "file"]
fastqs = sorted(x["path"] for x in tree if ".fastq" in x["path"] or ".fq" in x["path"])
r1 = next(p for p in fastqs if "_L001_R1" in p)
r2 = r1.replace("_L001_R1", "_L001_R2")
assert r2 in fastqs, "lane 1 R2 not found"
summary["fastq"] = {"lane": "L001", "files_in_repo": len(fastqs),
                    "sizes_gb": {k: round(next(x["size"] for x in tree if x["path"] == p) / 1e9, 2) for k, p in (("R1", r1), ("R2", r2))}}
save()

# Toolchain and inputs.
sh("mkdir -p /tmp/mm && curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /tmp/mm bin/micromamba")
sh(f"/tmp/mm/bin/micromamba create -y -q -p {ENV} -c conda-forge -c bioconda minimap2")
summary["minimap2"] = sh(f"{ENV}/bin/minimap2 --version").stdout.strip()
refs = os.path.dirname(next(iter(glob.glob("/kaggle/input/**/refs/build.json", recursive=True))))
FA = f"{refs}/GRCh38_no_alt_plus_hs38d1_analysis_set_no_chr.fa.gz"
hits = glob.glob("/kaggle/input/**/pipeline/lowpass_depth.py", recursive=True)
if not hits:
    for z in glob.glob("/kaggle/input/**/*.zip", recursive=True):
        zipfile.ZipFile(z).extractall("/tmp/code_zip")
    hits = glob.glob("/tmp/code_zip/**/pipeline/lowpass_depth.py", recursive=True)
code = os.path.dirname(os.path.dirname(hits[0]))
shutil.copytree(code, "/tmp/code", dirs_exist_ok=True)
penv = dict(os.environ, PYTHONPATH="/tmp/code", PYTHONDONTWRITEBYTECODE="1")

# 2. Stream the head of lane 1. head closes the pipe early, so curl and gzip exit with SIGPIPE: no pipefail.
for tag, path in (("R1", r1), ("R2", r2)):
    url = f"https://huggingface.co/datasets/{REPO}/resolve/main/{path}"
    sh(f'curl -sL -H "Authorization: Bearer $HF_TOKEN" "{url}" | gzip -dc 2>/dev/null | head -n {4 * N_PAIRS} > {TMPD}/{tag}.fq',
       env=henv, pipefail=False)
    n = int(sh(f"wc -l < {TMPD}/{tag}.fq").stdout.strip()) // 4
    summary.setdefault("pairs_streamed", {})[tag] = n
    save()
    if n < N_PAIRS // 2:
        raise SystemExit(f"only {n} reads streamed for {tag}")
del henv, token

# 3. Align and count (index built in memory; decoys kept so they absorb mismapped reads).
sh(f"{ENV}/bin/minimap2 -t 4 -ax sr {FA} {TMPD}/R1.fq {TMPD}/R2.fq 2> {TMPD}/mm2.log | "
   f"python3 -m pipeline.lowpass_depth count --bin 100000 > {TMPD}/counts.tsv", env=penv)
summary["minutes_align"] = round((time.time() - t0) / 60, 1)
save()

# 4. GC and ratios.
sh(f"python3 -m pipeline.lowpass_depth gc {FA} --bin 100000 > {TMPD}/gc.tsv", env=penv)
sh(f"python3 -m pipeline.lowpass_depth ratios {TMPD}/counts.tsv {TMPD}/gc.tsv {TMPD}/ratios_all.tsv", env=penv)
sys.path.insert(0, "/tmp/code")
from pipeline.common import is_autosome  # noqa: E402
import csv  # noqa: E402
from statistics import median  # noqa: E402

rows = [r for r in csv.DictReader(open(f"{TMPD}/ratios_all.tsv"), delimiter="\t") if is_autosome(r["chrom"])]
rows.sort(key=lambda r: int(r["chrom"].removeprefix("chr")))
with open(f"{W}/depth_by_chromosome.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
counted = sum(int(r["reads"]) for r in csv.DictReader(open(f"{TMPD}/counts.tsv"), delimiter="\t") if is_autosome(r["chrom"]))
ratios = [float(r["ratio"]) for r in rows]
med = median(ratios)
mad = median(abs(x - med) for x in ratios) * 1.4826
summary["result"] = {"autosomes": len(rows), "autosomal_reads_counted": counted,
                     "ratio_median": round(med, 4), "ratio_robust_sd": round(mad, 4),
                     "detectable_trisomy_fraction_3sd": round(2 * 3 * mad, 3),
                     "flagged_robust_z_ge_3": sum(abs(float(r.get("robust_z") or 0)) >= 3 for r in rows)}
summary["finished"] = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())
save()
shutil.rmtree(TMPD, ignore_errors=True)
print(json.dumps(summary["result"], indent=1))
