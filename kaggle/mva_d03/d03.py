"""Private Kaggle kernel santanubanerjee9/mva-d03: low-pass depth for mosaic aneuploidy (package D03).

PATIENT DATA (FASTQ). Runs only as a private mva-* kernel on santanubanerjee9 (CLAUDE.md, approvals of
2026-10-14 and 2026-10-18: FASTQ steps run only inside private Kaggle/HF compute, never on the laptop;
streaming the head of one FASTQ pair here approved by Santanu in a live chat on 2026-10-24).

Self-contained, so it runs from the Kaggle editor with no inputs attached: Kaggle secrets reach only
runs started from the website, not CLI pushes (seen 2026-10-24). The counting, GC and ratio functions
mirror pipeline/lowpass_depth.py (commit 703160b).

  1. Read the Kaggle secret HF_TOKEN (Santanu's own fine-grained, read-only token). It goes to curl
     through the environment only and is never printed or written anywhere.
  2. Stream the first N_PAIRS read pairs of lane 1 (R1 and R2) from the gated dataset into /tmp.
  3. GRCh38 analysis set from NCBI (md5 checked against NCBI's list), the file mva-refs-build uses.
  4. minimap2 -ax sr; primary reads with MAPQ >= 20 counted in 100-kb bins; GC per bin; GC-corrected
     per-chromosome depth ratios against the autosomes.
Kept in /kaggle/working/d03: ratios for the 22 autosomes only (X and Y are left out, because their
expected depth depends on sex, which this pipeline never infers) and summary.json. Reads, alignments
and bin counts stay in /tmp.
"""
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import time
import urllib.request
from collections import Counter, defaultdict
from statistics import median

REPO = "SageBio/mva-hackathon-2026-data"
NCBI = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
        "seqs_for_alignment_pipelines.ucsc_ids/")
FNA = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
N_PAIRS, BIN = 4_000_000, 100_000
W, TMPD, ENV = "/kaggle/working/d03", "/tmp/d03", "/tmp/env"
os.makedirs(W, exist_ok=True)
os.makedirs(TMPD, exist_ok=True)
t0 = time.time()
summary = {"started": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()), "n_pairs_requested": N_PAIRS}


def save():
    summary["minutes"] = round((time.time() - t0) / 60, 1)
    json.dump(summary, open(f"{W}/summary.json", "w"), indent=1)


def sh(cmd, check=True, env=None, pipefail=True):
    print("+", cmd[:160], flush=True)
    r = subprocess.run(["bash", "-c", ("set -o pipefail; " if pipefail else "") + cmd], text=True,
                       capture_output=True, env=env)
    if r.returncode:
        print(r.stderr[-1500:], flush=True)
        if check:
            summary["failed_at"] = cmd[:100]
            save()
            raise SystemExit(f"failed ({r.returncode})")
    return r


def is_autosome(chrom):
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return c.isdigit() and 1 <= int(c) <= 22


def gc_bins(lines, bin_size=BIN, min_called=0.5):
    chrom, pos, gc, called = None, 0, 0, 0
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if chrom is not None and pos % bin_size:
                yield chrom, pos // bin_size * bin_size, (gc / called if called >= min_called * (pos % bin_size) else None)
            chrom, pos, gc, called = line[1:].split()[0], 0, 0, 0
            continue
        line = line.upper()
        while line:
            take = min(len(line), bin_size - pos % bin_size)
            seg, line = line[:take], line[take:]
            g = seg.count("G") + seg.count("C")
            gc += g
            called += g + seg.count("A") + seg.count("T")
            pos += take
            if pos % bin_size == 0:
                yield chrom, pos - bin_size, (gc / called if called >= min_called * bin_size else None)
                gc, called = 0, 0
    if chrom is not None and pos % bin_size:
        yield chrom, pos // bin_size * bin_size, (gc / called if called >= min_called * (pos % bin_size) else None)


def chromosome_ratios(bins):
    strata = defaultdict(list)
    for b in bins:
        if b["gc"] is not None and is_autosome(b["chrom"]) and b["reads"] > 0:
            strata[round(b["gc"], 2)].append(b["reads"])
    medians = {g: median(v) for g, v in strata.items() if len(v) >= 3}
    per_chrom = defaultdict(list)
    for b in bins:
        m = medians.get(round(b["gc"], 2)) if b["gc"] is not None else None
        if m:
            per_chrom[b["chrom"]].append(b["reads"] / m)
    base = median([x for c, v in per_chrom.items() if is_autosome(c) for x in v])
    rows = [{"chrom": c, "bins": len(v), "ratio": round(median(v) / base, 4),
             "trisomy_fraction_if_gain": round(max(0.0, 2 * (median(v) / base - 1)), 3),
             "monosomy_fraction_if_loss": round(max(0.0, 2 * (1 - median(v) / base)), 3)}
            for c, v in per_chrom.items() if is_autosome(c)]
    med = median(r["ratio"] for r in rows)
    mad = median(abs(r["ratio"] - med) for r in rows) * 1.4826 or 1e-4
    for r in rows:
        r["robust_z"] = round((r["ratio"] - med) / mad, 2)
    return sorted(rows, key=lambda r: int(r["chrom"].removeprefix("chr"))), med, mad


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
summary["secret"] = {"attached": bool(token), "errors": errors}
save()
if not token:
    raise SystemExit("HF_TOKEN secret not available to this run")
henv = dict(os.environ, HF_TOKEN=token)
req = urllib.request.Request(f"https://huggingface.co/api/datasets/{REPO}/tree/main?recursive=true",
                             headers={"Authorization": f"Bearer {token}"})
with urllib.request.urlopen(req, timeout=60) as r:
    tree = [x for x in json.loads(r.read().decode()) if x.get("type") == "file"]
fastqs = sorted(x["path"] for x in tree if ".fastq" in x["path"] or ".fq" in x["path"])
r1 = next(p for p in fastqs if "_L001_R1" in p)
r2 = r1.replace("_L001_R1", "_L001_R2")
assert r2 in fastqs, "lane 1 R2 not found"
summary["fastq"] = {"lane": "L001", "files_in_repo": len(fastqs)}
save()

# 2. Stream the head of lane 1. head closes the pipe early (SIGPIPE), so no pipefail here.
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

# 3. Toolchain and reference.
sh("mkdir -p /tmp/mm && curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /tmp/mm bin/micromamba")
sh(f"/tmp/mm/bin/micromamba create -y -q -p {ENV} -c conda-forge -c bioconda minimap2")
summary["minimap2"] = sh(f"{ENV}/bin/minimap2 --version").stdout.strip()
FA = f"{TMPD}/{FNA}"
sh(f"curl -fsSL --retry 5 -o {FA} {NCBI}{FNA}")
want = [ln.split()[0] for ln in sh(f"curl -fsSL --retry 5 {NCBI}md5checksums.txt").stdout.splitlines() if ln.strip().endswith(FNA)]
h = hashlib.md5()
with open(FA, "rb") as fh:
    for chunk in iter(lambda: fh.read(1 << 22), b""):
        h.update(chunk)
summary["reference"] = {"source": NCBI + FNA, "md5_matches_ncbi": bool(want) and want[0] == h.hexdigest()}
save()
if not summary["reference"]["md5_matches_ncbi"]:
    raise SystemExit("reference md5 mismatch")

# 4. Align and count (index built in memory; decoys kept so they absorb mismapped reads).
counts = Counter()
p = subprocess.Popen(f"{ENV}/bin/minimap2 -t 4 -ax sr {FA} {TMPD}/R1.fq {TMPD}/R2.fq 2> {TMPD}/mm2.log",
                     shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
for line in p.stdout:
    if line[0] == "@":
        continue
    f = line.split("\t", 5)
    flag, chrom, pos, mapq = int(f[1]), f[2], int(f[3]), int(f[4])
    if flag & (4 | 256 | 1024 | 2048) or chrom == "*" or mapq < 20:
        continue
    counts[(chrom, (pos - 1) // BIN * BIN)] += 1
if p.wait():
    raise SystemExit("minimap2 failed: " + open(f"{TMPD}/mm2.log").read()[-1000:])
summary["minutes_align"] = round((time.time() - t0) / 60, 1)
save()
with gzip.open(FA, "rt") as fh:
    gc = {(c, s): g for c, s, g in gc_bins(fh) if is_autosome(c)}
bins = [{"chrom": c, "start": s, "reads": n, "gc": gc.get((c, s))} for (c, s), n in counts.items() if is_autosome(c)]
rows, med, mad = chromosome_ratios(bins)
with open(f"{W}/depth_by_chromosome.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
summary["result"] = {"autosomes": len(rows), "autosomal_reads_counted": sum(b["reads"] for b in bins),
                     "ratio_median": round(med, 4), "ratio_robust_sd": round(mad, 4),
                     "detectable_trisomy_fraction_3sd": round(2 * 3 * mad, 3),
                     "flagged_robust_z_ge_3": sum(abs(r["robust_z"]) >= 3 for r in rows)}
summary["finished"] = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())
save()
shutil.rmtree(TMPD, ignore_errors=True)
print(json.dumps(summary["result"], indent=1))
print("Done. Download /kaggle/working/d03/depth_by_chromosome.tsv and summary.json.")
