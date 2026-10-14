"""Private Kaggle kernel santanubanerjee9/mva-d02: the Track 1 fast path on the child's VCF.

PATIENT DATA. Approved by Santanu in a live chat on 2026-10-14: the VCF-first files may be
processed in private Kaggle datasets and kernels named mva-* on santanubanerjee9, and everything
is deleted by 2026-11-23. Keep this kernel private, and never attach its inputs or outputs to a
kernel outside mva-*.

Inputs:
  santanubanerjee9/mva-code         code only (pipeline/, tools/)
  santanubanerjee9/mva-refs-build   kernel output: references + toolchain checks (public data)
  santanubanerjee9/mva-patient-vcf  private dataset: the VCF + index (+ proband_id.txt)
Steps (same code as the laptop route, run natively):
  1. tools/annotate_vcf.sh: bcftools norm, csq on MANE Select, ClinVar, gnomAD v4.1 slices;
  2. pipeline/annotate_table.py; 3. pipeline/track1_rank.py (rule pre-registered in methods B9);
  4. tools/validate_track1_csv.py; 5. pipeline/baf_mosaic.py (per-chromosome BAF, D03 lead).
Only small outputs are kept in /kaggle/working/d02: the ranked CSV, the rows of the ten MVA
genes plus every ranked variant, the BAF table and summary.json. The whole-genome annotated
VCF and table stay in /tmp and vanish with the session. Logs print counts, never variants.
"""
import csv
import glob
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path

W = "/kaggle/working/d02"
TMPD = "/tmp/d02"
ENV = "/tmp/env"
TOOLS = f"{ENV}/bin"
HF_USER = "SantanuBanerjee"
MVA_GENES = ["BUB1B", "CEP57", "TRIP13", "CENATAC", "SLF2", "SMC5", "MAD1L1", "MAD2L1BP", "CEP192", "BUB1"]
t0 = time.time()
summary = {"started": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())}


def sh(cmd, check=True, show=True):
    print("+", cmd.split(" ;")[0][:200], flush=True)
    r = subprocess.run(["bash", "-c", "set -o pipefail; " + cmd], text=True, capture_output=True)
    if show and r.stdout.strip():
        print(r.stdout[-3000:], flush=True)
    if r.returncode:
        print(r.stderr[-3000:], flush=True)
        if check:
            raise SystemExit(f"failed ({r.returncode})")
    return r


def one(pattern, exclude=()):
    hits = [h for h in glob.glob(pattern, recursive=True) if not any(x in h for x in exclude)]
    if len(hits) != 1:
        raise SystemExit(f"expected one match for {pattern}, found {len(hits)}")
    return hits[0]


def code_root():
    hits = glob.glob("/kaggle/input/**/tools/annotate_vcf.sh", recursive=True)
    if not hits:
        for z in glob.glob("/kaggle/input/**/repo.zip", recursive=True):
            zipfile.ZipFile(z).extractall("/tmp/code_zip")
        hits = glob.glob("/tmp/code_zip/**/tools/annotate_vcf.sh", recursive=True)
    src = os.path.dirname(os.path.dirname(hits[0]))
    shutil.copytree(src, "/tmp/code", dirs_exist_ok=True)
    return "/tmp/code"


os.makedirs(W, exist_ok=True)
os.makedirs(TMPD, exist_ok=True)

# Toolchain (bioconda, pinned to the laptop's 1.24 where available); kept if a wrapper installed it.
if not os.path.exists(f"{TOOLS}/bcftools"):
    sh("mkdir -p /tmp/mm && curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /tmp/mm bin/micromamba")
    if sh(f"/tmp/mm/bin/micromamba create -y -q -p {ENV} -c conda-forge -c bioconda bcftools=1.24 samtools=1.24 htslib=1.24", check=False).returncode:
        sh(f"/tmp/mm/bin/micromamba create -y -q -p {ENV} -c conda-forge -c bioconda bcftools samtools htslib")
summary["bcftools"] = sh(f"{TOOLS}/bcftools --version | head -1").stdout.strip()

# Inputs.
refs = os.path.dirname(one("/kaggle/input/**/refs/build.json"))
summary["refs_build"] = json.load(open(f"{refs}/build.json"))
FA = f"{refs}/GRCh38_no_alt_plus_hs38d1_analysis_set_no_chr.fa.gz"
GFF = f"{refs}/GRCh38.116.mane_select.gff3.gz"
# Kaggle unpacks uploaded .gz files (sometimes into a folder of the same name), so accept a
# bgzipped VCF with its index or a plain one, and bgzip + index into /tmp when needed.
# MVA_D02_VCF overrides the search; only the synthetic dry run (kaggle/mva_d02_synth) sets it.
found = [os.environ["MVA_D02_VCF"]] if os.environ.get("MVA_D02_VCF") else sorted(
    {h for pat in ("*.vcf.gz", "*.vcf.bgz", "*.vcf")
     for h in glob.glob(f"/kaggle/input/**/{pat}", recursive=True)
     if os.path.isfile(h) and not any(x in h for x in ("/refs/", "clinvar", "gnomad"))})
if len(found) != 1:
    raise SystemExit(f"expected one patient VCF, found {len(found)}")
src = found[0]
if src.endswith((".gz", ".bgz")) and (os.path.exists(src + ".tbi") or os.path.exists(src + ".csi")):
    vcf = src
else:
    vcf = f"{TMPD}/input.vcf.gz"
    sh(f"{'gzip -dc' if src.endswith(('.gz', '.bgz')) else 'cat'} {src} | {TOOLS}/bgzip -@4 -c > {vcf} && {TOOLS}/tabix -p vcf {vcf}", show=False)
summary["input_recompressed"] = vcf != src
pid_files = [] if os.environ.get("MVA_D02_VCF") else glob.glob("/kaggle/input/**/proband_id.txt", recursive=True)
code = code_root()

# Header checks (counts only in the log).
samples = sh(f"{TOOLS}/bcftools query -l {vcf}", show=False).stdout.split()
contigs = sh(f"{TOOLS}/bcftools view -h {vcf} | grep -c '^##contig'", show=False).stdout.strip()
chr_prefixed = sh(f"{TOOLS}/bcftools view -h {vcf} | grep -m1 '^##contig' | grep -c 'ID=chr'", check=False, show=False).stdout.strip()
summary["input"] = {"samples": len(samples), "contigs": int(contigs), "contigs_chr_prefixed": chr_prefixed == "1",
                    "records": int(sh(f"{TOOLS}/bcftools index -n {vcf}", show=False).stdout.strip() or 0)}
proband_id = open(pid_files[0]).read().strip() if pid_files else samples[0]
summary["proband_id_source"] = "proband_id.txt" if pid_files else "VCF sample name"

# 1. Annotation (whole genome, kept in /tmp only).
annotated = f"{TMPD}/annotated.vcf.gz"
env = dict(os.environ, TOOLS=TOOLS, FA=FA, GFF=GFF, CLINVAR=f"{refs}/clinvar.vcf.gz",
           GNOMAD=f"{refs}/gnomad_v4.1_joint_slices.vcf.gz", GNOMAD_BED=f"{refs}/gnomad_v4.1_covered.bed", LOW="")
r = subprocess.run(["bash", f"{code}/tools/annotate_vcf.sh", vcf, annotated], env=env, text=True, capture_output=True)
print(r.stdout[-3000:], r.stderr[-3000:], flush=True)
if r.returncode:
    raise SystemExit("annotate_vcf.sh failed")
summary["annotation_counts"] = [ln for ln in r.stdout.splitlines() if ":" in ln][-8:]
summary["minutes_annotate"] = round((time.time() - t0) / 60, 1)

# 2. Table, 3. ranking, 4. validation.
penv = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
table = f"{TMPD}/annotated_table.tsv"
subprocess.run([sys.executable, "-m", "pipeline.annotate_table", annotated, GFF, table], cwd=code, env=penv, check=True)
sys.path.insert(0, code)
from pipeline import comphet, track1_rank  # noqa: E402
from tools.validate_track1_csv import validate  # noqa: E402

rows = comphet.load(table)
panel = [ln.strip() for ln in open(f"{refs}/panel_genes_mva.txt") if ln.strip() and not ln.startswith("#")]
ranked = track1_rank.score_candidates(rows, panel)
csv_path = f"{W}/{HF_USER}_track1.csv"
track1_rank.write_csv(csv_path, proband_id, ranked)
errors, warnings = validate(Path(csv_path), hf_user=HF_USER)
summary["track1"] = {"rows": len(ranked), "validator_errors": errors, "validator_warnings": warnings,
                     "panel_genes": len(panel)}
src_counts = {}
for row in rows:
    src_counts[row.get("pop_af_source", "")] = src_counts.get(row.get("pop_af_source", ""), 0) + 1
summary["table_rows"] = len(rows)
summary["pop_af_source_counts"] = src_counts

# Small extracts for the report: every row in the ten MVA genes, plus every ranked variant.
keys = set()
for _, a, b, _ in ranked:
    for v in (a, b):
        if v is not None:
            keys.add((v["chrom"], v["pos"], v["ref"], v["alt"]))
with open(table) as fh, open(f"{W}/mva_genes_and_ranked_rows.tsv", "w", newline="") as out:
    rd = csv.DictReader(fh, delimiter="\t")
    wr = csv.DictWriter(out, fieldnames=rd.fieldnames, delimiter="\t")
    wr.writeheader()
    n = 0
    for row in rd:
        if row["gene"] in MVA_GENES or (row["chrom"], row["pos"], row["ref"], row["alt"]) in keys:
            wr.writerow(row)
            n += 1
summary["extract_rows"] = n

# 5. Per-chromosome BAF (package D03 lead; aggregate statistics only).
r = subprocess.run([sys.executable, "-m", "pipeline.baf_mosaic", vcf, f"{W}/baf_by_chromosome.tsv"],
                   cwd=code, env=penv, text=True, capture_output=True)
summary["baf"] = {"exit": r.returncode, "log_tail": r.stdout.strip().splitlines()[-3:] if r.stdout else []}

summary["minutes_total"] = round((time.time() - t0) / 60, 1)
summary["finished"] = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())
json.dump(summary, open(f"{W}/summary.json", "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in summary.items() if k != "refs_build"}, indent=1, default=str))
shutil.rmtree(TMPD, ignore_errors=True)
