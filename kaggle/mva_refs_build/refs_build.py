"""Private Kaggle kernel santanubanerjee9/mva-refs-build: public references and toolchain checks.

Public data only. It builds what the patient-data kernel (kaggle/mva_d02) needs, so that kernel
can attach this one's output instead of downloading everything again:
  - GRCh38 no-alt analysis set with hs38d1 (NCBI), contigs renamed without 'chr', bgzipped and
    faidx-indexed (same recipe as tools/setup_core.sh on the laptop);
  - ClinVar GRCh38 VCF, the 2026-09-28 release used on the laptop when NCBI still archives it;
  - the gnomAD v4.1 slices from the mva-refs dataset, bgzipped and tabix-indexed (package R10);
  - the MANE Select GFF3 and the pre-registered PanelApp gene list, copied through.
Then it runs the pipeline unit tests and tools/test_annotate_e2e.py natively (MVA_NATIVE=1) on a
synthetic VCF built from public ClinVar BUB1B alleles, and records versions and results in
refs/build.json.
"""
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile

W = "/kaggle/working"
OUT = f"{W}/refs"
ENV = "/tmp/env"
TOOLS = f"{ENV}/bin"
NCBI = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
        "seqs_for_alignment_pipelines.ucsc_ids/")
FNA = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
FA = f"{OUT}/GRCh38_no_alt_plus_hs38d1_analysis_set_no_chr.fa.gz"
CLINVAR_ARCHIVE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20260928.vcf.gz"
CLINVAR_CURRENT = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
build = {"started": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())}


def sh(cmd, check=True):
    print("+", cmd, flush=True)
    r = subprocess.run(["bash", "-c", "set -o pipefail; " + cmd], text=True, capture_output=True)
    if r.stdout.strip():
        print(r.stdout[-4000:], flush=True)
    if r.returncode and r.stderr.strip():
        print(r.stderr[-4000:], flush=True)
    if check and r.returncode:
        raise SystemExit(f"failed ({r.returncode}): {cmd}")
    return r


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_input(*names):
    """First regular file matching any name. Kaggle unpacks uploaded .gz files (sometimes into a
    folder of the same name), so callers list both the .gz and the plain name."""
    for name in names:
        hits = [h for h in glob.glob(f"/kaggle/input/**/{name}", recursive=True) if os.path.isfile(h)]
        if hits:
            return hits[0]
    raise SystemExit(f"input not found: {names}")


def gz_copy(src, dst):
    """Copy src to dst as bgzip, whether src is plain or gzip-compressed."""
    cat = "gzip -dc" if src.endswith(".gz") else "cat"
    sh(f"{cat} {src} | {TOOLS}/bgzip -@4 -c > {dst}")


def code_root():
    """The mva-code dataset holds repo/ (or repo.zip if Kaggle kept the archive)."""
    hits = glob.glob("/kaggle/input/**/tools/annotate_vcf.sh", recursive=True)
    if not hits:
        for z in glob.glob("/kaggle/input/**/repo.zip", recursive=True):
            zipfile.ZipFile(z).extractall("/tmp/code_zip")
        hits = glob.glob("/tmp/code_zip/**/tools/annotate_vcf.sh", recursive=True)
    if not hits:
        raise SystemExit("mva-code not found")
    src = os.path.dirname(os.path.dirname(hits[0]))
    shutil.copytree(src, "/tmp/code", dirs_exist_ok=True)  # writable copy
    return "/tmp/code"


os.makedirs(OUT, exist_ok=True)

# 1. Toolchain: bcftools, samtools and htslib from bioconda, matching the laptop's 1.24 if possible.
sh("mkdir -p /tmp/mm && curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /tmp/mm bin/micromamba")
r = sh(f"/tmp/mm/bin/micromamba create -y -q -p {ENV} -c conda-forge -c bioconda bcftools=1.24 samtools=1.24 htslib=1.24", check=False)
if r.returncode:
    sh(f"/tmp/mm/bin/micromamba create -y -q -p {ENV} -c conda-forge -c bioconda bcftools samtools htslib")
build["bcftools"] = sh(f"{TOOLS}/bcftools --version | head -1").stdout.strip()
build["samtools"] = sh(f"{TOOLS}/samtools --version | head -1").stdout.strip()

# 2. Reference: download, check NCBI's md5, rename contigs, bgzip, index.
sh(f"curl -fsSL --retry 5 -o /tmp/ref.fna.gz {NCBI}{FNA}")
sums = sh(f"curl -fsSL --retry 5 {NCBI}md5checksums.txt").stdout
want = [ln.split()[0] for ln in sums.splitlines() if ln.strip().endswith(FNA)]
got = md5("/tmp/ref.fna.gz")
build["reference"] = {"source": NCBI + FNA, "md5": got, "md5_matches_ncbi": bool(want) and want[0] == got}
if want and want[0] != got:
    raise SystemExit("reference md5 mismatch")
sh(f"gzip -dc /tmp/ref.fna.gz | sed 's/^>chr/>/' | {TOOLS}/bgzip -@4 -c > {FA}")
os.remove("/tmp/ref.fna.gz")
sh(f"{TOOLS}/samtools faidx {FA}")
build["reference"]["contigs"] = int(sh(f"wc -l < {FA}.fai").stdout.strip())

# 3. ClinVar: the archived 2026-09-28 release (the laptop's), else the current one.
cv = f"{OUT}/clinvar.vcf.gz"
r = sh(f"curl -fsSL --retry 3 -o {cv} {CLINVAR_ARCHIVE} && curl -fsSL --retry 3 -o {cv}.tbi {CLINVAR_ARCHIVE}.tbi", check=False)
src = CLINVAR_ARCHIVE
if r.returncode:
    src = CLINVAR_CURRENT
    sh(f"curl -fsSL --retry 3 -o {cv} {CLINVAR_CURRENT} && curl -fsSL --retry 3 -o {cv}.tbi {CLINVAR_CURRENT}.tbi")
sh(f"touch {cv}.tbi")  # the index must not look older than the VCF
r = sh(f"curl -fsSL --retry 3 {src}.md5", check=False)
build["clinvar"] = {"source": src, "fileDate": sh(f"{TOOLS}/bcftools view -h {cv} | grep -m1 '^##fileDate'").stdout.strip(),
                    "md5_matches_ncbi": (r.stdout.split() or [""])[0] == md5(cv) if not r.returncode else None}

# 4. gnomAD v4.1 slices (R10): bgzip + tabix; copy the BED, manifest, MANE GFF3 and panel genes.
g = f"{OUT}/gnomad_v4.1_joint_slices.vcf.gz"
gz_copy(find_input("gnomad_v4.1_joint_slices.vcf.gz", "gnomad_v4.1_joint_slices.vcf"), g)
sh(f"{TOOLS}/tabix -f -p vcf {g}")
gz_copy(find_input("GRCh38.116.mane_select.gff3.gz", "GRCh38.116.mane_select.gff3"), f"{OUT}/GRCh38.116.mane_select.gff3.gz")
for name in ("gnomad_v4.1_covered.bed", "manifest.json", "panel_genes_mva.txt"):
    shutil.copy(find_input(name), OUT)
build["gnomad"] = {"records": int(sh(f"{TOOLS}/bcftools index -n {g}").stdout.strip()),
                   "compressed_and_indexed": os.path.exists(g + ".tbi")}

# 5. Tests, natively: the pipeline unit tests and the synthetic end-to-end annotation check.
code = code_root()
env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTEST_DISABLE_PLUGIN_AUTOLOAD="1", MVA_NATIVE="1",
           TOOLS=TOOLS, MVA_GFF=f"{OUT}/GRCh38.116.mane_select.gff3.gz", GFF=f"{OUT}/GRCh38.116.mane_select.gff3.gz",
           MVA_FA=FA, FA=FA, CLINVAR=cv, LOW="", MVA_E2E_TMP="/tmp/mva-e2e")
r = subprocess.run([sys.executable, "-m", "pytest", "pipeline/tests", "-q", "-p", "no:cacheprovider"],
                   cwd=code, env=env, capture_output=True, text=True)
print(r.stdout[-3000:], r.stderr[-2000:], flush=True)
build["unit_tests"] = {"exit": r.returncode, "summary": (r.stdout.strip().splitlines() or [""])[-1]}
r = subprocess.run([sys.executable, "tools/test_annotate_e2e.py"], cwd=code, env=env, capture_output=True, text=True)
print(r.stdout[-6000:], r.stderr[-3000:], flush=True)
build["e2e"] = {"exit": r.returncode,
                "checks": [ln for ln in r.stdout.splitlines() if ln.startswith(("PASS ", "FAIL ", "SKIP"))]}

build["finished"] = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())
with open(f"{OUT}/build.json", "w") as fh:
    json.dump(build, fh, indent=1)
print(json.dumps(build, indent=1))
sh(f"du -sh {OUT}; ls -la {OUT}")
