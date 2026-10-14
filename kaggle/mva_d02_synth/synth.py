"""Private Kaggle kernel santanubanerjee9/mva-d02-synth: dry run of kaggle/mva_d02/d02.py.

Public data only. It builds the synthetic VCF of tools/test_annotate_e2e.py (public ClinVar
BUB1B alleles with invented genotypes, plus two invented SNVs) and runs the real d02.py on it,
so the patient-data kernel is exercised end to end before any patient data exists on Kaggle.
"""
import glob
import os
import runpy
import shutil
import subprocess
import sys

ENV = "/tmp/env"
TOOLS = f"{ENV}/bin"


def sh(cmd):
    print("+", cmd, flush=True)
    subprocess.run(["bash", "-c", "set -o pipefail; " + cmd], check=True)


sh("mkdir -p /tmp/mm && curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /tmp/mm bin/micromamba")
r = subprocess.run(["bash", "-c", f"/tmp/mm/bin/micromamba create -y -q -p {ENV} -c conda-forge -c bioconda "
                    "bcftools=1.24 samtools=1.24 htslib=1.24"])
if r.returncode:
    sh(f"/tmp/mm/bin/micromamba create -y -q -p {ENV} -c conda-forge -c bioconda bcftools samtools htslib")

refs = os.path.dirname(glob.glob("/kaggle/input/**/refs/build.json", recursive=True)[0])
src = os.path.dirname(os.path.dirname(glob.glob("/kaggle/input/**/tools/annotate_vcf.sh", recursive=True)[0]))
shutil.copytree(src, "/tmp/code_synth", dirs_exist_ok=True)
os.environ.update(MVA_NATIVE="1", TOOLS=TOOLS, MVA_FA=f"{refs}/GRCh38_no_alt_plus_hs38d1_analysis_set_no_chr.fa.gz",
                  MVA_GFF=f"{refs}/GRCh38.116.mane_select.gff3.gz", MVA_E2E_TMP="/tmp/synth")
sys.path.insert(0, "/tmp/code_synth")
from pathlib import Path  # noqa: E402
from tools import test_annotate_e2e as e2e  # noqa: E402  (reads the env above at import)

os.makedirs("/tmp/synth", exist_ok=True)
e2e.build_vcf(Path("/tmp/synth/synth.vcf"))
sh(f"cd /tmp/synth && {TOOLS}/bgzip -f synth.vcf && {TOOLS}/tabix -p vcf synth.vcf.gz")

os.environ["MVA_D02_VCF"] = "/tmp/synth/synth.vcf.gz"
d02 = glob.glob("/kaggle/input/**/kaggle/mva_d02/d02.py", recursive=True)[0]
runpy.run_path(d02, run_name="__main__")
