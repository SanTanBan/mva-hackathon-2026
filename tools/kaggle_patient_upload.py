#!/usr/bin/env python3
"""Upload the child's VCF and its index to the PRIVATE Kaggle dataset santanubanerjee9/mva-patient-vcf.

Approved by Santanu in a live chat on 2026-10-14 (CLAUDE.md, second exception): the VCF-first
files may be processed in private Kaggle datasets and kernels named mva-*, and are deleted by
2026-11-23. Run after D00 (tools/fetch_vcf_first.py) and D01 (tools/d01_intake.py).

It also writes proband_id.txt, the identifier the Track 1 scorer looks up. The candidates are the
VCF sample name, the VCF file stem and any proband ID that D01's pattern finds in the phenotype
document. The one named by most sources wins. Nothing identifying is printed: only whether the
sources agree.

The upload folder is D:\\mva-hackathon\\work\\kaggle_patient_vcf (outside the repository, hard
links where possible). TEMP and TMP point to D: because the kaggle CLI stages uploads there.

Run: python tools/kaggle_patient_upload.py [--dry-run]
"""
import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.common import assert_outside_repo  # noqa: E402

WORK = Path(r"D:\mva-hackathon\work")
FOLDER = WORK / "kaggle_patient_vcf"
DATASET = "santanubanerjee9/mva-patient-vcf"
KTMP = Path(r"D:\DevCache\tmp\claude\Child-with-rare-Disease\ktmp")


def kaggle(*args):
    KTMP.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TEMP=str(KTMP), TMP=str(KTMP))
    r = subprocess.run(["kaggle", *args], env=env, capture_output=True, text=True)
    out = "\n".join(ln for ln in (r.stdout + r.stderr).splitlines()
                    if "warn" not in ln.lower() and "%|" not in ln)
    return r.returncode, out


def vcf_sample(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    return []


def doc_ids():
    """Proband IDs found in the phenotype document by D01's own pattern."""
    try:
        from tools import d01_intake
    except Exception:
        return []
    inventory = json.loads((WORK / "d00_inventory.json").read_text(encoding="utf-8"))
    ids = []
    for f in inventory["files"]:
        if f.get("local") and f["kind"] == "phenotype-or-other":
            text = d01_intake.document_text(Path(f["local"]))
            ids += [t.rstrip("._-") for t in d01_intake.PROBAND.findall(text)]
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true", help="build the folder, don't upload")
    args = ap.parse_args()
    inventory = json.loads((WORK / "d00_inventory.json").read_text(encoding="utf-8"))
    local = [f for f in inventory["files"] if f.get("local")]
    vcfs = [Path(f["local"]) for f in local if f["kind"] == "vcf"]
    idx = [Path(f["local"]) for f in local if f["kind"] == "vcf-index"]
    if len(vcfs) != 1:
        raise SystemExit(f"expected one VCF in the D00 inventory, found {len(vcfs)}")
    if not all(f.get("local_sha256_ok") for f in local):
        raise SystemExit("D00 checksums are not all OK; rerun tools/fetch_vcf_first.py first")

    folder = assert_outside_repo(FOLDER)
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    for src in vcfs + idx:
        dst = folder / src.name
        try:
            os.link(src, dst)
        except OSError:
            shutil.copy2(src, dst)

    samples = vcf_sample(vcfs[0])
    stem = vcfs[0].name.split(".")[0]
    votes = Counter([s for s in samples[:1]] + [stem] + doc_ids())
    proband, n = votes.most_common(1)[0]
    (folder / "proband_id.txt").write_text(proband + "\n", encoding="utf-8")
    print(f"proband_id: chosen by {n} of {sum(votes.values())} sources "
          f"(VCF sample, file stem, phenotype document); all agree: {len(votes) == 1}")
    (folder / "dataset-metadata.json").write_text(json.dumps(
        {"title": "mva-patient-vcf", "id": DATASET, "licenses": [{"name": "other"}]}), encoding="utf-8")
    print(f"upload folder: {folder} ({sum(p.stat().st_size for p in folder.iterdir()) / 1e6:.0f} MB)")
    if args.dry_run:
        return 0

    rc, out = kaggle("datasets", "status", DATASET)
    if rc == 0 and "ready" in out:
        rc, out = kaggle("datasets", "version", "-p", str(folder), "-m", "VCF-first files (private)")
    else:
        rc, out = kaggle("datasets", "create", "-p", str(folder))  # private unless --public is given
    print(out.strip()[-800:])
    if rc:
        raise SystemExit(f"kaggle upload failed ({rc})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
