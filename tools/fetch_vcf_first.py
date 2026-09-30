#!/usr/bin/env python3
"""Package D00: download the VCF-first files from the gated dataset and verify their checksums.

Approved by Santanu in a live chat on 2026-09-30 (laptop, VCF-first route): the VCF, its index and
the phenotype document go to D:\\mva-hackathon\\data. FASTQ pairs are fetched only with --fastq-pair,
one pair at a time, for an overnight job that deletes them afterwards. Nothing is written inside the
repository or on C: (HF_HOME points to D:\\DevCache\\huggingface).

The inventory of every dataset file (path, size, SHA-256 from the Hub's LFS metadata) and the result
of each local checksum go to D:\\mva-hackathon\\work\\d00_inventory.json.

Run: python tools/fetch_vcf_first.py                 (VCF, index, phenotype document)
     python tools/fetch_vcf_first.py --fastq-pair 1  (first FASTQ pair only, for an overnight job)
"""
import argparse
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = "SageBio/mva-hackathon-2026-data"
DATA = Path(r"D:\mva-hackathon\data")
WORK = Path(r"D:\mva-hackathon\work")
os.environ.setdefault("HF_HOME", r"D:\DevCache\huggingface")

FASTQ = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
VCF = re.compile(r"\.vcf(\.gz)?$", re.I)
INDEX = re.compile(r"\.(tbi|csi)$", re.I)
ALIGN = re.compile(r"\.(bam|bai|cram|crai)$", re.I)


def classify(path):
    name = path.rsplit("/", 1)[-1]
    if name.startswith(".") or name.lower() in ("readme.md", "license", "license.md"):
        return "repo-metadata"
    if FASTQ.search(name):
        return "fastq"
    if VCF.search(name):
        return "vcf"
    if INDEX.search(name):
        return "vcf-index"
    if ALIGN.search(name):
        return "alignment"
    return "phenotype-or-other"


def fastq_pairs(paths):
    """Group FASTQ files into read pairs by replacing the R1/R2 (or _1/_2) mate tag."""
    groups = {}
    for p in sorted(paths):
        key = re.sub(r"([._-])(R?)([12])([._-])", r"\1\2#\4", p, count=1)
        groups.setdefault(key, []).append(p)
    return [sorted(v) for _, v in sorted(groups.items())]


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while b := fh.read(chunk):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fastq-pair", type=int, help="also fetch this FASTQ pair (1-based), for an overnight job")
    args = ap.parse_args()
    from huggingface_hub import HfApi, hf_hub_download

    api = HfApi()
    try:
        api.auth_check(REPO, repo_type="dataset")
    except Exception as e:
        sys.exit(f"dataset access not granted ({type(e).__name__}); nothing downloaded")
    files = []
    for f in api.list_repo_tree(REPO, repo_type="dataset", recursive=True):
        if getattr(f, "size", None) is None:
            continue
        lfs = getattr(f, "lfs", None)
        files.append({"path": f.path, "bytes": f.size, "kind": classify(f.path),
                      "sha256": getattr(lfs, "sha256", None) if lfs else None})
    wanted = [f for f in files if f["kind"] in ("vcf", "vcf-index", "phenotype-or-other")]
    pairs = fastq_pairs([f["path"] for f in files if f["kind"] == "fastq"])
    if args.fastq_pair:
        pair = set(pairs[args.fastq_pair - 1])
        wanted += [f for f in files if f["path"] in pair]
    DATA.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    for f in wanted:
        local = Path(hf_hub_download(REPO, f["path"], repo_type="dataset", local_dir=str(DATA)))
        f["local"] = str(local)
        f["local_sha256_ok"] = (sha256(local) == f["sha256"]) if f["sha256"] else None
        print(f"{f['kind']:<20} {f['bytes'] / 1e9:7.2f} GB  checksum {'ok' if f['local_sha256_ok'] else f['local_sha256_ok']}")
    report = {"checked_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "repo": REPO,
              "files": files, "fastq_pairs": pairs, "downloaded": [f["path"] for f in wanted]}
    (WORK / "d00_inventory.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    bad = [f["path"] for f in wanted if f.get("local_sha256_ok") is False]
    if bad:
        sys.exit(f"checksum mismatch: {bad}")
    print(f"inventory written to {WORK / 'd00_inventory.json'}")


if __name__ == "__main__":
    main()
