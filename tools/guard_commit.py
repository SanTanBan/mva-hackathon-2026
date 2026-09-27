#!/usr/bin/env python3
"""Pre-commit guard that keeps gated and patient-derived data out of git.

.githooks/pre-commit runs this script on every commit. If it can't run, the
commit fails. Never bypass it with --no-verify. See CLAUDE.md for the rules.
"""
import re
import subprocess
import sys
from pathlib import Path

BLOCKED_SUFFIXES = (".vcf", ".vcf.gz", ".bcf", ".gvcf", ".gvcf.gz", ".fastq", ".fastq.gz",
                    ".fq", ".fq.gz", ".bam", ".bai", ".cram", ".crai", ".sam", ".tbi", ".csi")
BLOCKED_PREFIXES = ("private/", "research/", "data/", "work/", "submissions/")
MAX_BYTES = 5_000_000
# A VCF body line, with or without the chr prefix.
VCF_LINE = re.compile(rb"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)\t\d+\t", re.M)
# A predictions-style row: chrN, position, REF, ALT (comma- or tab-separated).
VARIANT_ROW = re.compile(rb"\bchr([1-9]|1[0-9]|2[0-2]|X|Y|M)[,\t]\d{4,}[,\t][ACGTN]+[,\t][ACGTN]+\b")
# Identifiers and names that must never be published. The list lives outside the repo.
PRIVATE_TERMS_FILE = Path(r"D:\mva-hackathon\private-terms.txt")


def staged_files():
    out = subprocess.run(["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR", "-z"],
                         capture_output=True, check=True).stdout
    return [p for p in out.decode("utf-8").split("\0") if p]


def staged_blob(path):
    return subprocess.run(["git", "show", f":{path}"], capture_output=True, check=True).stdout


def private_terms():
    if not PRIVATE_TERMS_FILE.is_file():
        return []
    terms = []
    for line in PRIVATE_TERMS_FILE.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            # Letters and digits bound the term but underscores don't, so the ID in
            # WGS_<ID>_<flowcell> still matches while a longer word like "Georgetown" doesn't.
            term = re.escape(line.encode("utf-8"))
            terms.append(re.compile(rb"(?<![A-Za-z0-9])" + term + rb"(?![A-Za-z0-9])", re.I))
    return terms


def main():
    problems = []
    terms = private_terms()
    for path in staged_files():
        low = path.lower().replace("\\", "/")
        if low.endswith(BLOCKED_SUFFIXES):
            problems.append(f"{path}: genomic data file type is never allowed in git")
            continue
        if low.startswith(BLOCKED_PREFIXES):
            problems.append(f"{path}: this folder is private and must stay out of git")
            continue
        blob = staged_blob(path)
        if len(blob) > MAX_BYTES:
            problems.append(f"{path}: {len(blob) / 1e6:.1f} MB exceeds the {MAX_BYTES / 1e6:.0f} MB limit")
            continue
        if b"\0" in blob[:8000]:
            continue  # binary file (image or PDF): the size check above is all that applies
        if VCF_LINE.search(blob):
            problems.append(f"{path}: contains VCF-style variant lines")
        if VARIANT_ROW.search(blob):
            problems.append(f"{path}: contains variant coordinate rows (chrom,pos,ref,alt); "
                            "cite public variants by ClinVar ID or HGVS instead")
        for term in terms:
            if term.search(blob):
                problems.append(f"{path}: contains a protected identifier or name "
                                f"(listed in {PRIVATE_TERMS_FILE})")
                break
    if problems:
        print("Commit blocked by tools/guard_commit.py (see CLAUDE.md):", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        print("Remove these from the commit. Do not bypass with --no-verify.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
