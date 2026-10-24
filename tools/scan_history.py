#!/usr/bin/env python3
"""Pre-release scan: search every commit, blob and commit message for patient-derived strings and secrets.

Run before the repository is made public. Terms are read from files outside the repository:
  D:\\mva-hackathon\\private-terms.txt        identifiers also enforced by the commit guard
  D:\\mva-hackathon\\history-scan-terms.txt   extra terms for this scan only (for example allele names)
  D:\\mva-hackathon\\history-scan-allow.txt   accepted "term<TAB>path" pairs, e.g. a public worked example
Prints only the term number, the commit or path and a category, never a matched value. Also flags
VCF-style lines, chrom,pos,ref,alt rows and credential patterns. PDFs are scanned on their extracted text
too, because their streams are compressed. Exits 1 if anything not allowed is found.

Run: python tools/scan_history.py [--repo PATH]
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from guard_commit import VARIANT_ROW, VCF_LINE  # noqa: E402

ROOT = Path(r"D:\mva-hackathon")
SECRETS = {"hf token": rb"hf_[A-Za-z0-9]{30,}", "github token": rb"gh[pousr]_[A-Za-z0-9]{30,}",
           "kaggle key": rb'"key"\s*:\s*"[0-9a-f]{32}"', "private key": rb"-----BEGIN [A-Z ]*PRIVATE KEY-----"}


def read_lines(path):
    return [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines()
            if ln.strip() and not ln.startswith("#")] if path.is_file() else []


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repo", default=str(Path(__file__).resolve().parents[1]))
    a = ap.parse_args()

    def git(*args, data=None):
        return subprocess.run(["git", "-C", a.repo, *args], input=data, capture_output=True).stdout

    terms = read_lines(ROOT / "private-terms.txt") + read_lines(ROOT / "history-scan-terms.txt")
    allow = {tuple(ln.split("\t", 1)) for ln in read_lines(ROOT / "history-scan-allow.txt") if "\t" in ln}
    pats = [re.compile(r"(?<![A-Za-z0-9])" + re.escape(t) + r"(?![A-Za-z0-9])", re.I) for t in terms]
    found, allowed = [], []

    def report(i, where):
        (allowed if (terms[i], where) in allow else found).append(f"term #{i} in {where}")

    for entry in git("log", "--all", "--format=%H%x00%B%x01").decode("utf-8", "replace").split("\x01"):
        if "\x00" in entry:
            sha, msg = entry.strip().split("\x00", 1)
            for i, p in enumerate(pats):
                if p.search(msg):
                    report(i, f"commit message {sha[:8]}")
    paths = {}
    for line in git("rev-list", "--all", "--objects").decode("utf-8", "replace").splitlines():
        oid, _, path = line.partition(" ")
        if path:
            paths.setdefault(oid, path)
    kinds = git("cat-file", "--batch-check", data="\n".join(paths).encode()).decode().splitlines()
    blobs = [k.split()[0] for k in kinds if len(k.split()) > 1 and k.split()[1] == "blob"]
    for oid in blobs:
        data = git("cat-file", "-p", oid)
        if data[:5] == b"%PDF-":  # compressed streams: scan the extracted text as well
            try:
                import io
                import pypdf
                data += "\n".join(pg.extract_text() or "" for pg in pypdf.PdfReader(io.BytesIO(data)).pages).encode()
            except Exception as e:  # an unreadable PDF must not pass silently
                found.append(f"PDF text not extractable ({type(e).__name__}) in {paths[oid]}")
        if VCF_LINE.search(data) or VARIANT_ROW.search(data):
            found.append(f"VCF-style or variant row in {paths[oid]}")
        for name, pat in SECRETS.items():
            if re.search(pat, data):
                found.append(f"{name} pattern in {paths[oid]}")
        text = data.decode("utf-8", "replace")
        for i, p in enumerate(pats):
            if p.search(text):
                report(i, paths[oid])
    commits = len(git("rev-list", "--all").split())
    print(f"scanned {commits} commits, {len(blobs)} blobs, {len(terms)} terms")
    for f in sorted(set(allowed)):
        print("allowed:", f)
    for f in sorted(set(found)):
        print("FOUND:", f)
    print("clean" if not found else f"{len(set(found))} finding(s): do not publish")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
