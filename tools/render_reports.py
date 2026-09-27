#!/usr/bin/env python3
"""Render a report as a redacted public version and a full submission version.

The source Markdown marks patient-level content in one of two ways:
    <!-- private:start --> ... <!-- private:end -->    blocks
    [[private: ...]]                                     inline spans

- The public version replaces each marked part with a "withheld" note. It is written to
  reports/<track>/ (or --public-dir).
- The full version keeps the content, with the markers removed. It is written to
  D:\\mva-hackathon\\submissions\\<track>\\<date>\\ and never inside the repository.

Sources that contain patient data live outside the repository, for example
D:\\mva-hackathon\\work\\reports\\<track>\\<name>.src.md. Both versions are written as
Markdown and, when Microsoft Edge is installed, as PDF (Edge's headless print, so nothing
has to be downloaded).

The public version is refused if it still contains anything the commit guard would block
(VCF-style lines, variant rows, protected identifiers) or a chrom-pos-ref-alt token.
See CLAUDE.md.

Run: python tools/render_reports.py SOURCE --track track1 --name methods_report
"""
import argparse
import html
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import date
from pathlib import Path

import markdown

sys.path.insert(0, str(Path(__file__).parent))
from guard_commit import PRIVATE_TERMS_FILE, VARIANT_ROW, VCF_LINE  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
SUBMISSIONS = Path(os.environ.get("MVA_SUBMISSIONS_DIR", r"D:\mva-hackathon\submissions"))
TERMS_FILE = Path(os.environ.get("MVA_PRIVATE_TERMS", str(PRIVATE_TERMS_FILE)))
WITHHELD = "*[Withheld from the public version: patient-level detail. It appears in the submitted report.]*"

BLOCK = re.compile(r"<!--\s*private:start\s*-->(.*?)<!--\s*private:end\s*-->", re.S)
INLINE = re.compile(r"\[\[private:\s*(.*?)\]\]", re.S)
# chrom-pos-ref-alt in any common spelling: 15-40000000-A-G, chr15:40000000:A>G, chr15_40000000_A_G
VARIANT_TOKEN = re.compile(
    r"\b(?:chr|Chr|CHR)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)[:_-]\d{5,}[:_-]?[ACGTN]+(?:>|[:_/-])[ACGTN]+\b")
EDGE_PATHS = (r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
              r"C:\Program Files\Microsoft\Edge\Application\msedge.exe")
CSS = """body{font-family:Segoe UI,Arial,sans-serif;max-width:48em;margin:2em auto;line-height:1.45;color:#111}
table{border-collapse:collapse;margin:1em 0;font-size:.9em}th,td{border:1px solid #bbb;padding:.3em .5em;vertical-align:top}
th{background:#f0f0f0}code{font-size:.9em}blockquote{border-left:4px solid #ccc;margin-left:0;padding-left:1em;color:#333}"""


def split_versions(text):
    """Return (public, full) Markdown from a source with private markers."""
    public = INLINE.sub(WITHHELD, BLOCK.sub("\n" + WITHHELD + "\n", text))
    full = INLINE.sub(r"\1", BLOCK.sub(r"\1", text))
    return public, full


def private_terms():
    if not TERMS_FILE.is_file():
        return []
    lines = TERMS_FILE.read_text(encoding="utf-8").splitlines()
    return [re.compile(r"(?<![A-Za-z0-9])" + re.escape(t.strip()) + r"(?![A-Za-z0-9])", re.I)
            for t in lines if t.strip() and not t.strip().startswith("#")]


def public_problems(text):
    """Everything that makes a public version unsafe to write."""
    problems = []
    blob = text.encode("utf-8")
    if VCF_LINE.search(blob):
        problems.append("contains a VCF-style variant line")
    if VARIANT_ROW.search(blob):
        problems.append("contains a variant coordinate row (chrom,pos,ref,alt)")
    if VARIANT_TOKEN.search(text):
        problems.append("contains a chrom-pos-ref-alt variant token")
    if any(t.search(text) for t in private_terms()):
        problems.append(f"contains a protected identifier or name (listed in {TERMS_FILE})")
    if "private:start" in text or "[[private:" in text:
        problems.append("contains an unbalanced private marker")
    return problems


def find_edge():
    found = shutil.which("msedge")
    if found:
        return found
    return next((p for p in EDGE_PATHS if Path(p).is_file()), None)


def write_pdf(md_text, pdf_path, title):
    """Markdown -> HTML -> PDF with headless Edge. Returns False if Edge is unavailable."""
    edge = find_edge()
    if not edge:
        return False
    body = markdown.markdown(md_text, extensions=["tables", "fenced_code", "sane_lists"])
    page = (f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title>"
            f"<style>{CSS}</style></head><body>{body}</body></html>")
    # Scratch files sit next to the output, so a full version's content never lands on C:.
    work = Path(tempfile.mkdtemp(prefix=".render-", dir=pdf_path.parent))
    try:
        src = work / "page.html"
        src.write_text(page, encoding="utf-8")
        subprocess.run([edge, "--headless", "--disable-gpu", "--no-first-run", "--no-pdf-header-footer",
                        f"--user-data-dir={work / 'profile'}", f"--print-to-pdf={pdf_path}", src.as_uri()],
                       capture_output=True, timeout=180, check=False)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return pdf_path.is_file() and pdf_path.read_bytes()[:5] == b"%PDF-"


def render(source, track, name, day=None, public_dir=None, pdf=True):
    text = Path(source).read_text(encoding="utf-8")
    public, full = split_versions(text)
    problems = public_problems(public)
    if problems:
        raise SystemExit("Public version refused:\n  - " + "\n  - ".join(problems))
    full_dir = (SUBMISSIONS / track / (day or date.today().isoformat())).resolve()
    if full_dir == REPO or REPO in full_dir.parents:
        raise SystemExit(f"Full version must be written outside the repository, not {full_dir}")
    public_dir = Path(public_dir) if public_dir else REPO / "reports" / track
    public_dir.mkdir(parents=True, exist_ok=True)
    full_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for folder, content, label in ((public_dir, public, "public"), (full_dir, full, "full")):
        md_path = folder / f"{name}.md"
        md_path.write_text(content, encoding="utf-8", newline="\n")
        outputs.append(md_path)
        if pdf:
            pdf_path = folder / f"{name}.pdf"
            if write_pdf(content, pdf_path, f"{name} ({label})"):
                outputs.append(pdf_path)
            else:
                print(f"warning: no PDF for the {label} version (Edge not found or printing failed)",
                      file=sys.stderr)
    return outputs


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("source", help="Markdown source with private markers (keep patient data outside the repo)")
    ap.add_argument("--track", required=True, choices=["track1", "track2"])
    ap.add_argument("--name", required=True, help="output file stem, e.g. methods_report")
    ap.add_argument("--date", help="submission folder date (default: today)")
    ap.add_argument("--public-dir", help="override the public output folder (default: reports/<track>)")
    ap.add_argument("--no-pdf", action="store_true")
    args = ap.parse_args()
    for path in render(args.source, args.track, args.name, args.date, args.public_dir, not args.no_pdf):
        print(path)


if __name__ == "__main__":
    main()
