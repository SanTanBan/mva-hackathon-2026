#!/usr/bin/env python3
"""Package D01: data intake and integrity, run locally after D00 (tools/fetch_vcf_first.py).

Reads the downloaded VCF and phenotype document under D:\\mva-hackathon\\data and writes, outside the
repository:
- D:\\mva-hackathon\\work\\d01_intake.json: the VCF header summary (format, reference, contig naming,
  sample count, callers, INFO/FORMAT fields, record and PASS counts) and the checksum status from D00.
  No genotypes are written.
- D:\\mva-hackathon\\work\\hpo_terms.tsv: the HPO IDs found in the phenotype document, with a label
  when the document gives one on the same line.
It also appends the VCF sample names and any proband ID found in the phenotype document to
D:\\mva-hackathon\\private-terms.txt, so the pre-commit guard blocks them from git.

Run: python tools/d01_intake.py
"""
import argparse
import json
import re
import subprocess
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.common import assert_outside_repo, open_text  # noqa: E402
from pipeline.vcf_qc import check_header, read_header  # noqa: E402

DATA = Path(r"D:\mva-hackathon\data")
WORK = Path(r"D:\mva-hackathon\work")
PRIVATE_TERMS = Path(r"D:\mva-hackathon\private-terms.txt")
HPO_LINE = re.compile(r"(HP:\d{7})\s*[-:|\t]*\s*([A-Za-z][^\t\r\n|;]{2,80})?")
# Only an explicit "proband ID" label counts, and the ID must contain a digit, so ordinary words never
# become protected terms (the guard would then block every commit that uses them).
PROBAND = re.compile(r"proband[\s_-]*id\s*[:=#]?\s*([A-Za-z0-9_.-]*\d[A-Za-z0-9_.-]*)", re.I)
TEXT_SUFFIXES = {".txt", ".tsv", ".csv", ".json", ".md", ".yaml", ".yml", ".obo", ".phenopacket", ".xml"}


def document_text(path):
    """Plain text from a phenotype document: text formats, .docx, .xlsx or .pdf (via pdftotext)."""
    suffix = path.suffix.lower()
    if suffix in TEXT_SUFFIXES or suffix == "":
        return path.read_text(encoding="utf-8", errors="replace")
    if suffix == ".docx":
        with zipfile.ZipFile(path) as z:
            xml = z.read("word/document.xml").decode("utf-8", errors="replace")
        xml = re.sub(r"</w:p>", "\n", xml)
        return re.sub(r"<[^>]+>", "", xml)
    if suffix == ".xlsx":
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True)
        return "\n".join("\t".join("" if v is None else str(v) for v in row)
                         for ws in wb.worksheets for row in ws.iter_rows(values_only=True))
    if suffix == ".pdf":
        for cmd in (["pdftotext", "-layout", str(path), "-"],):
            try:
                return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
            except (OSError, subprocess.CalledProcessError):
                pass
        raise SystemExit(f"{path.name}: PDF text extraction needs pdftotext (poppler); install it on D: first")
    raise SystemExit(f"{path.name}: unsupported phenotype document type {suffix}")


def hpo_terms(text):
    terms = {}
    for m in HPO_LINE.finditer(text):
        label = (m.group(2) or "").strip()
        if m.group(1) not in terms or (label and not terms[m.group(1)]):
            terms[m.group(1)] = label
    return terms


def append_private_terms(terms, path=PRIVATE_TERMS):
    existing = set()
    if path.exists():
        existing = {l.strip() for l in path.read_text(encoding="utf-8").splitlines() if l.strip() and not l.startswith("#")}
    new = [t for t in dict.fromkeys(terms) if t and any(ch.isdigit() for ch in t) and t not in existing]
    if new:
        with open(path, "a", encoding="utf-8") as fh:
            fh.write("".join(f"{t}\n" for t in new))
    return new


def vcf_summary(path):
    with open_text(path) as fh:
        meta, _, info = read_header(fh)
        records = passed = 0
        for line in fh:
            records += 1
            cols = line.split("\t", 7)
            if len(cols) > 6 and cols[6] in ("PASS", "."):
                passed += 1
    return {
        "path": path.name, "fileformat": info["fileformat"], "reference": info["reference"],
        "n_contigs": len(info["contigs"]), "chr_prefix": any(c.startswith("chr") for c in info["contigs"]),
        "n_samples": len(info["samples"]), "header_problems": check_header(info),
        "sources": [m for m in meta if m.startswith(("##source", "##Sentieon", "##GATKCommandLine", "##DRAGEN"))][:8],
        "info_fields": sorted(info["info_number"]), "records": records, "pass_or_missing_filter": passed,
    }, info["samples"]


def main():
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    out_json = assert_outside_repo(WORK / "d01_intake.json")
    inventory = json.loads((WORK / "d00_inventory.json").read_text(encoding="utf-8"))
    local = {f["path"]: f for f in inventory["files"] if f.get("local")}
    vcfs = [Path(f["local"]) for f in local.values() if f["kind"] == "vcf"]
    docs = [Path(f["local"]) for f in local.values() if f["kind"] == "phenotype-or-other"]
    summary = {"checksums_ok": all(f.get("local_sha256_ok") for f in local.values()), "vcf": [], "phenotype": []}
    private = []
    for v in vcfs:
        s, samples = vcf_summary(v)
        summary["vcf"].append(s)
        private += samples
    all_terms = {}
    for d in docs:
        text = document_text(d)
        terms = hpo_terms(text)
        all_terms.update({k: v for k, v in terms.items() if k not in all_terms or v})
        private += [t.rstrip("._-") for t in PROBAND.findall(text)]
        summary["phenotype"].append({"path": d.name, "chars": len(text), "hpo_terms": len(terms)})
    with open(assert_outside_repo(WORK / "hpo_terms.tsv"), "w", encoding="utf-8") as fh:
        fh.write("hpo_id\tlabel\n" + "".join(f"{k}\t{v}\n" for k, v in sorted(all_terms.items())))
    summary["private_terms_added"] = len(append_private_terms(private))
    out_json.write_text(json.dumps(summary, indent=1), encoding="utf-8")
    print(f"VCFs {len(vcfs)}, phenotype documents {len(docs)}, HPO terms {len(all_terms)}, "
          f"private terms added {summary['private_terms_added']}; summary in {out_json}")


if __name__ == "__main__":
    main()
