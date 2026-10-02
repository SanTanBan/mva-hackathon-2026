#!/usr/bin/env python3
"""Cut the Ensembl GFF3 down to MANE Select transcripts, for `bcftools csq`.

Consequences in the Track 1 report are pinned to MANE Select (methods B3), and readthrough
transcripts must not give a gene a spurious consequence (locus dossier §1). `bcftools csq`
annotates every transcript it is given, so the simplest way to pin it is to give it only the
MANE Select transcripts. This keeps:
  - transcript lines whose 'tag' attribute contains MANE_Select (not MANE_Plus_Clinical);
  - the gene lines that are parents of a kept transcript;
  - exon, CDS, five_prime_UTR and three_prime_UTR lines whose parent is a kept transcript;
  - the ##gff-version and ##sequence-region headers.
Everything else (other transcripts, ncRNA genes, pseudogenes) is dropped.

The input and output are public reference data, so the output goes to
D:\\mva-hackathon\\refs\\ensembl next to the input (outside the repository).

Run: python -m pipeline.mane_gff IN.gff3.gz OUT.gff3.gz
"""
import argparse
import gzip
import json
import re

from pipeline.common import assert_outside_repo, open_text

CHILD_TYPES = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR"}
ATTR = re.compile(r"(?:^|;)([^=;]+)=([^;]*)")


def attrs(col9):
    return dict(ATTR.findall(col9.rstrip("\n")))


def is_mane_select(a):
    return "MANE_Select" in a.get("tag", "").split(",")


def mane_ids(lines):
    """First pass: the IDs of MANE Select transcripts and of their parent genes."""
    transcripts, genes = set(), set()
    for line in lines:
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        a = attrs(cols[8])
        if "Parent" in a and a.get("ID", "").startswith("transcript:") and is_mane_select(a):
            transcripts.add(a["ID"])
            genes.add(a["Parent"])
    return transcripts, genes


def keep_lines(lines, transcripts, genes):
    """Second pass: yield the kept lines in input order."""
    for line in lines:
        if line.startswith("##gff-version") or line.startswith("##sequence-region"):
            yield line
            continue
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        a = attrs(cols[8])
        if a.get("ID") in genes or a.get("ID") in transcripts:
            yield line
        elif cols[2] in CHILD_TYPES and a.get("Parent") in transcripts:
            yield line


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gff_in")
    ap.add_argument("gff_out")
    args = ap.parse_args()
    out = assert_outside_repo(args.gff_out)
    with open_text(args.gff_in) as fh:
        transcripts, genes = mane_ids(fh)
    n = 0
    with open_text(args.gff_in) as fh, gzip.open(out, "wt", encoding="utf-8", newline="") as oh:
        for line in keep_lines(fh, transcripts, genes):
            oh.write(line)
            n += 1
    print(json.dumps({"genes": len(genes), "transcripts": len(transcripts), "lines": n}))


if __name__ == "__main__":
    main()
