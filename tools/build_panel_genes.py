#!/usr/bin/env python3
"""Build the --panel-genes file for pipeline/track1_rank.py from public PanelApp data.

Written on 2026-10-04, before any access to the child's genome, so the panel set is fixed
before the variants or the phenotype document are seen. The panels are the ten Genomics
England PanelApp panels on which BUB1B is green with biallelic inheritance (locus dossier §6):
the panel route by which an NHS clinician would reach the gene from MVA's phenotype classes
(microcephaly, childhood tumours, developmental disorder). They are chosen by panel, never
from the child's HPO terms.

A gene is kept if it is green (confidence level 3) on at least one of those panels and its mode
of inheritance there allows a biallelic or hemizygous genotype (BIALLELIC, BOTH, X-LINKED, or
blank). Monoallelic-only genes are dropped, because the ranker scores recessive genotypes.
PanelApp symbols can be out of date, so each gene is mapped to the symbol in the Ensembl 116
MANE Select GFF3 (the names bcftools csq writes) through its HGNC ID; the PanelApp symbol is
kept as a fallback when the HGNC ID is not in the GFF3.

Panel responses are cached in D:\\mva-hackathon\\refs\\panelapp, so a rerun with --offline
reproduces the file exactly. The output records the panel versions and the query date.

Run: python tools/build_panel_genes.py [--offline] [--out FILE]
"""
import argparse
import datetime
import gzip
import json
import re
import sys
import urllib.request
from pathlib import Path

API = "https://panelapp.genomicsengland.co.uk/api/v1/panels/{}/"
CACHE = Path(r"D:\mva-hackathon\refs\panelapp")
GFF = Path(r"D:\mva-hackathon\refs\ensembl\GRCh38.116.mane_select.gff3.gz")
OUT = CACHE / "panel_genes_mva.txt"
# Panels where BUB1B is green and biallelic (locus dossier §6, PanelApp API 2026-09-27).
PANELS = [259, 243, 290, 162, 179, 478, 484, 285, 1636, 1637]
RECESSIVE_MOI = ("BIALLELIC", "BOTH", "X-LINKED")


def fetch_panel(pid, offline=False):
    path = CACHE / f"panel_{pid}.json"
    if offline:
        return json.loads(path.read_text(encoding="utf-8"))
    req = urllib.request.Request(API.format(pid), headers={"Accept": "application/json", "User-Agent": "mva-hackathon-panel-builder/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = json.load(r)
    CACHE.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")
    return data


def hgnc_to_symbol(gff=GFF):
    acc = re.compile(r"Name=([^;]+);.*Acc:(HGNC:\d+)")
    out = {}
    with gzip.open(gff, "rt", encoding="utf-8") as fh:
        for line in fh:
            if "\tgene\t" in line:
                m = acc.search(line)
                if m:
                    out[m.group(2)] = m.group(1)
    return out


def allows_recessive(moi):
    moi = (moi or "").strip().upper()
    return not moi or moi.startswith(RECESSIVE_MOI)


def select_genes(panels, hgnc_map):
    """Return {symbol: [panel ids]} for green, recessive-compatible genes."""
    genes = {}
    for p in panels:
        for g in p.get("genes", []):
            if str(g.get("confidence_level")) != "3" or not allows_recessive(g.get("mode_of_inheritance")):
                continue
            gd = g.get("gene_data") or {}
            sym = hgnc_map.get(gd.get("hgnc_id"), gd.get("gene_symbol") or g["entity_name"])
            genes.setdefault(sym, set()).add(p["id"])
    return {s: sorted(ids) for s, ids in genes.items()}


def write(path, panels, genes, today):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# PanelApp green, recessive-compatible genes; built {today} by tools/build_panel_genes.py\n")
        for p in panels:
            fh.write(f"# panel {p['id']} v{p['version']}: {p['name']}\n")
        fh.write(f"# {len(genes)} genes\n")
        for sym in sorted(genes):
            fh.write(sym + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--offline", action="store_true", help="use cached panel responses only")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    panels = [fetch_panel(pid, args.offline) for pid in PANELS]
    genes = select_genes(panels, hgnc_to_symbol())
    today = datetime.date.today().isoformat()
    write(Path(args.out), panels, genes, today)
    missing = [p["id"] for p in panels if "BUB1B" not in genes or p["id"] not in genes["BUB1B"]]
    print(f"{len(genes)} genes from {len(panels)} panels -> {args.out}")
    if missing:
        print("warning: BUB1B is not green and recessive on panels", missing)
    return 0


if __name__ == "__main__":
    sys.exit(main())
