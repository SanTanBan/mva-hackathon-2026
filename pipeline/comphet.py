#!/usr/bin/env python3
"""Rare-recessive filter and compound-heterozygous pairing with read-backed phase.

Input is an annotated, tab-separated table with one row per (normalised) variant and these
columns: chrom, pos, ref, alt, gene, consequence, impact, pop_af, gt, pgt, pid.
Annotation (for example VEP pinned to MANE Select) runs upstream; this script only
reasons about inheritance.

- Keeps rare variants (pop_af <= --max-af; a missing frequency counts as rare) whose impact
  is in --impacts. In the genes given by --keep-all-in (default BUB1B), every rare variant
  is kept, because MVA1 second alleles are often non-coding.
- Reports homozygous candidates, and every pair of heterozygous variants in the same gene.
- Phase: two heterozygous variants with the same read-backed phase set (PID) are 'cis' when
  their phased genotypes (PGT) match and 'trans' when they differ. Otherwise phase is
  'unknown'. With no parental data, 'unknown' must be scored as such (PM3 at reduced strength).

Run: python -m pipeline.comphet ANNOTATED.tsv OUT.tsv   (output outside the repo)
"""
import argparse
import csv
from collections import defaultdict
from itertools import combinations

from pipeline.common import assert_outside_repo

COLUMNS = ("chrom", "pos", "ref", "alt", "gene", "consequence", "impact", "pop_af", "gt", "pgt", "pid")
HET = {"0/1", "1/0", "0|1", "1|0"}
HOM = {"1/1", "1|1"}
PHASE_RANK = {"trans": 0, "unknown": 1, "cis": 2}
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def load(path):
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    missing = [c for c in COLUMNS if rows and c not in rows[0]]
    if missing:
        raise SystemExit(f"missing columns: {missing}")
    return rows


def allele_freq(v):
    try:
        return float(v.get("pop_af") or 0.0)
    except ValueError:
        return 0.0


def keep(v, max_af, impacts, keep_all_in):
    if allele_freq(v) > max_af:
        return False
    return v["gene"] in keep_all_in or v["impact"] in impacts


def phase(a, b):
    pa, pb = a.get("pgt") or "", b.get("pgt") or ""
    if a.get("pid") and a.get("pid") == b.get("pid") and "|" in pa and "|" in pb:
        return "cis" if pa == pb else "trans"
    return "unknown"


def label(v):
    return f"{v['chrom']}:{v['pos']}:{v['ref']}>{v['alt']}"


def candidates(rows, max_af=0.01, impacts=("HIGH", "MODERATE"), keep_all_in=("BUB1B",)):
    kept = [v for v in rows if keep(v, max_af, set(impacts), set(keep_all_in))]
    out = []
    for v in kept:
        if v["gt"] in HOM:
            out.append({"gene": v["gene"], "model": "homozygous", "phase": "n/a",
                        "variant_1": label(v), "variant_2": "", "worst_impact": v["impact"]})
    by_gene = defaultdict(list)
    for v in kept:
        if v["gt"] in HET:
            by_gene[v["gene"]].append(v)
    for gene, vs in by_gene.items():
        for a, b in combinations(sorted(vs, key=lambda x: int(x["pos"])), 2):
            worst = min((a["impact"], b["impact"]), key=lambda i: IMPACT_RANK.get(i, 9))
            out.append({"gene": gene, "model": "compound_het", "phase": phase(a, b),
                        "variant_1": label(a), "variant_2": label(b), "worst_impact": worst})
    out.sort(key=lambda r: (PHASE_RANK.get(r["phase"], 1), IMPACT_RANK.get(r["worst_impact"], 9), r["gene"]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("table")
    ap.add_argument("out")
    ap.add_argument("--max-af", type=float, default=0.01)
    ap.add_argument("--impacts", default="HIGH,MODERATE")
    ap.add_argument("--keep-all-in", default="BUB1B")
    args = ap.parse_args()
    assert_outside_repo(args.out)
    rows = candidates(load(args.table), args.max_af, args.impacts.split(","), args.keep_all_in.split(","))
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "model", "phase", "variant_1", "variant_2", "worst_impact"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} candidate genotypes written")


if __name__ == "__main__":
    main()
