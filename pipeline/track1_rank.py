#!/usr/bin/env python3
"""Rank candidate genotypes and write the Track 1 submission CSV (package D02).

Written on 2026-10-02, before any access to the child's genome, so the scoring rule is fixed
before the variants are seen. Input is the annotated table from pipeline/annotate_table.py.
Candidates come from pipeline/comphet.py (homozygous and compound-het genotypes, phase from
reads). In BUB1B (the single_allele_genes argument), every rare heterozygous allele with
allele evidence of at least 0.3 also gets a single-variant row, alongside any pair it is part
of: the hypothesis that its partner was missed. MVA second alleles are often non-coding or
structural and can be missed by short-read calling.

Score of a candidate = gene prior x allele evidence x phase factor.
  gene prior     MVA genes: BUB1B 1.0, CEP57 0.6, TRIP13 0.6 (BUB1B is by far the most
                 common cause; locus dossier §7). Genes in --panel-genes (for example the
                 green genes of the PanelApp panels for the child's HPO terms) 0.1. Others 0.01.
  allele         ClinVar pathogenic or likely pathogenic 1.0; otherwise by impact:
                 HIGH 0.8, MODERATE 0.3, LOW 0.15, MODIFIER 0.05. ClinVar benign or likely
                 benign 0.01 whatever the impact. Conflicting or VUS entries fall back to impact.
                 A call that fails a VCF filter scores half.
  genotype       a compound het is the product of its two allele scores; a homozygote is
                 the allele score once (one allele, seen twice); a single allele whose partner
                 was not seen is the allele score x 0.1. That places it below a seen
                 moderate-impact partner of unknown phase (x 0.15) and above a seen low-impact
                 one (x 0.075).
  phase factor   trans 1.0, unknown 0.5, cis 0.05; 1.0 for homozygous and single rows.
Rows are sorted by score; ties go to the row with the stronger best allele, then by position.
A genotype listed under two overlapping genes is written once, with its higher score. The top
10 are written, with EPCR = score / (sum of the written rows' scores + 0.25). Normalising over
the written rows, not over every candidate, stops thousands of weak candidates (for example
intronic variants whose frequency is unknown without gnomAD) from diluting the strong ones.
The 0.25 is the probability that the cause is in none of the written rows, so no row can claim
certainty: the most any row can get is 0.8. EPCR is rounded to 3 decimals with a floor of
0.001, and a tied value is lowered by 0.001 so that the file's order is explicit.

Output is patient-derived and goes outside the repository; validate it with
tools/validate_track1_csv.py before staging.

Run: python -m pipeline.track1_rank TABLE.tsv OUT.csv --proband-id ID [--panel-genes FILE]
"""
import argparse
import csv

from pipeline import comphet
from pipeline.common import assert_outside_repo

FIELDS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
          "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]
GENE_PRIOR = {"BUB1B": 1.0, "CEP57": 0.6, "TRIP13": 0.6}
PANEL_PRIOR, OTHER_PRIOR = 0.1, 0.01
IMPACT_SCORE = {"HIGH": 0.8, "MODERATE": 0.3, "LOW": 0.15, "MODIFIER": 0.05}
PHASE = {"trans": 1.0, "unknown": 0.5, "cis": 0.05, "n/a": 1.0}
UNSEEN_PARTNER = 0.1
SINGLE_MIN = 0.3
FILTER_FAIL = 0.5
RESERVE = 0.25
MAX_ROWS = 10


def clinvar_class(clnsig):
    s = (clnsig or "").lower()
    if "conflicting" in s:
        return "conflicting"
    if "pathogenic" in s and "benign" not in s:
        return "plp"
    if "benign" in s and "pathogenic" not in s:
        return "blb"
    return ""


def allele_score(v):
    c = clinvar_class(v.get("clnsig"))
    if c == "plp":
        s = 1.0
    elif c == "blb":
        s = 0.01
    else:
        s = IMPACT_SCORE.get(v.get("impact"), 0.05)
    return s if v.get("filter", "PASS") in ("PASS", ".", "") else s * FILTER_FAIL


def gene_prior(gene, panel):
    return GENE_PRIOR.get(gene, PANEL_PRIOR if gene in panel else OTHER_PRIOR)


def describe(v):
    parts = [v["consequence"]]
    if v.get("aa_change"):
        parts.append(v["aa_change"])
    if clinvar_class(v.get("clnsig")):
        parts.append("ClinVar " + v["clnsig"].replace("_", " "))
    return " ".join(parts)


def score_candidates(rows, panel=(), max_af=0.01, single_allele_genes=("BUB1B",)):
    panel = set(panel)
    keep_all = tuple(set(GENE_PRIOR) | set(single_allele_genes))
    # A variant near two genes has one row per gene, so rows are looked up by gene and variant.
    by_key = {(v["gene"], comphet.label(v)): v for v in rows}
    scored = []
    for c in comphet.candidates(rows, max_af=max_af, keep_all_in=keep_all):
        a = by_key[(c["gene"], c["variant_1"])]
        prior = gene_prior(c["gene"], panel)
        if c["model"] == "homozygous":
            scored.append((prior * allele_score(a), a, None, f"{c['gene']} homozygous: {describe(a)}"))
        else:
            b = by_key[(c["gene"], c["variant_2"])]
            s = prior * allele_score(a) * allele_score(b) * PHASE[c["phase"]]
            scored.append((s, a, b, f"{c['gene']} compound het, phase {c['phase']}: {describe(a)}; {describe(b)}"))
    for v in rows:
        if (v["gene"] in single_allele_genes and v["gt"] in comphet.HET and allele_score(v) >= SINGLE_MIN
                and comphet.keep(v, max_af, {"HIGH", "MODERATE"}, set(keep_all))):
            s = gene_prior(v["gene"], panel) * allele_score(v) * UNSEEN_PARTNER
            scored.append((s, v, None, f"{v['gene']} one allele, second allele not seen in the VCF: {describe(v)}"))

    def order(x):
        s, a, b, _ = x
        best = max(allele_score(a), allele_score(b) if b is not None else 0.0)
        return (-s, -best, a["chrom"], int(a["pos"]), int(b["pos"]) if b is not None else 0)

    top, seen = [], set()
    for x in sorted(scored, key=order):
        key = tuple(sorted(comphet.label(v) for v in x[1:3] if v is not None))
        if key not in seen:
            seen.add(key)
            top.append(x)
        if len(top) == MAX_ROWS:
            break
    total = sum(s for s, *_ in top) + RESERVE
    out, last = [], None
    for s, a, b, note in top:
        epcr = max(0.001, round(s / total, 3))
        if last is not None and epcr >= last:
            epcr = max(0.001, round(last - 0.001, 3))
        out.append((epcr, a, b, note))
        last = epcr
    return out


def write_csv(path, proband_id, ranked):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for epcr, a, b, note in ranked:
            row = {"proband_id": proband_id, "chrom_1": a["chrom"], "pos_1": a["pos"], "ref_1": a["ref"],
                   "alt_1": a["alt"], "epcr": f"{epcr:.3f}", "finding_type": "primary", "notes": note}
            if b is not None:
                row.update({"chrom_2": b["chrom"], "pos_2": b["pos"], "ref_2": b["ref"], "alt_2": b["alt"]})
            w.writerow(row)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table")
    ap.add_argument("out")
    ap.add_argument("--proband-id", required=True)
    ap.add_argument("--panel-genes", help="file with one gene symbol per line")
    ap.add_argument("--max-af", type=float, default=0.01)
    args = ap.parse_args()
    out = assert_outside_repo(args.out)
    panel = []
    if args.panel_genes:
        with open(args.panel_genes, encoding="utf-8") as fh:
            panel = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    ranked = score_candidates(comphet.load(args.table), panel, args.max_af)
    write_csv(out, args.proband_id, ranked)
    print(f"{len(ranked)} rows -> {out}; top EPCR {ranked[0][0] if ranked else 'n/a'}")


if __name__ == "__main__":
    main()
