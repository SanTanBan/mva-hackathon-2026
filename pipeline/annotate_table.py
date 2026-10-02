#!/usr/bin/env python3
"""Turn the output of tools/annotate_vcf.sh into the table pipeline/comphet.py reads.

Input is a normalised VCF (one ALT per record) with INFO/BCSQ from `bcftools csq` on MANE
Select transcripts and the ClinVar fields added by `bcftools annotate`. The first sample is
used. One row is written per (variant, gene):
  - gene, consequence and impact come from the most severe BCSQ entry for that gene. BCSQ is
    Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change; entries that start
    with '@' point to another record (a joint codon change) and are skipped, and a leading '*'
    (a consequence on the other haplotype) is stripped.
  - Variants with no BCSQ entry are placed by gene interval from the same MANE Select GFF3:
    'intron_variant' inside a gene, 'upstream_gene_variant' or 'downstream_gene_variant'
    within --flank bases of it (strand-aware), all MODIFIER. The BUB1B upstream interval
    matters (locus dossier §4), so these rows are kept for comphet's keep-all genes.
  - pop_af is the highest of ClinVar's AF_EXAC, AF_TGP and AF_ESP, or empty. Empty means
    unknown, and comphet treats unknown as rare. A local gnomAD subset would replace this.
  - Chromosomes are written with the 'chr' prefix, as the Track 1 CSV requires.

Output is patient-derived when the input is, so it must be written outside the repository.

Run: python -m pipeline.annotate_table ANNOTATED.vcf.gz MANE.gff3.gz OUT.tsv [--flank 5000]
"""
import argparse
import bisect
import csv
from collections import defaultdict

from pipeline.common import add_chr, assert_outside_repo, open_text
from pipeline.comphet import COLUMNS, IMPACT_RANK
from pipeline.mane_gff import attrs

EXTRA = ("transcript", "aa_change", "dna_change", "clnsig", "clnrevstat", "clinvar_alleleid",
         "filter", "dp", "gq", "ad")

# bcftools csq consequence terms, mapped to VEP-style impact classes.
IMPACT = {
    "stop_gained": "HIGH", "frameshift": "HIGH", "splice_acceptor": "HIGH", "splice_donor": "HIGH",
    "start_lost": "HIGH", "stop_lost": "HIGH", "transcript_ablation": "HIGH",
    "missense": "MODERATE", "inframe_deletion": "MODERATE", "inframe_insertion": "MODERATE",
    "inframe_altering": "MODERATE", "protein_altering": "MODERATE",
    "splice_region": "LOW", "synonymous": "LOW", "stop_retained": "LOW", "start_retained": "LOW",
    "5_prime_utr": "MODIFIER", "3_prime_utr": "MODIFIER", "non_coding": "MODIFIER",
    "intron": "MODIFIER", "coding_sequence": "MODIFIER", "feature_elongation": "MODIFIER",
}


def impact_of(consequence):
    """Most severe impact among '&'-joined terms; unknown terms count as MODIFIER."""
    terms = consequence.lower().split("&")
    return min((IMPACT.get(t, "MODIFIER") for t in terms), key=lambda i: IMPACT_RANK[i])


def short(aa, keep=20):
    """A frameshift's amino-acid change lists the whole lost tail; keep both ends of each side."""
    sides = aa.split(">")
    return ">".join(x if len(x) <= 2 * keep + 3 else f"{x[:keep]}...{x[-keep:]}" for x in sides)


def parse_bcsq(value):
    """Return one entry per gene: the most severe consequence for that gene."""
    best = {}
    if not value or value == ".":
        return best
    for item in value.split(","):
        if item.startswith("@"):
            continue
        f = item.lstrip("*").split("|")
        if len(f) < 2 or not f[1]:
            continue
        f += [""] * (7 - len(f))
        entry = {"gene": f[1], "consequence": f[0], "impact": impact_of(f[0]),
                 "transcript": f[2], "aa_change": short(f[5]), "dna_change": f[6]}
        old = best.get(f[1])
        if old is None or IMPACT_RANK[entry["impact"]] < IMPACT_RANK[old["impact"]]:
            best[f[1]] = entry
    return best


class GeneIndex:
    """Gene intervals by chromosome from the MANE Select GFF3, for variants with no BCSQ."""

    def __init__(self, lines):
        by_chrom = defaultdict(list)
        for line in lines:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            a = attrs(cols[8])
            by_chrom[cols[0]].append((int(cols[3]), int(cols[4]), cols[6], a.get("Name", a.get("ID", ""))))
        self.genes = {c: sorted(v) for c, v in by_chrom.items()}
        self.starts = {c: [g[0] for g in v] for c, v in self.genes.items()}
        self.longest = max((g[1] - g[0] for v in self.genes.values() for g in v), default=0)

    def near(self, chrom, pos, flank):
        genes = self.genes.get(chrom, [])
        hi = bisect.bisect_right(self.starts.get(chrom, []), pos + flank)
        lo = bisect.bisect_left(self.starts.get(chrom, []), pos - flank - self.longest)
        out = {}
        for start, end, strand, name in genes[lo:hi]:
            if start <= pos <= end:
                out[name] = "intron_variant"
            elif start - flank <= pos < start:
                out.setdefault(name, "upstream_gene_variant" if strand == "+" else "downstream_gene_variant")
            elif end < pos <= end + flank:
                out.setdefault(name, "downstream_gene_variant" if strand == "+" else "upstream_gene_variant")
        return out


def info_dict(col):
    out = {}
    for kv in col.split(";"):
        k, _, v = kv.partition("=")
        out[k] = v
    return out


def pop_af(info):
    vals = []
    for key in ("AF_EXAC", "AF_TGP", "AF_ESP"):
        try:
            vals.append(float(info[key]))
        except (KeyError, ValueError):
            pass
    return f"{max(vals):g}" if vals else ""


def rows(vcf_lines, genes, flank=5000):
    for line in vcf_lines:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 10 or "," in c[4] or c[4] in ("*", ".") or c[4].startswith("<"):
            continue  # unsplit, spanning-deletion and symbolic alleles cannot be submitted
        info = info_dict(c[7])
        sample = dict(zip(c[8].split(":"), c[9].split(":")))
        base = {
            "chrom": add_chr(c[0]), "pos": c[1], "ref": c[3], "alt": c[4], "pop_af": pop_af(info),
            "gt": sample.get("GT", ""), "pgt": sample.get("PGT", "") if sample.get("PGT") != "." else "",
            "pid": sample.get("PID", "") if sample.get("PID") != "." else "",
            "clnsig": info.get("CLNSIG", ""), "clnrevstat": info.get("CLNREVSTAT", ""),
            "clinvar_alleleid": info.get("ALLELEID", ""), "filter": c[6],
            "dp": sample.get("DP", ""), "gq": sample.get("GQ", ""), "ad": sample.get("AD", ""),
        }
        per_gene = parse_bcsq(info.get("BCSQ"))
        for name, consequence in genes.near(c[0], int(c[1]), flank).items():
            per_gene.setdefault(name, {"gene": name, "consequence": consequence, "impact": "MODIFIER",
                                       "transcript": "", "aa_change": "", "dna_change": ""})
        if not per_gene:
            per_gene = {"": {"gene": "", "consequence": "intergenic_variant", "impact": "MODIFIER",
                             "transcript": "", "aa_change": "", "dna_change": ""}}
        for entry in per_gene.values():
            yield {**base, **entry}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf")
    ap.add_argument("gff")
    ap.add_argument("out")
    ap.add_argument("--flank", type=int, default=5000)
    args = ap.parse_args()
    out = assert_outside_repo(args.out)
    with open_text(args.gff) as fh:
        genes = GeneIndex(fh)
    n = 0
    with open_text(args.vcf) as fh, open(out, "w", encoding="utf-8", newline="") as oh:
        w = csv.DictWriter(oh, fieldnames=COLUMNS + EXTRA, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows(fh, genes, args.flank):
            w.writerow(r)
            n += 1
    print(f"{n} rows -> {out}")


if __name__ == "__main__":
    main()
