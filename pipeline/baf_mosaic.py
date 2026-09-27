#!/usr/bin/env python3
"""Per-chromosome B-allele-fraction (BAF) test for mosaic aneuploidy, from VCF allele depths.

A fraction f of trisomic cells splits heterozygous-site allele fractions into two bands at
c - delta and c + delta, where c is the genome-wide centre (it absorbs reference bias) and
delta = f / (2 (2 + f)). For monosomy with cell fraction m, delta = m / (2 (2 - m)).
For each chromosome the script fits a symmetric two-binomial mixture by maximum
likelihood over delta. It reports delta, the likelihood-ratio statistic against delta = 0,
the implied trisomy and monosomy fractions, and a robust z-score of delta against the
autosomes.

Caveats, which must be stated with every result:
- Variegated mosaicism puts different chromosomes in different cells, often each at a low
  fraction. It may therefore fall below this test's limit even when the karyotype shows it.
- Copy-neutral loss of heterozygosity, mapping artefacts and local copy-number changes can
  also shift BAF. Treat a signal as a lead to check against read depth, not a diagnosis.
- simulate() gives the detection limit for the observed depth and site count.

Only biallelic SNVs with a heterozygous GT, a PASS filter and min_depth <= DP <= max_depth
are used. Memory use is small: the fit works on counts of (alt, depth) pairs.

Run: python -m pipeline.baf_mosaic NORMALISED.vcf OUT.tsv   (output outside the repo)
"""
import argparse
import csv
import math
import random
from collections import Counter, defaultdict
from statistics import median

from pipeline.common import assert_outside_repo, is_autosome, open_text

HET = {"0/1", "1/0", "0|1", "1|0"}


def collect(vcf_path, min_depth=20, max_depth=200):
    """Return {chrom: Counter((alt, depth) -> n)} for usable heterozygous SNVs."""
    by_chrom = defaultdict(Counter)
    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 10 or len(f[3]) != 1 or len(f[4]) != 1 or f[6] not in ("PASS", "."):
                continue
            s = dict(zip(f[8].split(":"), f[9].split(":")))
            if s.get("GT") not in HET or "AD" not in s:
                continue
            ad = s["AD"].split(",")
            if len(ad) != 2 or "." in ad:
                continue
            ref, alt = int(ad[0]), int(ad[1])
            depth = ref + alt
            if min_depth <= depth <= max_depth:
                by_chrom[f[0]][(alt, depth)] += 1
    return by_chrom


def loglik(counts, centre, delta):
    lo, hi = max(centre - delta, 1e-6), min(centre + delta, 1 - 1e-6)
    total = 0.0
    for (a, d), n in counts.items():
        l1 = a * math.log(lo) + (d - a) * math.log(1 - lo)
        l2 = a * math.log(hi) + (d - a) * math.log(1 - hi)
        m = max(l1, l2)
        total += n * (m + math.log(0.5 * math.exp(l1 - m) + 0.5 * math.exp(l2 - m)))
    return total


def fit_delta(counts, centre):
    """Coarse-to-fine grid search for the ML delta in [0, 0.25]."""
    grid = [i * 0.01 for i in range(26)]
    best = max(grid, key=lambda d: loglik(counts, centre, d))
    fine = [max(0.0, best - 0.01 + i * 0.0005) for i in range(41)]
    best = max(fine, key=lambda d: loglik(counts, centre, d))
    return best, 2 * (loglik(counts, centre, best) - loglik(counts, centre, 0.0))


def trisomy_fraction(delta):
    return 4 * delta / (1 - 2 * delta)


def monosomy_fraction(delta):
    return 4 * delta / (1 + 2 * delta)


def analyse(by_chrom):
    autosomal = [a / d for c, cnt in by_chrom.items() if is_autosome(c) for (a, d), n in cnt.items()
                 for _ in range(n)]
    centre = median(autosomal) if autosomal else 0.5
    rows = []
    for chrom, cnt in by_chrom.items():
        n = sum(cnt.values())
        depth_med = median([d for (a, d), k in cnt.items() for _ in range(k)])
        delta, lr = fit_delta(cnt, centre)
        rows.append({"chrom": chrom, "n_het": n, "median_depth": depth_med, "centre": round(centre, 4),
                     "delta": round(delta, 4), "lr_stat": round(lr, 2),
                     "trisomy_fraction_if_gain": round(trisomy_fraction(delta), 3),
                     "monosomy_fraction_if_loss": round(monosomy_fraction(delta), 3)})
    autos = [r["delta"] for r in rows if is_autosome(r["chrom"])]
    if len(autos) >= 3:
        med = median(autos)
        mad = median([abs(x - med) for x in autos]) * 1.4826 or 1e-4
        for r in rows:
            r["robust_z"] = round((r["delta"] - med) / mad, 2)
    return rows


def simulate(n_sites, depth, delta, centre=0.5, seed=1):
    """Synthetic heterozygous sites: Counter((alt, depth) -> n). Used for tests and detection limits."""
    rng = random.Random(seed)
    cnt = Counter()
    for _ in range(n_sites):
        p = centre - delta if rng.random() < 0.5 else centre + delta
        alt = sum(1 for _ in range(depth) if rng.random() < p)
        cnt[(alt, depth)] += 1
    return cnt


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("vcf")
    ap.add_argument("out")
    ap.add_argument("--min-depth", type=int, default=20)
    ap.add_argument("--max-depth", type=int, default=200)
    args = ap.parse_args()
    assert_outside_repo(args.out)
    rows = analyse(collect(args.vcf, args.min_depth, args.max_depth))
    fields = list(rows[0].keys()) if rows else ["chrom"]
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} chromosomes written")


if __name__ == "__main__":
    main()
