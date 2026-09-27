#!/usr/bin/env python3
"""Low-pass read depth per chromosome, with GC correction, as a second test for mosaic aneuploidy.

Full-genome alignment does not fit on this laptop, so reads from the start of the FASTQs are
aligned (minimap2 -ax sr in WSL, run by the caller) and counted in fixed bins:

  count  SAM on stdin -> per-bin counts of primary, mapped reads with MAPQ >= 20 (TSV)
  gc     reference FASTA -> GC fraction per bin (N bases excluded; streamed, one line at a time)
  ratios counts + GC tables -> per-chromosome GC-corrected depth ratio vs the autosomes

GC correction divides each bin by the median of autosomal bins in its GC stratum (0.01 steps).
A chromosome's ratio r implies a trisomic cell fraction 2(r - 1) or a monosomic fraction
2(1 - r). Variegated mosaicism spreads small gains and losses over many chromosomes, so
state the detection limit (the spread of autosomal ratios) with every result.

Run (in WSL, once the toolchain is installed):
  minimap2 -t 4 -ax sr REF.mmi R1_head.fq R2_head.fq | python3 -m pipeline.lowpass_depth count --bin 100000 > counts.tsv
  python3 -m pipeline.lowpass_depth gc REF.fa --bin 100000 > gc.tsv
  python3 -m pipeline.lowpass_depth ratios counts.tsv gc.tsv OUT.tsv     (outputs outside the repo)
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict
from statistics import median

from pipeline.common import assert_outside_repo, is_autosome, open_text


def count_sam(lines, bin_size=100_000, min_mapq=20):
    counts = Counter()
    for line in lines:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        flag, chrom, pos, mapq = int(f[1]), f[2], int(f[3]), int(f[4])
        if flag & (4 | 256 | 1024 | 2048) or chrom == "*" or mapq < min_mapq:
            continue
        counts[(chrom, (pos - 1) // bin_size * bin_size)] += 1
    return counts


def gc_corrected(bins):
    """bins: list of dicts with chrom, start, gc, reads. Returns list with 'corrected'."""
    strata = defaultdict(list)
    for b in bins:
        if b["gc"] is not None and is_autosome(b["chrom"]) and b["reads"] > 0:
            strata[round(b["gc"], 2)].append(b["reads"])
    medians = {g: median(v) for g, v in strata.items() if len(v) >= 3}
    out = []
    for b in bins:
        m = medians.get(round(b["gc"], 2)) if b["gc"] is not None else None
        if m:
            out.append(dict(b, corrected=b["reads"] / m))
    return out


def chromosome_ratios(bins):
    per_chrom = defaultdict(list)
    for b in gc_corrected(bins):
        per_chrom[b["chrom"]].append(b["corrected"])
    autos = [x for c, v in per_chrom.items() if is_autosome(c) for x in v]
    base = median(autos) if autos else 1.0
    rows = []
    for chrom, vals in per_chrom.items():
        r = median(vals) / base
        rows.append({"chrom": chrom, "bins": len(vals), "ratio": round(r, 4),
                     "trisomy_fraction_if_gain": round(max(0.0, 2 * (r - 1)), 3),
                     "monosomy_fraction_if_loss": round(max(0.0, 2 * (1 - r)), 3)})
    auto_ratios = [r["ratio"] for r in rows if is_autosome(r["chrom"])]
    if len(auto_ratios) >= 3:
        med = median(auto_ratios)
        mad = median(abs(x - med) for x in auto_ratios) * 1.4826 or 1e-4
        for r in rows:
            r["robust_z"] = round((r["ratio"] - med) / mad, 2)
    return rows


def gc_bins(lines, bin_size=100_000, min_called=0.5):
    """Yield (chrom, start, gc) per bin from FASTA lines; gc is None when under half the bin is called."""
    chrom, pos, gc, called = None, 0, 0, 0

    def flush():
        if chrom is not None and pos % bin_size:
            yield chrom, pos // bin_size * bin_size, (gc / called if called >= min_called * (pos % bin_size) else None)

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            yield from flush()
            chrom, pos, gc, called = line[1:].split()[0], 0, 0, 0
            continue
        for base in line.upper():
            if base in "GC":
                gc += 1
                called += 1
            elif base in "AT":
                called += 1
            pos += 1
            if pos % bin_size == 0:
                yield chrom, pos - bin_size, (gc / called if called >= min_called * bin_size else None)
                gc, called = 0, 0
    yield from flush()


def join_bins(counts_path, gc_path):
    with open(gc_path, encoding="utf-8", newline="") as fh:
        gc = {(r["chrom"], int(r["start"])): (float(r["gc"]) if r["gc"] not in ("", "NA") else None)
              for r in csv.DictReader(fh, delimiter="\t")}
    with open(counts_path, encoding="utf-8", newline="") as fh:
        return [{"chrom": r["chrom"], "start": int(r["start"]), "reads": int(r["reads"]),
                 "gc": gc.get((r["chrom"], int(r["start"])))}
                for r in csv.DictReader(fh, delimiter="\t")]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("count")
    c.add_argument("--bin", type=int, default=100_000)
    c.add_argument("--min-mapq", type=int, default=20)
    g = sub.add_parser("gc")
    g.add_argument("fasta")
    g.add_argument("--bin", type=int, default=100_000)
    r = sub.add_parser("ratios")
    r.add_argument("counts")
    r.add_argument("gc")
    r.add_argument("out")
    a = ap.parse_args()
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    if a.cmd == "count":
        w.writerow(["chrom", "start", "reads"])
        for (chrom, start), n in sorted(count_sam(sys.stdin, a.bin, a.min_mapq).items()):
            w.writerow([chrom, start, n])
    elif a.cmd == "gc":
        w.writerow(["chrom", "start", "gc"])
        with open_text(a.fasta) as fh:
            for chrom, start, frac in gc_bins(fh, a.bin):
                w.writerow([chrom, start, "NA" if frac is None else round(frac, 4)])
    else:
        assert_outside_repo(a.out)
        rows = chromosome_ratios(join_bins(a.counts, a.gc))
        with open(a.out, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
