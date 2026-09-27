#!/usr/bin/env python3
"""Synthetic tests for k-mer baiting and low-pass depth. All sequences are random; none is patient data.

Run: python -m unittest discover -s pipeline/tests -t .
"""
import gzip
import random
import tempfile
import unittest
from pathlib import Path

from pipeline import kmer_bait, lowpass_depth


def random_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


class KmerBaitTest(unittest.TestCase):
    def test_keeps_target_pairs_and_masks_low_complexity(self):
        rng = random.Random(7)
        target = random_seq(rng, 1500) + "A" * 80 + random_seq(rng, 1500)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "target.fa").write_text(">locus\n" + target + "\n", encoding="utf-8")
            r1, r2 = [], []

            def add(s1, s2, name):
                r1.append(f"@{name}/1\n{s1}\n+\n{'I' * len(s1)}\n")
                r2.append(f"@{name}/2\n{s2}\n+\n{'I' * len(s2)}\n")
            for i in range(5):                        # true locus pairs, mate 2 on the reverse strand
                start = rng.randrange(0, 1300)
                add(target[start:start + 100], kmer_bait.revcomp(target[start + 150:start + 250]), f"on{i}")
            for i in range(5):                        # unrelated pairs
                add(random_seq(rng, 100), random_seq(rng, 100), f"off{i}")
            add("A" * 100, "T" * 100, "polyA")        # matches only the masked homopolymer
            for name, recs in (("R1.fq.gz", r1), ("R2.fq.gz", r2)):
                with gzip.open(tmp / name, "wt") as fh:
                    fh.writelines(recs)
            stats = kmer_bait.bait(tmp / "target.fa", tmp / "R1.fq.gz", tmp / "R2.fq.gz", tmp / "baited")
            self.assertEqual(stats["pairs_seen"], 11)
            self.assertEqual(stats["pairs_kept"], 5)
            self.assertGreater(stats["masked_kmers"], 0)
            with gzip.open(tmp / "baited_R1.fq.gz", "rt") as fh:
                names = [l.strip() for l in fh if l.startswith("@")]
            self.assertTrue(all(n.startswith("@on") for n in names))

    def test_low_complexity(self):
        self.assertTrue(kmer_bait.low_complexity("A" * 31))
        self.assertTrue(kmer_bait.low_complexity(("CA" * 16)[:31]))
        self.assertFalse(kmer_bait.low_complexity("ACGTTGCAAGCTTAGCCGATCGGATCCAGTA"))


class LowpassDepthTest(unittest.TestCase):
    def test_gc_corrected_ratio_recovers_mosaic_gain(self):
        rng = random.Random(3)
        bins = []
        for c in range(1, 23):
            copy_ratio = 1.2 if c == 8 else 1.0      # 40% trisomic cells: ratio (2 + 0.4) / 2
            for i in range(60):
                gc = 0.35 + 0.01 * (i % 26)
                expected = 1000 * (1 + 2 * (gc - 0.45)) * copy_ratio
                bins.append({"chrom": f"chr{c}", "start": i * 100_000, "gc": gc,
                             "reads": max(1, int(rng.gauss(expected, expected ** 0.5)))})
        rows = {r["chrom"]: r for r in lowpass_depth.chromosome_ratios(bins)}
        self.assertAlmostEqual(rows["chr8"]["ratio"], 1.2, delta=0.03)
        self.assertAlmostEqual(rows["chr8"]["trisomy_fraction_if_gain"], 0.4, delta=0.06)
        self.assertAlmostEqual(rows["chr1"]["ratio"], 1.0, delta=0.03)
        self.assertEqual(max(rows.values(), key=lambda r: r["robust_z"])["chrom"], "chr8")

    def test_gc_bins_and_join(self):
        fasta = [">c1 synthetic", "GGGGCCCCAT", "ATATATATNN", "GC", ">c2", "NNNNNNNNNN"]
        bins = list(lowpass_depth.gc_bins(fasta, bin_size=10))
        self.assertEqual(bins[0], ("c1", 0, 0.8))          # 8 G/C of 10 called
        self.assertEqual(bins[1], ("c1", 10, 0.0))         # 8 called (A/T), 2 N
        self.assertEqual(bins[2], ("c1", 20, 1.0))         # partial last bin "GC"
        self.assertEqual(bins[3], ("c2", 0, None))         # all N: no GC value
        with tempfile.TemporaryDirectory() as tmp:
            counts, gc = Path(tmp) / "counts.tsv", Path(tmp) / "gc.tsv"
            counts.write_text("chrom\tstart\treads\nc1\t0\t5\nc1\t10\t7\n", encoding="utf-8")
            gc.write_text("chrom\tstart\tgc\nc1\t0\t0.8\nc1\t10\tNA\n", encoding="utf-8")
            joined = lowpass_depth.join_bins(counts, gc)
        self.assertEqual([(b["start"], b["reads"], b["gc"]) for b in joined], [(0, 5, 0.8), (10, 7, None)])

    def test_count_sam_filters(self):
        def sam(flag, chrom, pos, mapq):
            return "\t".join(["read", str(flag), chrom, str(pos), str(mapq), "100M", "*", "0", "0", "A", "I"])
        lines = ["@HD\tVN:1.6", sam(0, "chr1", 150, 60), sam(16, "chr1", 99_999, 60), sam(0, "chr1", 100_001, 60),
                 sam(256, "chr1", 10, 60), sam(4, "*", 0, 0), sam(0, "chr1", 20, 5), sam(1024, "chr1", 30, 60)]
        counts = lowpass_depth.count_sam(lines, bin_size=100_000)
        self.assertEqual(dict(counts), {("chr1", 0): 2, ("chr1", 100_000): 1})


if __name__ == "__main__":
    unittest.main(verbosity=2)
