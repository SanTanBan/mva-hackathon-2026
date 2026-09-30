#!/usr/bin/env python3
"""Tests for the local splice-map lookup, on a synthetic map. No patient data.

Run: python -m unittest pipeline.tests.test_splice_map
"""
import gzip
import tempfile
import unittest
from pathlib import Path

from pipeline import splice_map

HEADER = ["chrom", "pos", "ref", "alt", "hgvs_c", "region", "exon_or_intron", "gain", "gain_offset", "loss", "loss_offset"]
ROWS = [
    ["chr15", 1000, "A", "G", "c.2386-11A>G", "intron", 18, 0.7577, 1, -0.7145, 11],
    ["chr15", 2000, "T", "C", "c.1402-250T>C", "intron", 10, 0.62, 3, 0.0, 0],
    ["chr15", 3000, "G", "T", "c.1402-1G>T", "intron", 10, 0.3046, 5, -0.857, 1],
    ["chr15", 4000, "C", "T", "c.1500C>T", "exon", 11, 0.01, 0, -0.02, 2],
]


class SpliceMapTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        path = Path(self.tmp.name) / "map.tsv.gz"
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write("\t".join(HEADER) + "\n")
            for r in ROWS:
                fh.write("\t".join(str(x) for x in r) + "\n")
        self.index = splice_map.load(path)

    def tearDown(self):
        self.tmp.cleanup()

    def test_lookup_by_hgvs_and_position(self):
        a = splice_map.lookup("c.2386-11A>G", self.index)
        b = splice_map.lookup("chr15:1000:a:g", self.index)
        self.assertEqual(a, b)
        self.assertEqual((a["distance_class"], a["gain_call"], a["loss_call"]), ("5-14 nt", "strong", "strong"))
        self.assertIn("too close", a["note"])

    def test_deep_intronic_candidate(self):
        r = splice_map.lookup("c.1402-250T>C", self.index)
        self.assertEqual(r["distance_class"], "deep intronic")
        self.assertIn("candidate for splice correction", r["note"])

    def test_canonical_and_exonic(self):
        self.assertEqual(splice_map.lookup("c.1402-1G>T", self.index)["distance_class"], "canonical")
        r = splice_map.lookup("15:4000:C:T", self.index)
        self.assertEqual((r["distance_class"], r["gain_call"], r["loss_call"], r["note"]), ("exonic", "none", "none", None))

    def test_missing_variant(self):
        self.assertIsNone(splice_map.lookup("c.1500C>G", self.index))

    def test_branch_window(self):
        self.assertEqual(splice_map.distance_class("c.1402-30A>G", "intron"), "branch-point window")
        self.assertEqual(splice_map.distance_class("c.1288+30A>G", "intron"), "15-100 nt")


if __name__ == "__main__":
    unittest.main()
