#!/usr/bin/env python3
"""Synthetic tests for pipeline/track1_rank.py and the same-site phase rule. Nothing here is patient data.

Run: python -m unittest discover -s pipeline/tests -t .
"""
import csv
import tempfile
import unittest
from pathlib import Path

from pipeline import comphet, track1_rank
from tools.validate_track1_csv import validate


def v(pos, alt, gene="BUB1B", impact="HIGH", gt="0/1", clnsig="", pgt="", pid="", af="", ref="C"):
    return {"chrom": "chr15", "pos": str(pos), "ref": ref, "alt": alt, "gene": gene,
            "consequence": "stop_gained" if impact == "HIGH" else "missense", "impact": impact,
            "pop_af": af, "gt": gt, "pgt": pgt, "pid": pid, "clnsig": clnsig, "aa_change": ""}


class PhaseTest(unittest.TestCase):
    def test_same_site_is_trans(self):
        self.assertEqual(comphet.phase(v(100, "T"), v(100, "A")), "trans")
        self.assertEqual(comphet.phase(v(100, "T"), v(200, "A")), "unknown")


class ScoreTest(unittest.TestCase):
    def test_clinvar_classes(self):
        self.assertEqual(track1_rank.clinvar_class("Pathogenic/Likely_pathogenic"), "plp")
        self.assertEqual(track1_rank.clinvar_class("Benign/Likely_benign"), "blb")
        self.assertEqual(track1_rank.clinvar_class("Conflicting_classifications_of_pathogenicity"), "conflicting")
        self.assertEqual(track1_rank.clinvar_class("Uncertain_significance"), "")

    def test_ranking(self):
        rows = [
            v(100, "T", clnsig="Pathogenic", pgt="0|1", pid="P1"),
            v(300, "A", impact="MODERATE", pgt="1|0", pid="P1"),
            v(500, "G", gene="GENEX", gt="1/1"),
            v(700, "G", gene="GENEY", impact="MODERATE"),
            v(900, "G", gene="GENEY", impact="MODERATE"),
            v(1100, "G", clnsig="Benign", af="0.2"),            # common: filtered out
        ]
        ranked = track1_rank.score_candidates(rows, panel=["GENEY"])
        top = ranked[0]
        self.assertEqual((top[1]["pos"], top[2]["pos"]), ("100", "300"))
        self.assertIn("phase trans", top[3])
        epcrs = [r[0] for r in ranked]
        self.assertEqual(epcrs, sorted(epcrs, reverse=True))
        self.assertEqual(len(set(epcrs)), len(epcrs))
        self.assertLess(sum(epcrs), 1.0)
        self.assertNotIn("1100", [r[1]["pos"] for r in ranked])
        # BUB1B alleles also appear as single-allele rows only when unpaired: both are paired here
        self.assertFalse(any("second allele not seen" in r[3] for r in ranked))

    def test_single_allele_row_and_csv(self):
        rows = [v(100, "T", clnsig="Likely_pathogenic"), v(500, "G", gene="GENEX", gt="1/1", impact="MODERATE")]
        ranked = track1_rank.score_candidates(rows)
        self.assertIn("second allele not seen", ranked[0][3])
        with tempfile.TemporaryDirectory(dir=r"D:\DevCache\tmp") as d:
            path = Path(d) / "SantanuBanerjee_synthetic.csv"
            track1_rank.write_csv(path, "SYNTH-1", ranked)
            errors, _ = validate(path, hf_user="SantanuBanerjee")
            self.assertEqual(errors, [])
            with open(path, encoding="utf-8") as fh:
                first = next(csv.DictReader(fh))
            self.assertEqual(first["chrom_2"], "")


if __name__ == "__main__":
    unittest.main()
