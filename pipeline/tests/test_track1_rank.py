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
        # The pathogenic allele also gets a single-allele row, ranked below its read-backed trans pair
        self.assertIn("second allele not seen", ranked[1][3])
        self.assertEqual((ranked[1][1]["pos"], ranked[1][2]), ("100", None))

    def test_rows_looked_up_by_gene(self):
        # Each BUB1B allele also has a MODIFIER row for a neighbouring gene, listed after it.
        rows = [
            v(100, "T", pgt="0|1", pid="P1"), v(300, "A", pgt="1|0", pid="P1"),
            v(100, "T", gene="GENEZ", impact="MODIFIER"), v(300, "A", gene="GENEZ", impact="MODIFIER"),
        ]
        ranked = track1_rank.score_candidates(rows)
        self.assertIn("BUB1B compound het, phase trans: stop_gained", ranked[0][3])
        self.assertGreater(ranked[0][0], 0.5)             # 0.64 / (0.64 + 0.08 + 0.08 + 0.25)

    def test_overlapping_genes_written_once(self):
        rows = [v(100, "T", gene="GENEA"), v(300, "A", gene="GENEA"),
                v(100, "T", gene="GENEB"), v(300, "A", gene="GENEB")]
        self.assertEqual(len(track1_rank.score_candidates(rows)), 1)

    def test_filter_fail_halves(self):
        self.assertEqual(track1_rank.allele_score(v(100, "T")), 0.8)
        self.assertEqual(track1_rank.allele_score({**v(100, "T"), "filter": "LowQual"}), 0.4)

    def test_weak_candidates_do_not_dilute(self):
        strong = [v(100, "T", clnsig="Pathogenic", pgt="0|1", pid="P1"),
                  v(300, "A", clnsig="Pathogenic", pgt="1|0", pid="P1")]
        noise = [v(1000 + i, "G", impact="MODIFIER") for i in range(150)]  # unknown-frequency intronic calls
        ranked = track1_rank.score_candidates(strong + noise)
        self.assertEqual(len(ranked), 10)
        self.assertGreater(ranked[0][0], 0.6)

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
