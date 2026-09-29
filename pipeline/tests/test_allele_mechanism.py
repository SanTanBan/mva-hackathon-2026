#!/usr/bin/env python3
"""Tests for the pre-registered D08 allele classifier.

Named variants are published ClinVar or literature alleles, used as worked examples (for the
literature substitutions, the c. name is derived from NM_001211.6); the rest are invented. None is
patient data, and nothing here says which allele the child carries.

Run: python -m unittest pipeline.tests.test_allele_mechanism
"""
import unittest

from pipeline import allele_mechanism as am
from pipeline.common import REPO


class ExonTableTest(unittest.TestCase):
    def test_exons_tile_the_coding_sequence(self):
        self.assertEqual(am.EXONS[0][1], 1)
        self.assertEqual(am.EXONS[-1][2], am.CDS_LEN)
        for (_, _, end, _), (_, start, _, _) in zip(am.EXONS, am.EXONS[1:]):
            self.assertEqual(start, end + 1)

    def test_in_frame_internal_exons(self):
        in_frame = {e for e, s, t, _ in am.EXONS[1:-1] if (t - s + 1) % 3 == 0}
        self.assertEqual(in_frame, {2, 3, 17, 19})
        self.assertEqual(set(am.IN_FRAME_SKIP), in_frame)

    def test_table_renders(self):
        table = am.exon_table()
        self.assertIn("| 17 | 141 | in-frame | 715-762 | unlikely eligible |", table)
        self.assertIn("| 23 | 196 | n/a | 986-1050 |", table)

    def test_no_bub1b_exon_is_likely_eligible_for_skipping(self):
        # Pre-registered conclusion: exon 17 is the best candidate but R727C argues against it.
        self.assertNotIn("likely eligible", [v[0] for v in am.IN_FRAME_SKIP.values()])
        self.assertIn("R727C", am.IN_FRAME_SKIP[17][1])


class ParseTest(unittest.TestCase):
    def test_locations(self):
        self.assertEqual(am.parse_c("c.2210T>G"), {"region": "exon", "c": 2210, "exon": 17})
        loc = am.parse_c("c.2386-11A>G")
        self.assertEqual((loc["intron"], loc["side"], loc["distance"]), (18, "acceptor", 11))
        loc = am.parse_c("c.1288+5G>A")
        self.assertEqual((loc["intron"], loc["side"], loc["distance"]), (9, "donor", 5))
        self.assertEqual(am.parse_c("c.-60C>T")["region"], "5' UTR or upstream")
        with self.assertRaises(ValueError):
            am.parse_c("c.2200+3A>G")  # not an exon boundary

    def test_protein_changes(self):
        self.assertEqual(am.protein_change("p.Ile909Thr"), ("I", 909, "T"))
        self.assertEqual(am.protein_change("L1012P"), ("L", 1012, "P"))
        self.assertEqual(am.protein_change("p.Leu737Ter"), ("L", 737, "*"))
        self.assertEqual(am.protein_change("p.Arg97ProfsTer23"), ("fs", 97, 119))
        self.assertEqual(am.protein_change("p.Ser738fs"), ("fs", 738, None))


class NmdTest(unittest.TestCase):
    def test_rules(self):
        self.assertEqual(am.nmd(2209), "predicted")
        self.assertEqual(am.nmd(2930), "escapes (within 50 nt of the last exon junction)")
        self.assertEqual(am.nmd(3000), "escapes (last exon)")
        self.assertEqual(am.nmd(100), "likely escapes (first 100 coding nt)")
        self.assertEqual(am.nmd(169), "partial (101-200 nt from the start codon)")


class AlleleTest(unittest.TestCase):
    def test_exon17_stop_is_null(self):
        r = am.classify("c.2210T>G", "stop_gained", "p.Leu737Ter")
        self.assertEqual(r["nmd"], "predicted")
        self.assertEqual(r["mechanism_class"], "abundance")
        self.assertEqual(r["aso"]["exon_skipping"][0], "unlikely eligible")
        self.assertEqual(r["aso"]["splice_correction"][0], "not applicable")
        self.assertEqual(r["aso"]["knockdown"][0], "not eligible")

    def test_out_of_frame_and_boundary_exons(self):
        self.assertEqual(am.classify("c.1833del", "frameshift_variant", "p.Phe611fs")["aso"]["exon_skipping"][0], "not eligible")
        self.assertEqual(am.classify("c.100C>T", "stop_gained", "p.Gln34Ter")["aso"]["exon_skipping"][0], "not eligible")
        self.assertEqual(am.classify("c.3010C>T", "stop_gained", "p.Gln1004Ter")["aso"]["exon_skipping"][0], "not eligible")

    def test_frameshift_without_stop_position_is_flagged(self):
        r = am.classify("c.2208_2211dup", "frameshift_variant", "p.Ser738fs")
        self.assertIn("approximate", r["note"])
        self.assertEqual(r["aso"]["exon_skipping"][0], "unlikely eligible")

    def test_tested_substitutions(self):
        r = am.classify("c.3035T>C", "missense_variant", "p.Leu1012Pro")
        self.assertEqual(r["mechanism_class"], "abundance")
        self.assertEqual(r["aso"]["exon_skipping"][0], "not eligible")  # last exon
        self.assertEqual(r["aso"]["upregulation_untiered"][0], "candidate")
        self.assertEqual(am.classify("c.2441G>A", "missense_variant", "p.Arg814His")["aso"]["exon_skipping"][0], "not eligible")
        r = am.classify("c.464A>G", "missense_variant", "p.Tyr155Cys")
        self.assertEqual(r["mechanism_class"], "functional")
        self.assertEqual(r["aso"]["upregulation_untiered"][0], "unlikely to help")

    def test_untested_substitutions_by_position(self):
        self.assertEqual(am.classify("c.2995G>A", "missense_variant", "p.Gly999Arg")["mechanism_class"], "abundance (provisional)")
        self.assertEqual(am.classify("c.898A>G", "missense_variant", "p.Ser300Gly")["mechanism_class"],
                         "functional or no defect (provisional)")
        self.assertEqual(am.classify("c.1798G>A", "missense_variant", "p.Gly600Arg")["mechanism_class"], "unclassified")
        # A missense in exon 17 argues that the exon matters, so skipping it is downgraded.
        self.assertEqual(am.classify("c.2200G>A", "missense_variant", "p.Glu734Lys")["aso"]["exon_skipping"][0], "unlikely eligible")

    def test_splice_correction_tiers(self):
        self.assertEqual(am.classify("c.1402-1G>T", "splice_acceptor_variant")["aso"]["splice_correction"][0], "not eligible")
        # Published new acceptor 11 nt upstream (Rio Frio 2010): under the recommended 15-nt distance.
        self.assertEqual(am.classify("c.2386-11A>G", "splice_region_variant", rna="cryptic_site")["aso"]["splice_correction"][0],
                         "unlikely eligible")
        deep = am.classify("c.1402-250A>G", "intron_variant", rna="pseudoexon")["aso"]["splice_correction"]
        self.assertEqual(deep[0], "likely eligible")
        pending = am.classify("c.1402-250A>G", "intron_variant")["aso"]["splice_correction"]
        self.assertEqual(pending[0], "unable to assess")
        self.assertIn("would be 'likely eligible'", pending[1])
        self.assertEqual(am.classify("c.1402-30A>G", "intron_variant", rna="cryptic_site")["aso"]["splice_correction"][0],
                         "unlikely eligible")  # branch-point window
        self.assertEqual(am.classify("c.1288+5G>A", "splice_region_variant", rna="exon_skip")["aso"]["splice_correction"][0],
                         "not eligible")
        self.assertEqual(am.classify("c.1740C>T", "synonymous_variant", rna="cryptic_site")["aso"]["splice_correction"][0],
                         "unlikely eligible")  # 6 nt into a 275-nt exon
        self.assertEqual(am.classify("c.1850C>T", "synonymous_variant", rna="cryptic_site")["aso"]["splice_correction"][0],
                         "likely eligible")
        with self.assertRaises(ValueError):
            am.classify("c.1402-250A>G", "intron_variant", rna="maybe")


class GenotypeTest(unittest.TestCase):
    def test_exon17_null_with_kinase_missense(self):
        a = am.classify("c.2210T>G", "stop_gained", "p.Leu737Ter")
        b = am.classify("c.2726T>C", "missense_variant", "p.Ile909Thr")
        rules = am.genotype(a, b)
        self.assertIn("H4 stays rank 4; first test: is exon 17-skipped BUBR1 stable in patient cells (plasmid or ASO)?", rules)
        self.assertFalse(any(r.startswith("H4 -> rank 1") for r in rules))
        self.assertTrue(any(r.startswith("H1 and H2 keep") for r in rules))

    def test_null_with_functional_class(self):
        a = am.classify("c.1833del", "frameshift_variant", "p.Phe611fs")
        b = am.classify("c.464A>G", "missense_variant", "p.Tyr155Cys")
        rules = am.genotype(a, b)
        self.assertTrue(any(r.startswith("goal A weakened") for r in rules))
        self.assertFalse(any(r.startswith("H4 -> rank 1") for r in rules))

    def test_deep_intronic_awaiting_rna(self):
        a = am.classify("c.1833del", "frameshift_variant", "p.Phe611fs")
        b = am.classify("c.1402-250A>G", "intron_variant")
        self.assertIn("H4 stays rank 4: RNA from patient cells becomes the first experiment", am.genotype(a, b))

    def test_two_nulls_flagged(self):
        a = am.classify("c.1833del", "frameshift_variant", "p.Phe611fs")
        b = am.classify("c.1150C>T", "stop_gained", "p.Gln384Ter")
        self.assertTrue(any(r.startswith("both alleles look null") for r in am.genotype(a, b)))


class OutputGuardTest(unittest.TestCase):
    def test_tsv_output_refuses_repo(self):
        with self.assertRaises(SystemExit):
            am.assert_outside_repo(REPO / "reports" / "d08_alleles.json")


if __name__ == "__main__":
    unittest.main()
