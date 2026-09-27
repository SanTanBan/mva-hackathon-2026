#!/usr/bin/env python3
"""Synthetic tests for the ACMG SF v3.3 screen. Every row is invented; none is patient data.

Run: python -m unittest discover -s pipeline/tests -t .
"""
import unittest

from pipeline import sf_screen


def v(gene, cls="Pathogenic", gt="0/1", hgvsp=""):
    return {"gene": gene, "classification": cls, "gt": gt, "hgvsp": hgvsp}


class SfScreenTest(unittest.TestCase):
    def status(self, rows, sex="unknown"):
        return {f["gene"]: f["status"] for f in sf_screen.screen(rows, sex)}

    def test_gene_list_matches_published_count(self):
        self.assertEqual(len(sf_screen.SF_V33), 84)

    def test_rules(self):
        s = self.status([
            v("BRCA2"),                                        # AD: reportable
            v("MUTYH"),                                        # AR single: carrier
            v("GAA"), v("GAA", cls="Likely pathogenic"),       # AR two: reportable, phase to check
            v("HFE", gt="0/1", hgvsp="p.Cys282Tyr"),           # HFE het: not reportable
            v("ATP7B", gt="1/1"),                              # AR homozygote: reportable
            v("GLA"),                                          # XL with unknown sex
            v("TTN"),                                          # DCM: policy-text review
            v("BUB1B"),                                        # not an SF gene: ignored
            v("APC", cls="Uncertain significance"),            # VUS: ignored
        ])
        self.assertEqual(s["BRCA2"], "reportable")
        self.assertEqual(s["MUTYH"], "carrier (not reportable)")
        self.assertTrue(s["GAA"].startswith("reportable (biallelic"))
        self.assertTrue(s["HFE"].startswith("not reportable"))
        self.assertTrue(s["ATP7B"].startswith("reportable"))
        self.assertEqual(s["GLA"], "depends on stated sex (X-linked)")
        self.assertIn("review", s["TTN"])
        self.assertNotIn("BUB1B", s)
        self.assertNotIn("APC", s)

    def test_hfe_c282y_homozygote(self):
        self.assertEqual(self.status([v("HFE", gt="1/1", hgvsp="p.Cys282Tyr")])["HFE"], "reportable")

    def test_x_linked_needs_stated_sex(self):
        self.assertEqual(self.status([v("ABCD1")], sex="female")["ABCD1"], "carrier (not reportable)")
        self.assertEqual(self.status([v("ABCD1")], sex="male")["ABCD1"], "reportable")
        with self.assertRaises(ValueError):
            sf_screen.screen([v("GLA")], sex="inferred")

    def test_public_counts_only(self):
        counts = sf_screen.public_counts(sf_screen.screen([v("BRCA2"), v("MUTYH")]))
        self.assertEqual(counts, {"sf_genes_screened": 84, "genes_with_plp": 2,
                                  "reportable_or_review": 1, "carrier_or_not_reportable": 1})


if __name__ == "__main__":
    unittest.main(verbosity=2)
