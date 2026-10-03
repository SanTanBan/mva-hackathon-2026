#!/usr/bin/env python3
"""Synthetic tests for tools/build_panel_genes.py. No network, no patient data.

Run: python -m unittest tools.test_build_panel_genes
"""
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import build_panel_genes as b  # noqa: E402


def gene(name, level, moi, hgnc=None, symbol=None):
    return {"entity_name": name, "confidence_level": level, "mode_of_inheritance": moi,
            "gene_data": {"hgnc_id": hgnc, "gene_symbol": symbol or name}}


class PanelGenesTest(unittest.TestCase):
    def test_moi(self):
        self.assertTrue(b.allows_recessive("BIALLELIC, autosomal or pseudoautosomal"))
        self.assertTrue(b.allows_recessive("BOTH monoallelic and biallelic, autosomal or pseudoautosomal"))
        self.assertTrue(b.allows_recessive("X-LINKED: hemizygous mutation in males, biallelic mutations in females"))
        self.assertTrue(b.allows_recessive(""))
        self.assertTrue(b.allows_recessive(None))
        self.assertFalse(b.allows_recessive("MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted"))

    def test_select(self):
        panels = [
            {"id": 1, "genes": [gene("BUB1B", "3", "BIALLELIC, autosomal"),
                                gene("AARS", "3", "BIALLELIC, autosomal", hgnc="HGNC:20"),
                                gene("TP53", "3", "MONOALLELIC, autosomal"),
                                gene("AMBER1", "2", "BIALLELIC, autosomal")]},
            {"id": 2, "genes": [gene("BUB1B", 3, "BIALLELIC, autosomal")]},
        ]
        got = b.select_genes(panels, {"HGNC:20": "AARS1"})
        self.assertEqual(got, {"BUB1B": [1, 2], "AARS1": [1]})

    def test_hgnc_map_and_write(self):
        with tempfile.TemporaryDirectory(dir=r"D:\DevCache\tmp") as d:
            gff = Path(d) / "x.gff3.gz"
            with gzip.open(gff, "wt", encoding="utf-8") as fh:
                fh.write("##gff-version 3\n")
                fh.write("1\tensembl\tgene\t1\t9\t.\t+\t.\tID=gene:E1;Name=AARS1;description=x "
                         "[Source:HGNC Symbol%3BAcc:HGNC:20];gene_id=E1\n")
                fh.write("1\tensembl\tmRNA\t1\t9\t.\t+\t.\tID=t;Name=AARS1-201;Acc:HGNC:99\n")
            self.assertEqual(b.hgnc_to_symbol(gff), {"HGNC:20": "AARS1"})
            out = Path(d) / "genes.txt"
            b.write(out, [{"id": 1, "version": "2.0", "name": "Test panel"}], {"B": [1], "A": [1]}, "2026-10-04")
            lines = out.read_text(encoding="utf-8").splitlines()
            self.assertEqual([ln for ln in lines if not ln.startswith("#")], ["A", "B"])
            self.assertIn("# panel 1 v2.0: Test panel", lines)


if __name__ == "__main__":
    unittest.main()
