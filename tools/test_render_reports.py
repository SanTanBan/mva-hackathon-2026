#!/usr/bin/env python3
"""Tests for render_reports.py. All content is synthetic; none is patient data.

Variant-like strings are assembled at run time so this file never contains one.

Run: python tools/test_render_reports.py
"""
import importlib
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


def fake_variant(sep=":"):
    return sep.join(["chr" + "7", "1234" + "5678", "A", "G"])


SOURCE = """# Methods report

Public method text.

<!-- private:start -->
Patient-level result: SYNTHETIC-RESULT-42.
<!-- private:end -->

Inline: the pair was [[private: SYNTHETIC-PAIR-7]] and phase was scored.
"""


class RenderTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.submissions, self.public = root / "submissions", root / "public"
        self.terms = root / "terms.txt"
        self.terms.write_text("# synthetic\nSYNTHID9\n", encoding="utf-8")
        os.environ["MVA_SUBMISSIONS_DIR"] = str(self.submissions)
        os.environ["MVA_PRIVATE_TERMS"] = str(self.terms)
        import render_reports
        self.rr = importlib.reload(render_reports)

    def tearDown(self):
        for key in ("MVA_SUBMISSIONS_DIR", "MVA_PRIVATE_TERMS"):
            os.environ.pop(key, None)
        self.tmp.cleanup()

    def write_source(self, text):
        path = Path(self.tmp.name) / "report.src.md"
        path.write_text(text, encoding="utf-8")
        return path

    def run_render(self, text, pdf=False):
        return self.rr.render(self.write_source(text), "track1", "methods_report", "2026-01-01", self.public, pdf)

    def test_public_withholds_private_parts(self):
        self.run_render(SOURCE)
        public = (self.public / "methods_report.md").read_text(encoding="utf-8")
        self.assertNotIn("SYNTHETIC-RESULT-42", public)
        self.assertNotIn("SYNTHETIC-PAIR-7", public)
        self.assertEqual(public.count("Withheld from the public version"), 2)
        self.assertIn("Public method text.", public)

    def test_full_keeps_content_without_markers(self):
        self.run_render(SOURCE)
        full = (self.submissions / "track1" / "2026-01-01" / "methods_report.md").read_text(encoding="utf-8")
        self.assertIn("SYNTHETIC-RESULT-42", full)
        self.assertIn("the pair was SYNTHETIC-PAIR-7 and", full)
        self.assertNotIn("private:", full)

    def test_variant_token_outside_private_is_refused(self):
        for sep in (":", "-", "_"):
            with self.assertRaises(SystemExit):
                self.run_render(SOURCE + "\nLeaked " + fake_variant(sep) + "\n")

    def test_variant_token_inside_private_is_allowed(self):
        self.run_render(SOURCE + "\n[[private: " + fake_variant() + "]]\n")
        full = (self.submissions / "track1" / "2026-01-01" / "methods_report.md").read_text(encoding="utf-8")
        self.assertIn(fake_variant(), full)

    def test_protected_term_is_refused(self):
        with self.assertRaises(SystemExit):
            self.run_render(SOURCE + "\nSample SYNTHID9 was used.\n")

    def test_vcf_line_is_refused(self):
        line = "\t".join(["7", "12345678", ".", "A", "G"])
        with self.assertRaises(SystemExit):
            self.run_render(SOURCE + "\n" + line + "\n")

    def test_unbalanced_marker_is_refused(self):
        with self.assertRaises(SystemExit):
            self.run_render(SOURCE + "\n<!-- private:start -->\nforgot to close\n")

    def test_gene_coordinates_are_not_mistaken_for_variants(self):
        self.run_render(SOURCE + "\nGene span on chromosome 15: 40,160,984-40,221,137 (GRCh38).\n")

    def test_full_version_must_be_outside_repo(self):
        os.environ["MVA_SUBMISSIONS_DIR"] = str(Path(__file__).resolve().parents[1] / "submissions")
        rr = importlib.reload(self.rr)
        with self.assertRaises(SystemExit):
            rr.render(self.write_source(SOURCE), "track1", "x", "2026-01-01", self.public, False)

    @unittest.skipUnless(os.environ.get("MVA_TEST_PDF") == "1", "set MVA_TEST_PDF=1 to test Edge PDF printing")
    def test_long_table_tokens_get_break_points(self):
        body = ("<p>stop_gained&amp;splice_region_variant</p>"
                "<table><tr><th>Rank</th><td>stop_gained&amp;splice_region_variant</td><td>ACGTACGTACGTACGTACGT</td></tr></table>")
        out = self.rr.break_long_tokens(body)
        self.assertIn("<p>stop_gained&amp;splice_region_variant</p>", out)  # outside tables: untouched
        self.assertIn("<th>Rank</th>", out)  # short words stay whole
        self.assertIn("stop_<wbr>gained&amp;<wbr>splice_<wbr>region_<wbr>variant", out)
        self.assertIn("ACGTACGTACGT<wbr>ACGTACGT", out)
        self.assertEqual(out.replace("<wbr>", ""), body)
        wide = "<table><tr>" + "<th>c</th>" * 8 + "</tr></table>"
        self.assertIn('<table class="wide">', self.rr.break_long_tokens(wide))
        self.assertNotIn("wide", self.rr.break_long_tokens(wide.replace("<th>c</th>", "", 1)))

    def test_pdf_written_when_edge_available(self):
        if not self.rr.find_edge():
            self.skipTest("Microsoft Edge not found")
        outputs = self.run_render(SOURCE, pdf=True)
        pdfs = [p for p in outputs if p.suffix == ".pdf"]
        self.assertEqual(len(pdfs), 2)
        for p in pdfs:
            self.assertEqual(p.read_bytes()[:5], b"%PDF-")


if __name__ == "__main__":
    unittest.main(verbosity=2)
