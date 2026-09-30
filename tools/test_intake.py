#!/usr/bin/env python3
"""Synthetic tests for tools/fetch_vcf_first.py (D00) and tools/d01_intake.py (D01). No patient data.

Run: python -m unittest tools.test_intake
"""
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import d01_intake, fetch_vcf_first  # noqa: E402


class FetchTest(unittest.TestCase):
    def test_classify(self):
        c = fetch_vcf_first.classify
        self.assertEqual(c("x/SAMPLE_R1.fastq.gz"), "fastq")
        self.assertEqual(c("calls.vcf.gz"), "vcf")
        self.assertEqual(c("calls.vcf.gz.tbi"), "vcf-index")
        self.assertEqual(c("reads.cram"), "alignment")
        self.assertEqual(c("phenotype.pdf"), "phenotype-or-other")
        self.assertEqual(c(".gitattributes"), "repo-metadata")

    def test_fastq_pairs(self):
        files = ["A_L001_R1_001.fastq.gz", "A_L001_R2_001.fastq.gz", "A_L002_R1_001.fastq.gz", "A_L002_R2_001.fastq.gz"]
        self.assertEqual(fetch_vcf_first.fastq_pairs(files),
                         [["A_L001_R1_001.fastq.gz", "A_L001_R2_001.fastq.gz"], ["A_L002_R1_001.fastq.gz", "A_L002_R2_001.fastq.gz"]])
        self.assertEqual(fetch_vcf_first.fastq_pairs(["x_1.fq.gz", "x_2.fq.gz"]), [["x_1.fq.gz", "x_2.fq.gz"]])


class IntakeTest(unittest.TestCase):
    def test_hpo_terms_with_labels(self):
        text = "Findings:\nHP:0000252 Microcephaly\nHP:0001511 | Intrauterine growth retardation\nsee HP:0000252 again\n"
        self.assertEqual(d01_intake.hpo_terms(text),
                         {"HP:0000252": "Microcephaly", "HP:0001511": "Intrauterine growth retardation"})

    def test_docx_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "doc.docx"
            with zipfile.ZipFile(p, "w") as z:
                z.writestr("word/document.xml", "<w:document><w:p><w:t>HP:0000518 Cataract</w:t></w:p></w:document>")
            self.assertIn("HP:0000518 Cataract", d01_intake.document_text(p))

    def test_private_terms_dedupe(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "terms.txt"
            p.write_text("# protected\nSYNTH001\n", encoding="utf-8")
            self.assertEqual(d01_intake.append_private_terms(["SYNTH001", "SYNTH002", "SYNTH002"], p), ["SYNTH002"])
            self.assertEqual(p.read_text(encoding="utf-8").splitlines()[-1], "SYNTH002")
            # Words without digits never become protected terms.
            self.assertEqual(d01_intake.append_private_terms(["SAMPLE", "data"], p), [])

    def test_proband_pattern_needs_label_and_digit(self):
        found = [t.rstrip("._-") for t in d01_intake.PROBAND.findall("Proband ID: SYN-0042. The proband data were reviewed.")]
        self.assertEqual(found, ["SYN-0042"])
        self.assertEqual(d01_intake.PROBAND.findall("proband is 3 years"), [])

    def test_vcf_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "s.vcf"
            rec = lambda pos, flt: "\t".join(["15", str(pos), ".", "A", "G", "50", flt, ".", "GT", "0/1"])
            p.write_text("##fileformat=VCFv4.2\n##reference=synthetic_no_chr.fasta\n##contig=<ID=15,length=1000>\n"
                         '##INFO=<ID=DP,Number=1,Type=Integer,Description="x">\n##source=Synthetic\n'
                         + "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SYN1"]) + "\n"
                         + rec(10, "PASS") + "\n" + rec(20, "LowQual") + "\n", encoding="utf-8")
            s, samples = d01_intake.vcf_summary(p)
            self.assertEqual((s["records"], s["pass_or_missing_filter"], s["n_samples"], s["chr_prefix"]), (2, 1, 1, False))
            self.assertEqual(samples, ["SYN1"])


if __name__ == "__main__":
    unittest.main()
