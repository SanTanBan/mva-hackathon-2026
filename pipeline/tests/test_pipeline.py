#!/usr/bin/env python3
"""Synthetic tests for the pipeline scripts. Every record is invented; none is patient data.

VCF lines are assembled with "\\t".join at run time, so this file contains no VCF-style line.

Run: python -m unittest discover -s pipeline/tests -t .
"""
import csv
import tempfile
import unittest
from pathlib import Path

from pipeline import acmg_points, baf_mosaic, comphet, vcf_qc
from pipeline.common import REPO, add_chr, assert_outside_repo


def tsv(*fields):
    return "\t".join(str(f) for f in fields)


HEADER = [
    "##fileformat=VCFv4.2",
    "##reference=synthetic_reference_no_chr.fasta",
    "##contig=<ID=15,length=1000000>",
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="synthetic">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="synthetic">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="synthetic">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="synthetic">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="synthetic">',
    tsv("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SYNTH1"),
]
RECORDS = [
    tsv("15", 1000, ".", "A", "G,T", 50, "PASS", "AC=1,1;DP=25", "GT:AD:PL", "1/2:10,8,7:0,1,2,3,4,5"),
    tsv("15", 2000, ".", "C", "T", 50, "PASS", "AC=1;DP=30", "GT:AD", "0/1:15,15"),
    tsv("15", 3000, ".", "CAT", "CT", 50, "LowQual", "AC=1;DP=20", "GT:AD", "0/1:12,8"),
    tsv("15", 4000, ".", "G", "A,*", 50, "PASS", "AC=1,1;DP=22", "GT:AD", "1/2:9,7,6"),
]


class VcfQcTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write_vcf(self, header=HEADER, records=RECORDS):
        path = self.dir / "in.vcf"
        path.write_text("\n".join(header + records) + "\n", encoding="utf-8")
        return path

    def test_split_trim_and_chr_prefix(self):
        summary = vcf_qc.run(self.write_vcf(), self.dir / "out.vcf")
        body = [l.split("\t") for l in (self.dir / "out.vcf").read_text(encoding="utf-8").splitlines()
                if not l.startswith("#")]
        self.assertEqual(summary["header_problems"], [])
        self.assertTrue(all(r[0] == "chr15" for r in body))
        at_1000 = [r for r in body if r[1] == "1000"]
        self.assertEqual([r[4] for r in at_1000], ["G", "T"])
        self.assertEqual([r[9] for r in at_1000], ["1/0:10,8", "0/1:10,7"])
        self.assertEqual([r[7] for r in at_1000], ["AC=1;DP=25", "AC=1;DP=25"])
        self.assertEqual(at_1000[0][8], "GT:AD")  # PL (Number=G) dropped
        deletion = [r for r in body if r[1] == "3000"][0]
        self.assertEqual((deletion[3], deletion[4]), ("CA", "C"))
        self.assertEqual(summary["counts"]["spanning_deletion_dropped"], 1)
        self.assertEqual(summary["counts"]["filtered"], 1)
        self.assertIn("##contig=<ID=chr15", (self.dir / "out.vcf").read_text(encoding="utf-8"))

    def test_header_problems(self):
        header = [h.replace("ID=15", "ID=chr15") for h in HEADER]
        summary = vcf_qc.run(self.write_vcf(header), self.dir / "out.vcf")
        self.assertTrue(any("'chr' prefix" in p for p in summary["header_problems"]))

    def test_refuses_output_inside_repo(self):
        with self.assertRaises(SystemExit):
            assert_outside_repo(REPO / "pipeline" / "leak.vcf")

    def test_add_chr(self):
        self.assertEqual([add_chr("15"), add_chr("MT"), add_chr("chrX")], ["chr15", "chrM", "chrX"])


def row(gene, pos, impact, gt, af="0.0001", pgt="", pid=""):
    return {"chrom": "chr1", "pos": str(pos), "ref": "A", "alt": "G", "gene": gene, "consequence": "x",
            "impact": impact, "pop_af": af, "gt": gt, "pgt": pgt, "pid": pid}


class CompHetTest(unittest.TestCase):
    def test_pairs_phase_and_filters(self):
        rows = [
            row("G1", 100, "HIGH", "0/1", pgt="0|1", pid="100_A_G"),
            row("G1", 200, "MODERATE", "0/1", pgt="1|0", pid="100_A_G"),
            row("G2", 300, "HIGH", "0/1"),
            row("G2", 400, "MODERATE", "0/1"),
            row("G3", 500, "HIGH", "0/1", af="0.2"),            # too common
            row("G3", 600, "HIGH", "0/1"),
            row("BUB1B", 700, "MODIFIER", "0/1"),               # kept: non-coding allowed in BUB1B
            row("BUB1B", 800, "HIGH", "0/1"),
            row("G4", 900, "HIGH", "1/1"),
        ]
        out = comphet.candidates(rows)
        by = {(r["gene"], r["model"]): r for r in out}
        self.assertEqual(by[("G1", "compound_het")]["phase"], "trans")
        self.assertEqual(by[("G2", "compound_het")]["phase"], "unknown")
        self.assertNotIn(("G3", "compound_het"), by)
        self.assertIn(("BUB1B", "compound_het"), by)
        self.assertIn(("G4", "homozygous"), by)
        self.assertEqual(out[0]["phase"], "trans")

    def test_cis_detected(self):
        rows = [row("G1", 100, "HIGH", "0/1", pgt="0|1", pid="p"), row("G1", 200, "HIGH", "0/1", pgt="0|1", pid="p")]
        self.assertEqual(comphet.candidates(rows)[0]["phase"], "cis")


class BafMosaicTest(unittest.TestCase):
    def test_recovers_mosaic_trisomy(self):
        f = 0.4
        truth = f / (2 * (2 + f))
        data = {f"chr{i}": baf_mosaic.simulate(3000, 40, 0.0, seed=i) for i in range(1, 23)}
        data["chr8"] = baf_mosaic.simulate(3000, 40, truth, seed=99)
        rows = {r["chrom"]: r for r in baf_mosaic.analyse(data)}
        self.assertAlmostEqual(rows["chr8"]["delta"], truth, delta=0.02)
        self.assertAlmostEqual(rows["chr8"]["trisomy_fraction_if_gain"], f, delta=0.12)
        self.assertEqual(max(rows.values(), key=lambda r: r["robust_z"])["chrom"], "chr8")
        self.assertLess(rows["chr1"]["delta"], 0.02)

    def test_collect_filters_sites(self):
        lines = HEADER + [
            tsv("15", 10, ".", "A", "G", 50, "PASS", ".", "GT:AD", "0/1:20,18"),   # kept
            tsv("15", 20, ".", "A", "G", 50, "PASS", ".", "GT:AD", "0/1:5,4"),     # too shallow
            tsv("15", 30, ".", "A", "G", 50, "PASS", ".", "GT:AD", "1/1:0,30"),    # not heterozygous
            tsv("15", 40, ".", "AT", "A", 50, "PASS", ".", "GT:AD", "0/1:20,20"),  # not an SNV
            tsv("15", 50, ".", "A", "G", 50, "LowQual", ".", "GT:AD", "0/1:20,20"),  # filtered
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "in.vcf"
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            got = baf_mosaic.collect(path)
        self.assertEqual(dict(got), {"15": {(18, 38): 1}})


class AcmgPointsTest(unittest.TestCase):
    def test_thresholds(self):
        self.assertEqual(acmg_points.classify(["PVS1", "PM2_Supporting", "PM3_Supporting"])["class"], "Pathogenic")
        self.assertEqual(acmg_points.classify(["PVS1", "PM2_Supporting"])["class"], "Likely pathogenic")
        self.assertEqual(acmg_points.classify(["PS3", "PM3_Supporting", "PP3"])["points"], 6)
        self.assertEqual(acmg_points.classify(["PM2_Supporting", "PP3"])["class"], "Uncertain significance")
        self.assertEqual(acmg_points.classify(["BS1", "BP4"])["class"], "Likely benign")
        self.assertEqual(acmg_points.classify(["BS1", "BS2"])["class"], "Benign")
        self.assertTrue(acmg_points.classify(["BA1"])["class"].startswith("Benign"))

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            acmg_points.classify(["PM2", "PM2_Supporting"])
        with self.assertRaises(ValueError):
            acmg_points.classify(["PX9"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
