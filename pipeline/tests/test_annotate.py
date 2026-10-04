#!/usr/bin/env python3
"""Synthetic tests for pipeline/mane_gff.py and pipeline/annotate_table.py. Nothing here is patient data.

Lines are assembled with "\\t".join at run time, so this file contains no VCF- or GFF-style line.

Run: python -m unittest discover -s pipeline/tests -t .
"""
import unittest

from pipeline import annotate_table, comphet, mane_gff


def tsv(*fields):
    return "\t".join(str(f) for f in fields) + "\n"


GFF = [
    "##gff-version 3\n",
    "##sequence-region   7 1 1000000\n",
    tsv(7, "x", "gene", 1000, 5000, ".", "+", ".", "ID=gene:G1;Name=GENEA;biotype=protein_coding"),
    tsv(7, "x", "mRNA", 1000, 5000, ".", "+", ".", "ID=transcript:T1;Parent=gene:G1;tag=gencode_basic,MANE_Select"),
    tsv(7, "x", "mRNA", 1000, 4000, ".", "+", ".", "ID=transcript:T2;Parent=gene:G1;tag=gencode_basic"),
    tsv(7, "x", "exon", 1000, 1200, ".", "+", ".", "Parent=transcript:T1;Name=E1"),
    tsv(7, "x", "CDS", 1100, 1200, ".", "+", "0", "ID=CDS:P1;Parent=transcript:T1"),
    tsv(7, "x", "exon", 1000, 1200, ".", "+", ".", "Parent=transcript:T2;Name=E2"),
    tsv(7, "x", "gene", 20000, 30000, ".", "-", ".", "ID=gene:G2;Name=GENEB;biotype=protein_coding"),
    tsv(7, "x", "mRNA", 20000, 30000, ".", "-", ".", "ID=transcript:T3;Parent=gene:G2;tag=MANE_Select"),
    tsv(7, "x", "gene", 40000, 41000, ".", "+", ".", "ID=gene:G3;Name=GENEC;biotype=protein_coding"),
    tsv(7, "x", "mRNA", 40000, 41000, ".", "+", ".", "ID=transcript:T4;Parent=gene:G3;tag=MANE_Plus_Clinical"),
]


def kept_gff():
    transcripts, genes = mane_gff.mane_ids(GFF)
    return list(mane_gff.keep_lines(GFF, transcripts, genes))


class ManeGffTest(unittest.TestCase):
    def test_keeps_only_mane_select(self):
        kept = "".join(kept_gff())
        self.assertIn("transcript:T1", kept)
        self.assertIn("transcript:T3", kept)
        self.assertNotIn("transcript:T2;", kept)          # non-MANE transcript and its exon
        self.assertNotIn("Name=E2", kept)
        self.assertNotIn("GENEC", kept)                    # MANE Plus Clinical only
        self.assertIn("##sequence-region", kept)
        self.assertIn("CDS:P1", kept)


class BcsqTest(unittest.TestCase):
    def test_worst_per_gene_and_markers(self):
        value = ",".join([
            "synonymous|GENEA|T1|protein_coding|+|5A|100C>T",
            "*missense&splice_region|GENEA|T1|protein_coding|+|6R>6C|103C>T",
            "@12345",
            "intron|GENEB|T3|protein_coding|-",
        ])
        out = annotate_table.parse_bcsq(value)
        self.assertEqual(out["GENEA"]["consequence"], "missense&splice_region")
        self.assertEqual(out["GENEA"]["impact"], "MODERATE")
        self.assertEqual(out["GENEB"]["impact"], "MODIFIER")
        self.assertEqual(annotate_table.parse_bcsq("."), {})

    def test_impact_terms(self):
        self.assertEqual(annotate_table.impact_of("frameshift&stop_retained"), "HIGH")
        self.assertEqual(annotate_table.impact_of("splice_region&intron"), "LOW")
        self.assertEqual(annotate_table.impact_of("something_new"), "MODIFIER")

    def test_short_frameshift(self):
        long_tail = "94K" + "E" * 200 + "*"
        s = annotate_table.short(long_tail + ">94KE*")
        self.assertTrue(s.endswith(">94KE*"))
        self.assertIn("...", s)
        self.assertEqual(annotate_table.short("34Q>34*"), "34Q>34*")


class GeneIndexTest(unittest.TestCase):
    def setUp(self):
        self.idx = annotate_table.GeneIndex(kept_gff())

    def test_strand_aware_flanks(self):
        self.assertEqual(self.idx.near("7", 3000, 5000), {"GENEA": "intron_variant"})
        self.assertEqual(self.idx.near("7", 500, 1000), {"GENEA": "upstream_gene_variant"})
        self.assertEqual(self.idx.near("7", 30500, 1000), {"GENEB": "upstream_gene_variant"})
        self.assertEqual(self.idx.near("7", 19500, 1000), {"GENEB": "downstream_gene_variant"})
        self.assertEqual(self.idx.near("7", 12000, 1000), {})
        self.assertEqual(self.idx.near("8", 3000, 1000), {})


class RowsTest(unittest.TestCase):
    def test_rows(self):
        idx = annotate_table.GeneIndex(kept_gff())
        info = "BCSQ=stop_gained|GENEA|T1|protein_coding|+|40Q>40*|1150C>T;CLNSIG=Pathogenic;AF_EXAC=0.0001;AF_TGP=0.0004"
        lines = [
            "#header\n",
            tsv(7, 1150, ".", "C", "T", 50, "PASS", info, "GT:AD:DP:GQ:PGT:PID", "0/1:5,5:10:99:0|1:1150_C_T"),
            tsv(7, 3000, ".", "A", "G", 50, "PASS", ".", "GT:AD", "0/1:6,6"),
            tsv(7, 12000, ".", "A", "G", 50, "PASS", ".", "GT", "1/1"),
            tsv(7, 13000, ".", "A", "G,T", 50, "PASS", ".", "GT", "1/2"),
        ]
        out = list(annotate_table.rows(lines, idx, flank=1000))
        self.assertEqual(len(out), 3)                      # the unsplit multi-allelic is skipped
        first = out[0]
        self.assertEqual((first["chrom"], first["gene"], first["impact"]), ("chr7", "GENEA", "HIGH"))
        self.assertEqual(first["pop_af"], "0.0004")
        self.assertEqual((first["pgt"], first["pid"]), ("0|1", "1150_C_T"))
        self.assertEqual(out[1]["consequence"], "intron_variant")
        self.assertEqual((out[2]["consequence"], out[2]["gene"], out[2]["pop_af"]), ("intergenic_variant", "", ""))
        self.assertEqual((first["pop_af_source"], out[2]["pop_af_source"]), ("clinvar", "none"))


class GnomadFrequencyTest(unittest.TestCase):
    """R10: gnomAD inside the fetched slices, ClinVar outside them."""

    def table(self, *infos):
        idx = annotate_table.GeneIndex(kept_gff())
        bcsq = "BCSQ=missense|GENEA|T1|protein_coding|+|40R>40C|1150C>T"
        lines = [tsv(7, 1150, ".", "C", "T", 50, "PASS", f"{bcsq};{i}" if i else bcsq, "GT", "0/1")
                 for i in infos]
        return list(annotate_table.rows(lines, idx, flank=1000))

    def test_gnomad_common_variant_is_dropped(self):
        row = self.table("GNOMAD_COVERED;AF_joint=0.04;nhomalt_joint=120;fafmax_faf95_max_joint=0.05;AF_EXAC=0.0001")[0]
        self.assertEqual((row["pop_af"], row["pop_af_source"], row["nhomalt_joint"]), ("0.05", "gnomad_v4.1", "120"))
        self.assertFalse(comphet.keep(row, 0.01, {"HIGH", "MODERATE"}, set()))

    def test_absent_covered_variant_is_rare(self):
        row = self.table("GNOMAD_COVERED")[0]
        self.assertEqual((row["pop_af"], row["pop_af_source"]), ("0", "gnomad_v4.1"))
        self.assertTrue(comphet.keep(row, 0.01, {"HIGH", "MODERATE"}, set()))

    def test_uncovered_variant_falls_back_to_clinvar(self):
        row = self.table("AF_EXAC=0.02;AF_TGP=0.03")[0]
        self.assertEqual((row["pop_af"], row["pop_af_source"]), ("0.03", "clinvar"))
        self.assertFalse(comphet.keep(row, 0.01, {"HIGH", "MODERATE"}, set()))

    def test_af_joint_used_when_no_faf(self):
        row = self.table("GNOMAD_COVERED;AF_joint=0.002")[0]
        self.assertEqual((row["pop_af"], row["pop_af_source"]), ("0.002", "gnomad_v4.1"))


if __name__ == "__main__":
    unittest.main()
