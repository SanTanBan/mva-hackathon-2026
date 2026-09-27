#!/usr/bin/env python3
"""Tests for validate_track1_csv.py. Every row is synthetic; none is patient data.

Run: python tools/test_validate_track1_csv.py
"""
import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from validate_track1_csv import COLUMNS, validate  # noqa: E402


def pair(epcr=0.9, **over):
    row = dict(proband_id="P1", chrom_1="chr15", pos_1="1000001", ref_1="A", alt_1="G",
               chrom_2="chr15", pos_2="1000500", ref_2="C", alt_2="T",
               epcr=str(epcr), finding_type="primary", notes="")
    row.update(over)
    return row


def single(epcr=0.5, **over):
    row = dict(proband_id="P1", chrom_1="chr2", pos_1="2000001", ref_1="G", alt_1="A",
               chrom_2="", pos_2="", ref_2="", alt_2="",
               epcr=str(epcr), finding_type="secondary", notes="actionable, unrelated to phenotype")
    row.update(over)
    return row


class Track1CsvTest(unittest.TestCase):
    def run_rows(self, rows, header=COLUMNS, name="SantanuBanerjee_synthetic.csv", bom=False):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / name
            with open(path, "w", newline="", encoding="utf-8-sig" if bom else "utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)
            return validate(path, hf_user="SantanuBanerjee")

    def assert_error(self, rows, fragment, **kw):
        errors, _ = self.run_rows(rows, **kw)
        self.assertTrue(any(fragment in e for e in errors), f"expected '{fragment}' in {errors}")

    def test_valid_pair_plus_secondary(self):
        errors, warnings = self.run_rows([pair(), single()])
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])

    def test_bare_chromosome_gets_chr_hint(self):
        self.assert_error([pair(chrom_1="15")], "no 'chr' prefix")

    def test_unsorted_epcr(self):
        self.assert_error([single(epcr=0.2, finding_type="primary", notes=""), pair(epcr=0.9)], "not sorted by epcr")

    def test_half_filled_second_variant(self):
        self.assert_error([pair(ref_2="", alt_2="")], "all filled (a pair) or all blank")

    def test_finding_type_is_lowercase(self):
        self.assert_error([pair(finding_type="Primary")], "must be 'primary' or 'secondary'")

    def test_needs_a_primary(self):
        self.assert_error([single()], "no row has finding_type 'primary'")

    def test_row_limit(self):
        rows = [pair(epcr=round(1 - i * 0.05, 2), pos_1=str(1000001 + i)) for i in range(11)]
        self.assert_error(rows, "at most 10")

    def test_header_order(self):
        swapped = COLUMNS[:]
        swapped[1], swapped[2] = swapped[2], swapped[1]
        self.assert_error([pair()], "columns are out of order", header=swapped)

    def test_epcr_bounds(self):
        self.assert_error([pair(epcr=0)], "must be in (0, 1]")
        self.assert_error([pair(epcr=1.5)], "must be in (0, 1]")

    def test_duplicate_pair_in_either_order(self):
        flipped = pair(epcr=0.8, chrom_1="chr15", pos_1="1000500", ref_1="C", alt_1="T",
                       chrom_2="chr15", pos_2="1000001", ref_2="A", alt_2="G")
        self.assert_error([pair(), flipped], "duplicates an earlier row")

    def test_bom_is_only_a_warning(self):
        errors, warnings = self.run_rows([pair()], bom=True)
        self.assertEqual(errors, [])
        self.assertTrue(any("BOM" in w for w in warnings))

    def test_filename_should_carry_username(self):
        _, warnings = self.run_rows([pair()], name="predictions.csv")
        self.assertTrue(any("HF username" in w for w in warnings))


if __name__ == "__main__":
    unittest.main(verbosity=2)
