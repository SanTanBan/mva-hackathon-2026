#!/usr/bin/env python3
"""Validate a Track 1 predictions CSV against the MVA Hackathon 2026 spec.

The spec comes from the "Submit - Track 1" tab (template updated 28 Aug 2026).
Each row is one proposed causal variant or compound-het pair, with at most 10
rows, ranked by EPCR (estimated probability of causal relationship), in GRCh38
coordinates.

Usage:
    python tools/validate_track1_csv.py predictions.csv [--hf-user NAME] [--template template.csv]

Exits with 0 when the file is valid and 1 when errors are found. Warnings do not fail it.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

COLUMNS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
           "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]
SECOND = ["chrom_2", "pos_2", "ref_2", "alt_2"]
CHROM = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$")
BARE_CHROM = re.compile(r"^([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$")
ALLELE = re.compile(r"^[ACGTN]+$")
FINDING_TYPES = {"primary", "secondary"}
MAX_ROWS = 10
# The Space's scorer rejects any other ID ("Only PROBAND01 is accepted in this challenge", 2026-10-23),
# although the VCF sample name differs.
ACCEPTED_PROBAND = "PROBAND01"


def check_variant(row, n, suffix, errors):
    chrom, pos, ref, alt = (row[f"{k}_{suffix}"].strip() for k in ("chrom", "pos", "ref", "alt"))
    where = f"row {n}, variant {suffix}"
    if not CHROM.match(chrom):
        hint = ""
        if BARE_CHROM.match(chrom):
            hint = " (the VCF has no 'chr' prefix, but the template expects e.g. chr15)"
        errors.append(f"{where}: chrom '{chrom}' is not chr1-22/X/Y/M{hint}")
    if not pos.isdigit() or int(pos) < 1:
        errors.append(f"{where}: pos '{pos}' is not a positive integer")
    for name, allele in (("ref", ref), ("alt", alt)):
        if not ALLELE.match(allele):
            errors.append(f"{where}: {name} '{allele}' must be uppercase A/C/G/T/N")
    if ref and ref == alt:
        errors.append(f"{where}: ref and alt are identical")
    return (chrom, pos, ref, alt)


def validate(path, hf_user=None, template=None):
    errors, warnings = [], []
    raw = path.read_bytes()
    if raw.startswith(b"\xef\xbb\xbf"):
        warnings.append("file starts with a UTF-8 BOM (Excel adds one); save as plain UTF-8 CSV")
    text = raw.decode("utf-8-sig")

    expected = COLUMNS
    if template:
        with open(template, newline="", encoding="utf-8-sig") as fh:
            expected = [c.strip() for c in next(csv.reader(fh))]
        if expected != COLUMNS:
            warnings.append(f"official template header differs from the documented spec: {expected}")

    reader = csv.DictReader(text.splitlines())
    header = [c.strip() for c in (reader.fieldnames or [])]
    if header != expected:
        missing = [c for c in expected if c not in header]
        extra = [c for c in header if c not in expected]
        detail = []
        if missing:
            detail.append(f"missing {missing}")
        if extra:
            detail.append(f"unexpected {extra}")
        if not missing and not extra:
            detail.append("columns are out of order")
        errors.append(f"header must be exactly {expected}: " + "; ".join(detail))
        return errors, warnings

    rows = [r for r in reader if any((v or "").strip() for v in r.values())]
    if not rows:
        errors.append("no prediction rows")
        return errors, warnings
    if len(rows) > MAX_ROWS:
        errors.append(f"{len(rows)} rows; at most {MAX_ROWS} are accepted")

    probands = {r["proband_id"].strip() for r in rows}
    if "" in probands:
        errors.append("proband_id is blank on at least one row")
    if len(probands - {""}) > 1:
        errors.append(f"rows use more than one proband_id: {sorted(probands - {''})}")
    if probands - {"", ACCEPTED_PROBAND}:
        errors.append(f"proband_id must be {ACCEPTED_PROBAND}; the scorer rejects any other value")

    seen, epcrs = set(), []
    for n, row in enumerate(rows, start=1):
        first = check_variant(row, n, "1", errors)
        filled = [bool(row[c].strip()) for c in SECOND]
        key = (first,)
        if all(filled):
            second = check_variant(row, n, "2", errors)
            if second == first:
                errors.append(f"row {n}: variant 2 repeats variant 1")
            elif second[0] != first[0]:
                warnings.append(f"row {n}: the pair spans two chromosomes; a compound het is normally in one gene")
            key = tuple(sorted([first, second]))
        elif any(filled):
            errors.append(f"row {n}: chrom_2/pos_2/ref_2/alt_2 must be all filled (a pair) or all blank")
        if key in seen:
            errors.append(f"row {n}: duplicates an earlier row")
        seen.add(key)

        try:
            epcr = float(row["epcr"])
            if not 0 < epcr <= 1:
                errors.append(f"row {n}: epcr {epcr} must be in (0, 1]")
            epcrs.append(epcr)
        except ValueError:
            errors.append(f"row {n}: epcr '{row['epcr']}' is not a number")

        ftype = row["finding_type"].strip()
        if ftype not in FINDING_TYPES:
            errors.append(f"row {n}: finding_type '{ftype}' must be 'primary' or 'secondary' (lowercase)")
        elif ftype == "secondary" and not row["notes"].strip():
            warnings.append(f"row {n}: secondary finding with no notes; the organizers ask you to explain why it is flagged")

    if not any(r["finding_type"].strip() == "primary" for r in rows):
        errors.append("no row has finding_type 'primary'")
    if len(epcrs) == len(rows):
        if any(a < b for a, b in zip(epcrs, epcrs[1:])):
            errors.append("rows are not sorted by epcr, highest first; rank order must match EPCR")
        if len(set(epcrs)) < len(epcrs):
            warnings.append("tied epcr values; the tie-break order is unspecified")

    if hf_user and hf_user.lower() not in path.name.lower():
        warnings.append(f"filename should include the HF username, e.g. {hf_user}_<approach>.csv")
    return errors, warnings


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("csv", type=Path)
    ap.add_argument("--hf-user")
    ap.add_argument("--template", type=Path, help="official submission template CSV, once downloaded")
    args = ap.parse_args()
    if not args.csv.is_file():
        sys.exit(f"not found: {args.csv}")
    errors, warnings = validate(args.csv, args.hf_user, args.template)
    for w in warnings:
        print(f"WARNING  {w}")
    for e in errors:
        print(f"ERROR    {e}")
    print(f"{'INVALID' if errors else 'VALID'}: {args.csv.name} "
          f"({len(errors)} error(s), {len(warnings)} warning(s))")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
