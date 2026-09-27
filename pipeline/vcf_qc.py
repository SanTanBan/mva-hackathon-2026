#!/usr/bin/env python3
"""VCF quality control and normalisation.

Streams a VCF (plain or bgzip) one record at a time, so memory stays flat. It:
  1. checks the header: file format, reference name, contig naming (the hackathon VCF uses
     contigs without a 'chr' prefix) and the number of samples;
  2. splits multi-allelic records into one record per ALT. GT, PGT and AD are recoded, and
     other ALT alleles become reference (0), so a 1/2 genotype becomes 1/0 and 0/1: two
     heterozygous records at one position. Number=A and Number=R INFO values are subset;
     Number=G values are dropped;
  3. trims bases shared by REF and ALT (right end, then left end, keeping one anchor base);
  4. writes records with 'chr'-prefixed contigs, as the Track 1 CSV requires.

Only these FORMAT fields are kept: GT, AD, DP, GQ, PGT, PID, PS. Spanning deletions ('*')
are dropped. Left-aligning indels needs the reference FASTA, so it is left to
`bcftools norm -f` once the toolchain is installed; the summary counts the indels.

Run: python -m pipeline.vcf_qc IN.vcf.gz OUT.vcf --summary OUT.json   (outputs outside the repo)
"""
import argparse
import json
import re
from collections import Counter

from pipeline.common import add_chr, assert_outside_repo, open_text

KEEP_FORMAT = ("GT", "AD", "DP", "GQ", "PGT", "PID", "PS")
INFO_DEF = re.compile(r"^##INFO=<ID=([^,>]+),Number=([^,>]+)")
CONTIG_DEF = re.compile(r"^##contig=<ID=([^,>]+)")


def read_header(lines):
    """Consume header lines; return (meta_lines, column_line, header_info)."""
    meta, info = [], {"fileformat": None, "reference": None, "contigs": [], "info_number": {}, "samples": []}
    for line in lines:
        line = line.rstrip("\r\n")
        if line.startswith("##"):
            meta.append(line)
            if line.startswith("##fileformat="):
                info["fileformat"] = line.split("=", 1)[1]
            elif line.startswith("##reference="):
                info["reference"] = line.split("=", 1)[1]
            elif m := CONTIG_DEF.match(line):
                info["contigs"].append(m.group(1))
            elif m := INFO_DEF.match(line):
                info["info_number"][m.group(1)] = m.group(2)
        elif line.startswith("#CHROM"):
            info["samples"] = line.split("\t")[9:]
            return meta, line, info
    raise SystemExit("No #CHROM header line found")


def check_header(info, expect_chr=False, expect_samples=1):
    problems = []
    if not (info["fileformat"] or "").startswith("VCFv4"):
        problems.append(f"unexpected fileformat: {info['fileformat']}")
    if not info["reference"]:
        problems.append("no ##reference line")
    prefixed = [c for c in info["contigs"] if c.startswith("chr")]
    if not expect_chr and prefixed:
        problems.append(f"{len(prefixed)} contigs carry a 'chr' prefix but the dataset is documented without it")
    if expect_chr and len(prefixed) < len(info["contigs"]):
        problems.append("some contigs lack the expected 'chr' prefix")
    if len(info["samples"]) != expect_samples:
        problems.append(f"expected {expect_samples} sample(s), found {len(info['samples'])}")
    return problems


def trim(pos, ref, alt):
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def recode_gt(gt, keep):
    """Recode allele indices for the split record carrying ALT number `keep`."""
    if not gt or gt in (".", "./.", ".|."):
        return gt
    sep = "|" if "|" in gt else "/"
    out = []
    for allele in gt.split(sep):
        out.append(allele if allele == "." else ("1" if allele == str(keep) else "0"))
    return sep.join(out)


def split_info(info_field, keep, n_alt, numbers):
    if info_field in (".", "") or n_alt == 1:
        return info_field
    out = []
    for item in info_field.split(";"):
        key, _, val = item.partition("=")
        number = numbers.get(key)
        if not val or number not in ("A", "R", "G"):
            out.append(item)
        elif number == "A":
            vals = val.split(",")
            out.append(f"{key}={vals[keep - 1]}" if len(vals) == n_alt else item)
        elif number == "R":
            vals = val.split(",")
            out.append(f"{key}={vals[0]},{vals[keep]}" if len(vals) == n_alt + 1 else item)
    return ";".join(out) if out else "."


def normalise_record(fields, numbers, counts):
    """Yield normalised output field lists for one input record."""
    chrom, pos, vid, ref, alts, qual, filt, info = fields[:8]
    fmt_keys = fields[8].split(":") if len(fields) > 8 else []
    samples = [dict(zip(fmt_keys, s.split(":"))) for s in fields[9:]]
    alt_list = alts.split(",")
    if len(alt_list) > 1:
        counts["multiallelic_split"] += 1
    for i, alt in enumerate(alt_list, 1):
        if alt == "*":
            counts["spanning_deletion_dropped"] += 1
            continue
        p, r, a = trim(int(pos), ref, alt)
        kind = "snv" if len(r) == len(a) == 1 else ("mnv" if len(r) == len(a) else "indel")
        counts[kind] += 1
        keys = [k for k in KEEP_FORMAT if k in fmt_keys]
        sample_out = []
        for s in samples:
            vals = []
            for k in keys:
                v = s.get(k, ".")
                if k in ("GT", "PGT"):
                    v = recode_gt(v, i)
                elif k == "AD" and v != ".":
                    ad = v.split(",")
                    v = f"{ad[0]},{ad[i]}" if len(ad) == len(alt_list) + 1 else "."
                vals.append(v)
            sample_out.append(":".join(vals))
        yield [add_chr(chrom), str(p), vid, r, a, qual, filt, split_info(info, i, len(alt_list), numbers),
               ":".join(keys)] + sample_out


def run(in_path, out_path, summary_path=None, expect_chr=False, expect_samples=1):
    assert_outside_repo(out_path)
    counts = Counter()
    contigs_seen = Counter()
    with open_text(in_path) as fh, open(out_path, "w", encoding="utf-8", newline="\n") as out:
        meta, column_line, info = read_header(fh)
        problems = check_header(info, expect_chr, expect_samples)
        for line in meta:
            if m := CONTIG_DEF.match(line):
                line = line.replace(f"ID={m.group(1)}", f"ID={add_chr(m.group(1))}", 1)
            out.write(line + "\n")
        out.write("##pipeline_vcf_qc=split multi-allelic; trimmed; chr-prefixed; FORMAT kept: "
                  + ",".join(KEEP_FORMAT) + "\n")
        out.write(column_line + "\n")
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            counts["records_in"] += 1
            counts["pass" if fields[6] in ("PASS", ".") else "filtered"] += 1
            contigs_seen[fields[0]] += 1
            for rec in normalise_record(fields, info["info_number"], counts):
                out.write("\t".join(rec) + "\n")
                counts["records_out"] += 1
    summary = {"header": {k: info[k] for k in ("fileformat", "reference")} | {
                   "n_contigs": len(info["contigs"]), "n_samples": len(info["samples"])},
               "header_problems": problems, "counts": dict(counts),
               "contigs_with_records": len(contigs_seen)}
    if summary_path:
        assert_outside_repo(summary_path)
        with open(summary_path, "w", encoding="utf-8") as fh:
            json.dump(summary, fh, indent=2)
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("vcf")
    ap.add_argument("out")
    ap.add_argument("--summary")
    ap.add_argument("--expect-chr", action="store_true", help="input contigs should carry 'chr'")
    ap.add_argument("--samples", type=int, default=1)
    args = ap.parse_args()
    summary = run(args.vcf, args.out, args.summary, args.expect_chr, args.samples)
    print(json.dumps({"header_problems": summary["header_problems"], "counts": summary["counts"]}, indent=2))


if __name__ == "__main__":
    main()
