#!/usr/bin/env python3
"""Local lookup in the public BUB1B splice map (package R09), for D04.

The map scores every possible single-nucleotide change across the BUB1B MANE Select transcript
(NM_001211.6) with Pangolin. It was made on Kaggle from public data only
(kaggle/bub1b_splice_map/), and lives at D:\\mva-hackathon\\refs\\bub1b\\kaggle_out\\bub1b_pangolin_snv.tsv.gz.
Looking up the child's variants here keeps them on this machine. Indels are not in the map; they
need an individual Pangolin run in the same way.

Each hit is also placed in the distance classes of the D08 worksheet (Zardetto 2024 criteria), so the
splice-correction lane in pipeline/allele_mechanism.py can be given its conditional tier.

Run: python -m pipeline.splice_map c.2386-11A>G
     python -m pipeline.splice_map 15:40212488:A:G
"""
import argparse
import csv
import gzip
import json
import re
from pathlib import Path

DEFAULT_MAP = Path(r"D:\mva-hackathon\refs\bub1b\kaggle_out\bub1b_pangolin_snv.tsv.gz")
GAIN_STRONG, GAIN_WEAK = 0.5, 0.2


def load(path=DEFAULT_MAP):
    """Index the map by HGVS c. name and by (chrom without prefix, pos, ref, alt)."""
    by_hgvs, by_pos = {}, {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in ("gain", "loss"):
                r[k] = float(r[k])
            for k in ("pos", "gain_offset", "loss_offset"):
                r[k] = int(r[k])
            by_hgvs[r["hgvs_c"]] = r
            by_pos[(r["chrom"].removeprefix("chr"), r["pos"], r["ref"], r["alt"])] = r
    return by_hgvs, by_pos


def distance_class(hgvs_c, region):
    """The worksheet's distance classes for splice correction."""
    if region == "exon":
        return "exonic"
    m = re.match(r"^c\.[-*]?\d+([+-])(\d+)", hgvs_c)
    side, d = ("donor" if m.group(1) == "+" else "acceptor"), int(m.group(2))
    if d <= 2:
        return "canonical"
    if side == "acceptor" and 18 <= d <= 40:
        return "branch-point window"
    if d < 5:
        return "under 5 nt"
    if d < 15:
        return "5-14 nt"
    if d <= 100:
        return "15-100 nt"
    return "deep intronic"


def describe(row):
    cls = distance_class(row["hgvs_c"], row["region"])
    gain = "strong" if row["gain"] >= GAIN_STRONG else "weak" if row["gain"] >= GAIN_WEAK else "none"
    loss = "strong" if row["loss"] <= -GAIN_STRONG else "weak" if row["loss"] <= -GAIN_WEAK else "none"
    lane = None
    if gain == "strong" and cls in ("15-100 nt", "deep intronic"):
        lane = "candidate for splice correction if RNA confirms the cryptic site"
    elif gain == "strong" and cls in ("5-14 nt", "under 5 nt", "branch-point window"):
        lane = "cryptic site too close to the canonical site for splice correction"
    return {"hgvs_c": row["hgvs_c"], "region": row["region"], "distance_class": cls,
            "gain": row["gain"], "gain_offset": row["gain_offset"], "gain_call": gain,
            "loss": row["loss"], "loss_offset": row["loss_offset"], "loss_call": loss, "note": lane}


def lookup(query, index):
    by_hgvs, by_pos = index
    if query.startswith("c."):
        row = by_hgvs.get(query)
    else:
        chrom, pos, ref, alt = query.split(":")
        row = by_pos.get((chrom.removeprefix("chr"), int(pos), ref.upper(), alt.upper()))
    return describe(row) if row else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("queries", nargs="+", help="HGVS c. names on NM_001211.6, or CHROM:POS:REF:ALT")
    ap.add_argument("--map", default=str(DEFAULT_MAP))
    args = ap.parse_args()
    index = load(args.map)
    for q in args.queries:
        print(json.dumps({"query": q, "result": lookup(q, index)}))


if __name__ == "__main__":
    main()
