#!/usr/bin/env python3
"""ACMG/AMP classification with the naturally scaled point system.

Evidence strengths score Supporting 1, Moderate 2, Strong 4 and Very Strong 8 points;
benign evidence scores negative. Totals map to classes:
    Pathogenic >= 10, Likely pathogenic 6 to 9, VUS 0 to 5, Likely benign -1 to -6, Benign <= -7
(Tavtigian et al., Hum Mutat 2020;41:1734-1737, PMID 32720330). BA1 is stand-alone benign.

Codes use their default ACMG/AMP strength, or a ClinGen-style modifier, e.g. PVS1_Strong,
PM2_Supporting, PP3_Moderate or PM3_Supporting. The per-allele arithmetic is returned so the
methods report can show it.

Run: python -m pipeline.acmg_points PVS1 PM2_Supporting PM3_Supporting
"""
import argparse
import re

POINTS = {"Supporting": 1, "Moderate": 2, "Strong": 4, "VeryStrong": 8}
DEFAULT = {"PVS1": "VeryStrong"}
DEFAULT.update({f"PS{i}": "Strong" for i in range(1, 5)})
DEFAULT.update({f"PM{i}": "Moderate" for i in range(1, 7)})
DEFAULT.update({f"PP{i}": "Supporting" for i in range(1, 6)})
DEFAULT.update({f"BS{i}": "Strong" for i in range(1, 5)})
DEFAULT.update({f"BP{i}": "Supporting" for i in range(1, 8)})
CODE = re.compile(r"^(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(?:_(Supporting|Moderate|Strong|VeryStrong))?$")


def score(code):
    m = CODE.match(code.strip().replace("Very_Strong", "VeryStrong"))
    if not m:
        raise ValueError(f"unrecognised evidence code: {code}")
    base, strength = m.group(1), m.group(2)
    if base == "BA1":
        return base, "StandAlone", None
    strength = strength or DEFAULT[base]
    pts = POINTS[strength]
    return base, strength, -pts if base.startswith("B") else pts


def classify(codes):
    rows = [score(c) for c in codes]
    bases = [r[0] for r in rows]
    if len(set(bases)) != len(bases):
        raise ValueError("each criterion may be applied once per variant")
    if "BA1" in bases:
        return {"points": None, "class": "Benign (BA1 stand-alone)", "evidence": rows}
    total = sum(r[2] for r in rows)
    if total >= 10:
        cls = "Pathogenic"
    elif total >= 6:
        cls = "Likely pathogenic"
    elif total >= 0:
        cls = "Uncertain significance"
    elif total >= -6:
        cls = "Likely benign"
    else:
        cls = "Benign"
    return {"points": total, "class": cls, "evidence": rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("codes", nargs="+")
    result = classify(ap.parse_args().codes)
    for base, strength, pts in result["evidence"]:
        print(f"  {base:<5} {strength:<11} {'' if pts is None else f'{pts:+d}'}")
    print(f"total {result['points']} -> {result['class']}")


if __name__ == "__main__":
    main()
