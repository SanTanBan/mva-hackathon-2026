#!/usr/bin/env python3
"""Local lookup in the public BUBR1 stability map (package R13), with the pre-registered call.

The map comes from the Kaggle notebook kaggle/mva_bub1b_stability_map/ (public data only). The rules
for reading it were fixed before it was computed, in reports/track2/stability_map_rules.md:

- Validity gate: at least 3 of the 5 tested abundance-class substitutions (R727C, R814H, L844F,
  I909T, L1012P; Suijkerbuijk 2010) score a higher predicted ddG than both no-defect substitutions
  (R550Q, Q921H). If not, the map is not used for the child's substitution.
- Reference set: every substitution at residues 727-1050 with pLDDT >= 70.
- Call: ddG at or above the 75th percentile -> "supports"; at or below the 25th percentile at an
  exposed residue (relative accessibility > 0.25) -> "argues against"; otherwise "inconclusive".

The child's substitution is looked up here, on the laptop, so it never enters the Kaggle notebook.
Output that names it is patient-derived and goes outside the repository.

Run: python -m pipeline.stability_map MAP.tsv L1012P [--out D:\\mva-hackathon\\work\\...\\stability.json]
"""
import argparse
import csv
import json
import math
import re

from pipeline.common import assert_outside_repo

ABUNDANCE = ("R727C", "R814H", "L844F", "I909T", "L1012P")
NO_DEFECT = ("R550Q", "Q921H")
REF_FROM, MIN_PLDDT, EXPOSED = 727, 70.0, 0.25


def load(path):
    out = {}
    with open(path, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[f"{r['wt']}{r['pos']}{r['mut']}"] = {k: (float(v) if k not in ("wt", "mut") else v) for k, v in r.items()}
    return out


def _finite(x):
    return x is not None and not math.isnan(x)


def percentile(values, v):
    vals = [x for x in values if _finite(x)]
    return 100.0 * sum(x < v for x in vals) / len(vals) if vals and _finite(v) else None


def gate(m):
    nd = max(m[k]["ddg_destab"] for k in NO_DEFECT)
    above = sum(m[k]["ddg_destab"] > nd for k in ABUNDANCE)
    return {"abundance_controls_above_both_no_defect": above, "passes": _finite(nd) and above >= 3}


def call(m, name):
    if not re.match(r"^[A-Z]\d+[A-Z]$", name) or name not in m:
        raise ValueError(f"{name} is not a substitution in the map")
    ref = [r["ddg_destab"] for r in m.values() if r["pos"] >= REF_FROM and r["plddt"] >= MIN_PLDDT]
    ref_llr = [r["esm_llr"] for r in m.values() if r["pos"] >= REF_FROM and r["plddt"] >= MIN_PLDDT]
    v = m[name]
    g = gate(m)
    pct = percentile(ref, v["ddg_destab"])
    result = {"substitution": name, "ddg_destab": v["ddg_destab"], "ddg_pct_ref": pct, "esm_llr": v["esm_llr"],
              "esm_pct_ref_more_damaging": percentile([-x for x in ref_llr], -v["esm_llr"]),
              "plddt": v["plddt"], "rsa": v["rsa"], "gate": g, "reference_n": len(ref)}
    if not g["passes"]:
        result["call"] = "not interpreted (validity gate failed)"
    elif pct is None:
        result["call"] = "inconclusive (no ddG)"
    elif pct >= 75:
        result["call"] = "supports"
    elif pct <= 25 and v["rsa"] > EXPOSED:
        result["call"] = "argues against"
    else:
        result["call"] = "inconclusive"
    result["controls"] = {k: {"ddg_destab": m[k]["ddg_destab"], "ddg_pct_ref": percentile(ref, m[k]["ddg_destab"])}
                          for k in ABUNDANCE + NO_DEFECT}
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("map")
    ap.add_argument("substitution", help="one-letter form, e.g. L1012P")
    ap.add_argument("--out", help="JSON output; must be outside the repository")
    a = ap.parse_args()
    res = call(load(a.map), a.substitution)
    if a.out:
        assert_outside_repo(a.out)
        with open(a.out, "w", encoding="utf-8") as fh:
            json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
