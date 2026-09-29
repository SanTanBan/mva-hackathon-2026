#!/usr/bin/env python3
"""Pre-registered BUB1B allele classification for package D08: mechanism class and antisense lane.

Written on 2026-09-29, before any access to the child's genome, so the rules are fixed before the
alleles are seen (worksheet: reports/track2/d08_worksheet.md). Input is one allele as annotated by
package D04: an HGVS c. name on NM_001211.6, a consequence term, and optionally the HGVS p. name,
the premature stop position, and RNA evidence. For each allele it reports:

- NMD. A premature stop more than 50 nt upstream of the last exon junction triggers NMD; stops
  in the last exon or the last 50 nt of the penultimate exon escape; stops in the first 100
  coding nt are rarely targeted and those up to 200 nt only partly (Lindeboom 2016, PMID 27618451).
- Mechanism class (Suijkerbuijk 2010, PMID 20516114). "Abundance": no transcript, or a
  substitution in or near the kinase domain that speeds protein turnover. "Functional": a
  stable protein with a specific defect. The eight substitutions that paper tested are looked
  up; any other is provisional by position and needs the BUBR1 immunoblot (report §7, assay 1).
- A tier for each antisense (ASO) approach, in the five classes of the N=1 Collaborative's
  VARIANT guidelines (Cheerie 2025, PMID 40139194), applying the splice-switching and
  exon-skipping criteria of Zardetto 2024 (PMID 40225926). Upregulation is reported as a lane
  without a tier, because version 1.0 of the guidelines defines none for it.

"eligible" requires an ASO already shown to work for the variant or exon. PubMed lists none for
BUB1B (searched 2026-09-29), so "likely eligible" is the best tier this script can give.

Only public reference data are embedded: the exon table is from the RefSeq record NM_001211.6
(NCBI, 2026-09-29), with intron lengths from the locus dossier. Patient-derived input and output
stay under D:\\mva-hackathon\\ (see CLAUDE.md).

Run: python -m pipeline.allele_mechanism "c.2210T>G" stop_gained --hgvs-p p.Leu737Ter
     python -m pipeline.allele_mechanism --tsv D:\\mva-hackathon\\work\\d04_alleles.tsv --out D:\\mva-hackathon\\work\\d08_alleles.json
     python -m pipeline.allele_mechanism --exon-table
"""
import argparse
import csv
import json
import re

from pipeline.common import assert_outside_repo

TRANSCRIPT = "NM_001211.6"
# (exon, first coding c., last coding c., length of the following intron). Exon 1 also holds the
# 152-nt 5' UTR and exon 23 the stop codon and 3' UTR.
EXONS = [
    (1, 1, 35, 3797), (2, 36, 179, 4865), (3, 180, 239, 415), (4, 240, 384, 5795),
    (5, 385, 581, 7040), (6, 582, 751, 1281), (7, 752, 966, 171), (8, 967, 1058, 10902),
    (9, 1059, 1288, 2840), (10, 1289, 1401, 516), (11, 1402, 1517, 571), (12, 1518, 1567, 1424),
    (13, 1568, 1628, 123), (14, 1629, 1734, 3489), (15, 1735, 2009, 2178), (16, 2010, 2143, 864),
    (17, 2144, 2284, 334), (18, 2285, 2385, 2288), (19, 2386, 2535, 683), (20, 2536, 2678, 4021),
    (21, 2679, 2850, 788), (22, 2851, 2957, 2001), (23, 2958, 3153, None),
]
CDS_LEN = 3153
LAST_JUNCTION = 2957          # last coding nt before the final exon junction
NMD_BOUNDARY = LAST_JUNCTION - 50
KINASE_DOMAIN = (766, 1050)   # UniProt O60566 kinase-like domain

# Substitutions tested by Suijkerbuijk 2010 (PMID 20516114), by protein change.
TESTED = {
    "Y155C": ("functional", "stable; restores chromosome alignment but not checkpoint activity"),
    "R550Q": ("no defect detected", "stable; rescued checkpoint and alignment"),
    "R727C": ("abundance", "unstable; faster turnover in patient cells"),
    "R814H": ("abundance", "unstable in cells (kinase-domain substitution)"),
    "L844F": ("abundance", "poorly expressed; failed to rescue"),
    "I909T": ("abundance", "unstable and HSP90-dependent; fully functional when overexpressed"),
    "Q921H": ("no defect detected", "indistinguishable from wild type; passenger on the L844F allele"),
    "L1012P": ("abundance", "unstable and HSP90-dependent; fully functional when overexpressed"),
}
# Every tested substitution at or after residue 727 was unstable; those at 155 and 550 were stable.
STABLE_UP_TO, UNSTABLE_FROM = 550, 727

# Pre-registered verdicts for the in-frame internal exons (the only single-exon skips that keep
# the reading frame). Skipping any of them creates no stop codon at the new junction.
IN_FRAME_SKIP = {
    2: ("not eligible", "removes residues 12-60, including the first KEN box (residues 26-28) that "
        "BUBR1 needs to form the mitotic checkpoint complex (Di Fiore 2015, PMID 25669885); the "
        "exon 2 acceptor variant c.36-1G>A is likely pathogenic in ClinVar"),
    3: ("unlikely eligible", "removes residues 60-80 at the start of the TPR domain (62-226), which "
        "checkpoint signalling needs (Suijkerbuijk 2010); the exon 3 donor variant c.239+2T>C is "
        "likely pathogenic in ClinVar"),
    17: ("unlikely eligible", "in-frame, 47 residues (715-762, 4.5% of the protein), with no annotated "
         "domain or KEN box and no pathogenic splice-site variant at its boundaries; residues 1-715 "
         "alone supported checkpoint and attachment in cell lines (731X, Suijkerbuijk 2010). But the "
         "MVA substitution R727C in this exon destabilises BUBR1 (Suijkerbuijk 2010), so the segment "
         "probably supports kinase-domain folding, and skipping it may destabilise the protein too "
         "(Zardetto 2024). Still the only in-frame exon clear of the KEN boxes and the TPR and kinase "
         "domains, so a stability test of the shortened protein would settle it cheaply"),
    19: ("not eligible", "lies inside the kinase-like domain, whose integrity keeps BUBR1 stable "
         "(Suijkerbuijk 2010); R814H in this exon lowers BUBR1 abundance, and the exon 19 acceptor "
         "variant c.2386-2A>G is likely pathogenic in ClinVar"),
}

TRUNCATING = {"stop_gained", "frameshift_variant", "start_lost"}
MISSENSE = {"missense_variant"}
INFRAME = {"inframe_deletion", "inframe_insertion", "protein_altering_variant"}
SYNONYMOUS = {"synonymous_variant"}
RNA_EVIDENCE = {"none", "cryptic_site", "pseudoexon", "exon_skip", "intron_retention"}
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
       "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
       "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}


def exon_of(c):
    for exon, start, end, _ in EXONS:
        if start <= c <= end:
            return exon
    raise ValueError(f"c.{c} is outside the coding sequence of {TRANSCRIPT}")


def parse_c(hgvs):
    """Locate an HGVS c. name: exonic, intronic (with side and offset), 5' or 3' of the CDS."""
    m = re.match(r"^c\.(-?\*?\d+)([+-]\d+)?", hgvs.strip())
    if not m:
        raise ValueError(f"not an HGVS c. name: {hgvs}")
    base, offset = m.group(1), m.group(2)
    if base.startswith("*"):
        return {"region": "3' UTR or downstream"}
    pos = int(base)
    if pos < 1:
        return {"region": "5' UTR or upstream" if not offset else "5' UTR intron"}
    if not offset:
        return {"region": "exon", "c": pos, "exon": exon_of(pos)}
    off = int(offset)
    for i, (exon, start, end, intron) in enumerate(EXONS):
        if off > 0 and pos == end and intron:
            return {"region": "intron", "intron": exon, "side": "donor", "offset": off,
                    "distance": min(off, intron - off + 1)}
        if off < 0 and pos == start and i > 0:
            intron = EXONS[i - 1][3]
            return {"region": "intron", "intron": exon - 1, "side": "acceptor", "offset": -off,
                    "distance": min(-off, intron + off + 1)}
    raise ValueError(f"{hgvs}: c.{pos} is not an exon boundary of {TRANSCRIPT}")


def protein_change(hgvs_p):
    """'p.Ile909Thr' or 'I909T' -> ('I', 909, 'T'); 'p.Leu737Ter' -> ('L', 737, '*'); fs -> ('fs', first, ter)."""
    if not hgvs_p:
        return None
    s = hgvs_p.strip().removeprefix("p.").strip("()")
    fs = re.match(r"^([A-Z][a-z]{2}|[A-Z])(\d+)(?:[A-Z][a-z]{2}|[A-Z])?fs(?:\*|Ter)?(\d+)?", s)
    if fs:
        first = int(fs.group(2))
        return ("fs", first, first + int(fs.group(3)) - 1 if fs.group(3) else None)
    m = re.match(r"^([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z]|Ter|\*)$", s)
    if not m:
        return None
    one = lambda a: "*" if a in ("Ter", "*") else AA3.get(a, a)
    return (one(m.group(1)), int(m.group(2)), one(m.group(3)))


def nmd(ptc_c):
    if ptc_c is None:
        return "unknown"
    if ptc_c >= EXONS[-1][1]:
        return "escapes (last exon)"
    if ptc_c > NMD_BOUNDARY:
        return "escapes (within 50 nt of the last exon junction)"
    if ptc_c <= 100:
        return "likely escapes (first 100 coding nt)"
    if ptc_c <= 200:
        return "partial (101-200 nt from the start codon)"
    return "predicted"


def mechanism(consequence, loc, change, nmd_call):
    if consequence in TRUNCATING:
        if nmd_call == "predicted":
            return "abundance", "no transcript expected (NMD)"
        return "unclassified", "the truncated protein may be made; test BUBR1 protein level (assay 1)"
    if consequence in MISSENSE and change:
        key = f"{change[0]}{change[1]}{change[2]}"
        if key in TESTED:
            return TESTED[key][0], TESTED[key][1] + " (Suijkerbuijk 2010)"
        if change[1] >= UNSTABLE_FROM:
            return "abundance (provisional)", "in or near the kinase domain, where every tested substitution was unstable"
        if change[1] <= STABLE_UP_TO:
            return "functional or no defect (provisional)", "N-terminal to residue 551, where tested substitutions were stable"
        return "unclassified", "between residues 551 and 726, where no substitution has been tested"
    if loc["region"] == "intron" or consequence in {"splice_donor_variant", "splice_acceptor_variant"}:
        return "abundance (provisional)", "less normal transcript expected; RNA needed to confirm"
    if loc["region"].startswith("5'"):
        return "abundance (provisional)", "reduced expression, as for the published upstream allele (Ochiai 2013)"
    return "unclassified", "no rule applies; assess by hand"


def splice_correction(loc, consequence, rna):
    """Tier for blocking a cryptic splice site so that normal splicing returns."""
    if loc["region"] == "intron":
        d, side = loc["distance"], loc["side"]
        if d <= 2:
            return "not eligible", "canonical splice site destroyed; an ASO cannot restore it (Zardetto 2024)"
        if rna in ("exon_skip", "intron_retention"):
            return "not eligible", "the variant weakens the canonical site rather than creating a new one"
        if d < 5:
            tier, why = "not eligible", "under the 5-nt hard minimum distance from the splice site"
        elif d < 15:
            tier, why = "unlikely eligible", "under the recommended 15-nt distance from the splice site"
        elif side == "acceptor" and 18 <= loc["offset"] <= 40:
            tier, why = "unlikely eligible", "in the branch-point window 18-40 nt upstream of the acceptor"
        elif d > 100:
            tier, why = "likely eligible", "deep intronic (>100 nt), the ideal case"
        else:
            tier, why = "likely eligible", "at least 15 nt from the splice site"
        if rna in ("cryptic_site", "pseudoexon"):
            return tier, why
        return "unable to assess", f"needs RNA evidence of a cryptic site; would be '{tier}' ({why})"
    if loc["region"] == "exon" and consequence in SYNONYMOUS | MISSENSE and rna in ("cryptic_site", "pseudoexon"):
        exon = next(e for e in EXONS if e[0] == loc["exon"])
        d = min(loc["c"] - exon[1] + 1, exon[2] - loc["c"] + 1)
        if d < 15:
            return "unlikely eligible", "exonic cryptic site under 15 nt from the exon end"
        if consequence in MISSENSE:
            return "unlikely eligible", "the amino-acid change itself may be pathogenic"
        return "likely eligible", "synonymous exonic cryptic site at least 15 nt from the exon end"
    return "not applicable", "no mis-splicing mechanism to correct"


def exon_skipping(loc, consequence):
    if loc["region"] != "exon" or consequence not in TRUNCATING | MISSENSE | INFRAME:
        return "not applicable", "exon skipping targets coding variants"
    exon = loc["exon"]
    if exon in (EXONS[0][0], EXONS[-1][0]):
        return "not eligible", "first or last coding exon"
    if exon not in IN_FRAME_SKIP:
        return "not eligible", f"exon {exon} is out of frame; skipping it would shift the reading frame"
    tier, why = IN_FRAME_SKIP[exon]
    if consequence in MISSENSE | INFRAME and tier == "likely eligible":
        return "unlikely eligible", "a missense that inactivates BUBR1 implies the exon matters (Zardetto 2024); " + why
    return tier, why


def upregulation(consequence, mech):
    """Raising productive BUB1B mRNA (e.g. by blocking non-productive splicing). Untiered in N1C v1.0."""
    if mech.startswith("abundance") and consequence in MISSENSE:
        return "candidate", "more of a protein that works when abundant (I909T, L1012P overexpression; Suijkerbuijk 2010)"
    if mech.startswith("abundance") and consequence not in TRUNCATING:
        return "candidate", "if this allele still makes some normal transcript"
    if consequence in TRUNCATING:
        return "not for this allele", "its transcripts carry the premature stop; benefit would come from the other allele"
    if mech.startswith("functional"):
        return "unlikely to help", "more of a protein with a specific defect"
    return "unclassified", "decide after assay 1"


def classify(hgvs_c, consequence, hgvs_p=None, ptc_c=None, rna="none"):
    if rna not in RNA_EVIDENCE:
        raise ValueError(f"rna must be one of {sorted(RNA_EVIDENCE)}")
    loc = parse_c(hgvs_c)
    change = protein_change(hgvs_p)
    note = None
    if consequence in TRUNCATING and ptc_c is None:
        if change and change[0] != "fs" and change[2] == "*":
            ptc_c = 3 * change[1] - 2
        elif change and change[0] == "fs" and change[2]:
            ptc_c = 3 * change[2] - 2
        elif loc["region"] == "exon":
            ptc_c, note = loc["c"], "stop position unknown; used the variant position (approximate)"
    nmd_call = nmd(ptc_c) if consequence in TRUNCATING else "not applicable"
    mech, mech_why = mechanism(consequence, loc, change, nmd_call)
    up, up_why = upregulation(consequence, mech)
    result = {
        "hgvs_c": hgvs_c, "consequence": consequence, "location": loc, "ptc_c": ptc_c, "nmd": nmd_call,
        "mechanism_class": mech, "mechanism_basis": mech_why,
        "aso": {
            "splice_correction": splice_correction(loc, consequence, rna),
            "exon_skipping": exon_skipping(loc, consequence),
            "knockdown": ("not eligible", "MVA1 is caused by too little BUBR1; lowering it further is the harm"),
            "upregulation_untiered": (up, up_why),
        },
    }
    if note:
        result["note"] = note
    return result


def genotype(a, b):
    """Pre-registered consequences of an allele pair for the Track 2 ranking (worksheet step 4)."""
    out = []
    pair = (a, b)
    if any(x["aso"]["splice_correction"][0] == "likely eligible" for x in pair):
        out.append("H4 -> rank 1: splice correction of one allele (one working allele is expected to suffice in a recessive disease)")
    if any(x["aso"]["exon_skipping"][0] == "likely eligible" for x in pair):
        out.append("H4 -> rank 1, gated: the exon-skipped BUBR1 must be stable in assay 1 before further work")
    elif any(x["location"].get("exon") == 17 and x["consequence"] in TRUNCATING for x in pair):
        out.append("H4 stays rank 4; first test: is exon 17-skipped BUBR1 stable in patient cells (plasmid or ASO)?")
    if not out and any(x["aso"]["splice_correction"][0] == "unable to assess"
                       and "would be 'likely eligible'" in x["aso"]["splice_correction"][1] for x in pair):
        out.append("H4 stays rank 4: RNA from patient cells becomes the first experiment")
    classes = [x["mechanism_class"] for x in pair]
    if any(c.startswith("functional") for c in classes):
        out.append("goal A weakened for the stable allele: rank H3 first, then H1, H2, H5; H4 only if a splice or skip lane applies")
    if any(c.startswith("abundance") for c, x in zip(classes, pair) if x["consequence"] in MISSENSE):
        out.append("H1 and H2 keep their rationale (overexpressed kinase-domain mutants regain function); HSP90 inhibitors are doubly excluded")
    if all(x["consequence"] in TRUNCATING and x["nmd"] == "predicted" for x in pair):
        out.append("both alleles look null: re-check NMD escape and leaky splicing before interpreting")
    if not out:
        out.append("no pre-registered rule fires; H4 stays conditional (rank 4)")
    return out


def exon_table():
    lines = ["| Exon | Coding nt | Frame | Residues | Skip verdict |", "|---|---|---|---|---|"]
    for exon, start, end, _ in EXONS:
        n = end - start + 1
        frame = "in-frame" if n % 3 == 0 else "out of frame"
        res = f"{(start + 2) // 3}-{min((end + 2) // 3, CDS_LEN // 3 - 1)}"
        if exon in (1, 23):
            frame = "n/a"
            verdict = "not eligible (first or last coding exon)"
        elif exon in IN_FRAME_SKIP:
            verdict = IN_FRAME_SKIP[exon][0]
        else:
            verdict = "not eligible (out of frame)"
        lines.append(f"| {exon} | {n} | {frame} | {res} | {verdict} |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("hgvs_c", nargs="?")
    ap.add_argument("consequence", nargs="?")
    ap.add_argument("--hgvs-p")
    ap.add_argument("--ptc-c", type=int)
    ap.add_argument("--rna", default="none", choices=sorted(RNA_EVIDENCE))
    ap.add_argument("--tsv", help="D04 allele table: hgvs_c, consequence, hgvs_p, ptc_c, rna (tab-separated)")
    ap.add_argument("--out", help="JSON output for --tsv; must be outside the repository")
    ap.add_argument("--exon-table", action="store_true")
    args = ap.parse_args()
    if args.exon_table:
        print(exon_table())
        return
    if args.tsv:
        if not args.out:
            ap.error("--tsv needs --out (patient-derived output stays under D:\\mva-hackathon)")
        out = assert_outside_repo(args.out)
        with open(args.tsv, newline="", encoding="utf-8") as fh:
            rows = [classify(r["hgvs_c"], r["consequence"], r.get("hgvs_p") or None,
                             int(r["ptc_c"]) if r.get("ptc_c") else None, r.get("rna") or "none")
                    for r in csv.DictReader(fh, delimiter="\t")]
        report = {"alleles": rows, "genotype": genotype(*rows[:2]) if len(rows) >= 2 else []}
        out.write_text(json.dumps(report, indent=1), encoding="utf-8")
        print(f"wrote {out}")
        return
    if not (args.hgvs_c and args.consequence):
        ap.error("give an HGVS c. name and a consequence, or --tsv, or --exon-table")
    print(json.dumps(classify(args.hgvs_c, args.consequence, args.hgvs_p, args.ptc_c, args.rna), indent=1))


if __name__ == "__main__":
    main()
