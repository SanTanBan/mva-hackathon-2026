#!/usr/bin/env python3
"""Fill the Track 1 methods report from the mva-d02 outputs, for one submitted model.

Reads the public template reports/track1/methods_report.md and the outputs that
tools/kaggle_d02.py collected (D:\\mva-hackathon\\work\\d02\\<date>: the ranked CSV, the rows of
the ten MVA genes plus ranked rows, baf_by_chromosome.tsv and summary.json). Writes the source
with private markers to D:\\mva-hackathon\\work\\reports\\track1\\methods_report.src.md. Render it
with tools/render_reports.py, which withholds every private block from the public version.

Patient-level detail (coordinates, HGVS names, genotypes) goes only inside private markers.
Public text names genes, consequence classes and counts, never coordinates.

Run: python tools/draft_track1_report.py --d02 D:\\mva-hackathon\\work\\d02\\<date> --model 1
"""
import argparse
import csv
import json
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.common import assert_outside_repo  # noqa: E402

TEMPLATE = ROOT / "reports" / "track1" / "methods_report.md"
OUT = Path(r"D:\mva-hackathon\work\reports\track1\methods_report.src.md")
REPO_URL = "https://github.com/SanTanBan/mva-hackathon-2026"
MVA_GENES = ["BUB1B", "CEP57", "TRIP13", "CENATAC", "SLF2", "SMC5", "MAD1L1", "MAD2L1BP", "CEP192", "BUB1"]


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def section(text, heading, body):
    """Replace the body under '### heading' (up to the next ### or ## or ---) with body."""
    pat = re.compile(rf"(^### {re.escape(heading)}\n)(.*?)(?=^### |^## |^---)", re.S | re.M)
    if not pat.search(text):
        raise SystemExit(f"section not found: {heading}")
    return pat.sub(lambda m: m.group(1) + body.rstrip() + "\n\n", text, count=1)


def append_to(text, heading, extra):
    """Add extra at the end of the body under '### heading', keeping what is there."""
    pat = re.compile(rf"(^### {re.escape(heading)}\n.*?)(?=^### |^## |^---)", re.S | re.M)
    if not pat.search(text):
        raise SystemExit(f"section not found: {heading}")
    return pat.sub(lambda m: m.group(1).rstrip() + "\n\n" + extra.rstrip() + "\n\n", text, count=1)


# What each model has actually run, so the methods text never claims a step that is still pending.
STATUS_MODEL_1 = {
    "A3. Model or approach in detail":
        "**Status in model 1:** steps 1, 2 and the ranking in step 6 ran on the child's VCF, plus the "
        "B-allele-fraction part of step 4. The read-level and splice checks of step 3, low-pass depth, the "
        "coverage check and secondary-findings screen of step 5, and the ACMG points of step 6 are written "
        "and tested on synthetic data, and follow in later models.",
    "B5. *BUB1B* locus deep-dive":
        "**Status in model 1:** VCF level only. The Pangolin lookup of the child's variants, k-mer baiting "
        "and read-level phase follow in a later model.",
    "B6. Evidence for mosaic aneuploidy":
        "**Status in model 1:** B-allele fractions from the VCF only (results in §B10); low-pass depth follows.",
    "B7. Differential genes and secondary findings":
        "**Status in model 1:** variant counts for the ten MVA genes are in §B10; the coverage check and the "
        "ACMG SF v3.3 screen follow.",
    "B8. ACMG/AMP classification":
        "**Status in model 1:** not yet applied; the ranking uses the pre-registered allele-evidence weights in §B9.",
}


STATUS_MODEL_2 = {
    "A3. Model or approach in detail":
        "**Status in model 2:** the ranking and the submission file are model 1's, unchanged. Model 2 adds the "
        "allele-level work of steps 3, 5 and 6 on the child's data: the splice lookup and phase of the *BUB1B* "
        "alleles, their mechanism class, a structure-based stability prediction, the coverage check of the other "
        "MVA genes, the ACMG SF v3.3 screen and the per-allele ACMG points. Low-pass depth and k-mer baiting were "
        "not run: both need the FASTQ files inside the private compute, and two rare *BUB1B* alleles were already found.",
    "B6. Evidence for mosaic aneuploidy":
        "**Status in model 2:** B-allele fractions from the VCF only (results in §B10). Low-pass depth was not run.",
}
# Evidence codes, fixed with this script before model 2's report was drafted (2026-10-23). ClinGen SVI
# guidance: PM2 at Supporting; PM3 scored by points, 0.5 for a P/LP partner of unknown phase (Supporting).
# PP4 is not applied, because MVA has several causal genes. PP3 and PS3 are not applied: the stability
# and ESM scores are not calibrated for ACMG strength, and no functional data exist for the child's alleles.
PM2_MAX_AF = 1e-4


def acmg_codes(allele, partner_is_plp, pop_af):
    codes = []
    if allele["consequence"] in ("stop_gained", "frameshift", "frameshift_variant") and allele["nmd"] == "predicted":
        codes.append("PVS1")
    if pop_af is not None and pop_af <= PM2_MAX_AF:
        codes.append("PM2_Supporting")
    if partner_is_plp:
        codes.append("PM3_Supporting")
    return codes


def baf_leads(baf):
    """Choose which rows of baf_by_chromosome.tsv to test, and which of those to flag.

    baf: dicts with keys chrom, n_het, median_depth, centre, delta, lr_stat,
    trisomy_fraction_if_gain, monosomy_fraction_if_loss, robust_z (all strings). The table covers
    every contig in the VCF, including about 1,000 decoy and unplaced ones.
    Returns (tested, flagged): two lists of those dicts, flagged being a subset of tested.
    """
    # Autosomes only. Decoy and unplaced contigs are mostly collapsed repeats, so their allele balance
    # is wide for mapping reasons. X and Y are left out because their expected balance depends on sex,
    # which this pipeline never infers. A floor of 1,000 heterozygous sites keeps the fit stable.
    # Flag one-sided (robust z >= 3): a gain or a loss in a fraction of cells widens the spread around
    # 0.5, so a narrower-than-usual spread is not a lead. robust_z comes from pipeline/baf_mosaic.py
    # (median and MAD over the autosomes); the >= 3 cut-off was in this script before any results.
    autosomes = {str(i) for i in range(1, 23)}

    def num(r, key):
        try:
            return float(r.get(key) or "nan")
        except ValueError:
            return float("nan")

    tested = [r for r in baf if r.get("chrom", "").removeprefix("chr") in autosomes and num(r, "n_het") >= 1000]
    flagged = [r for r in tested if num(r, "robust_z") >= 3]
    return tested, flagged


def private(block):
    return "<!-- private:start -->\n" + block.rstrip() + "\n<!-- private:end -->\n"


def stability_context(stab):
    """Descriptive context beside the pre-registered call; it never changes the call."""
    ab = [stab["controls"][k]["ddg_destab"] for k in ("R727C", "R814H", "L844F", "I909T", "L1012P")]
    nd = max(stab["controls"][k]["ddg_destab"] for k in ("R550Q", "Q921H"))
    v = stab["ddg_destab"]
    where = ("within the range of" if min(ab) <= v <= max(ab) else "above" if v > max(ab) else "below")
    return (f"For context only, outside the pre-registered call: its predicted destabilisation is {where} the five tested "
            f"abundance-class substitutions and {'above' if v > nd else 'not above'} both no-defect controls.")


def model2_sections(t, d04, d02, extract):
    """Model 2: allele-level results (D04), other MVA genes and the SF screen (D05), ACMG points."""
    from pipeline.acmg_points import classify  # noqa: E402

    alleles = json.loads((d04 / "d08_alleles.json").read_text(encoding="utf-8"))["alleles"]
    splice = [json.loads(ln) for ln in (d04 / "splice_lookup.jsonl").read_text(encoding="utf-8-sig").splitlines() if ln.strip()]
    stab_path = d04 / "stability.json"
    stab = json.loads(stab_path.read_text(encoding="utf-8")) if stab_path.is_file() else None
    by_c = {s["result"]["hgvs_c"]: s for s in splice}
    bub = [r for r in extract if r["gene"] == "BUB1B"]

    def af(r):
        try:
            return float(r.get("pop_af") or 0)
        except ValueError:
            return 0.0

    rows = []
    for a in alleles:
        chrom, pos, ref, alt = by_c[a["hgvs_c"]]["query"].split(":")
        r = next(x for x in bub if x["pos"] == pos and x["ref"] == ref and x["alt"] == alt)
        rows.append((a, r))
    kind = {a["hgvs_c"]: ("truncating" if a["consequence"] == "stop_gained" else
                          "missense" if a["consequence"].startswith("missense") else a["consequence"]) for a, _ in rows}
    plp = {a["hgvs_c"]: "athogenic" in r.get("clnsig", "") and "onflicting" not in r.get("clnsig", "") for a, r in rows}
    allele_pos = {r["pos"] for _, r in rows}
    phased = any(r.get("pid") for _, r in rows) and len({r.get("pid") for _, r in rows}) == 1
    dist_kb = abs(int(rows[0][1]["pos"]) - int(rows[1][1]["pos"])) / 1000 if len(rows) == 2 else None
    others = [r for r in bub if r["pos"] not in allele_pos and af(r) < 0.01]
    others_coding = [r for r in others if r.get("impact") in ("HIGH", "MODERATE")]
    n_other_looked_up = len(splice) - len(rows)
    max_delta = max(max(abs(s["result"]["gain"]), abs(s["result"]["loss"])) for s in splice)
    upstream = [r for r in bub if r["consequence"].startswith("upstream")]

    # ACMG points per allele; PM3 only when the partner is P/LP (ClinVar) and only at Supporting (phase unknown).
    acmg = []
    for (a, r), (b, rb) in zip(rows, reversed(rows)):
        codes = acmg_codes(a, plp[b["hgvs_c"]], af(r) if r.get("pop_af") else None)
        res = classify(codes)
        acmg.append((a, r, codes, res))
    def class_words(m):  # "abundance (provisional)" -> "in the abundance class (provisional)"
        head, _, rest = m.partition(" (")
        return f"in the {head} class" + (f" ({rest}" if rest else "")
    mech = {a["hgvs_c"]: class_words(a["mechanism_class"]) for a, _ in rows}

    b5 = ["**Results in model 2.**", "",
          f"- **Both alleles found.** The VCF carries two rare *BUB1B* variants on MANE Select: a {kind[rows[0][0]['hgvs_c']]} "
          f"allele and a {kind[rows[1][0]['hgvs_c']]} allele, both heterozygous PASS calls. "
          f"{len(others)} further rare *BUB1B* variant(s) are called, "
          + (f"{len(others_coding)} of them protein-altering." if others_coding else "none protein-altering."),
          f"- **Phase: {'read-backed' if phased else 'unknown'}.** " + (
              f"The two variants are about {dist_kb:.0f} kb apart, beyond read-backed phasing, and no parental data were "
              "provided. Compound heterozygosity is assumed, and PM3 is scored at the reduced strength for unknown phase (§B8)."
              if not phased else "Both variants share a phase set."),
          f"- **Splicing.** In the public Pangolin map of every possible *BUB1B* single-nucleotide change, neither allele "
          f"nor the {n_other_looked_up} other rare variant(s) looked up changes a splice site (largest |Δ| {max_delta:.2f}; "
          "a strong call needs 0.5).",
          f"- **Upstream interval.** {len(upstream)} upstream variants are called, all common in gnomAD v4.1 "
          f"(lowest allele frequency {min(af(r) for r in upstream):.3f})." if upstream else "- **Upstream interval.** No variant is called.",
          f"- **Mechanism class** (rules fixed on 2026-09-29, `reports/track2/d08_worksheet.md`): the truncating allele is "
          f"{mech[rows[0][0]['hgvs_c']]} (nonsense-mediated decay predicted); the missense allele is {mech[rows[1][0]['hgvs_c']]}.",
          ]
    if stab:
        b5.append(f"- **Stability prediction** (rules fixed before the run, `reports/track2/stability_map_rules.md`): "
                  f"validity gate {'passed' if stab['gate']['passes'] else 'failed'} "
                  f"({stab['gate']['abundance_controls_above_both_no_defect']} of 5 abundance-class controls above both "
                  f"no-defect controls); call for the missense allele: **{stab['call']}**. "
                  + stability_context(stab))
    b5.append("- **Not run:** k-mer baiting for transposon insertions and structural variants, which needs the FASTQ files.")
    priv = ["| Allele | HGVS (NM_001211.6) | Protein | Genomic (GRCh38) | Exon | Reads ref,alt | gnomAD v4.1 AF | ClinVar |",
            "|---|---|---|---|---|---|---|---|"]
    for a, r in rows:
        priv.append(f"| {kind[a['hgvs_c']]} | {a['hgvs_c']} | {r.get('aa_change', '')} | chr{r['chrom'].removeprefix('chr')}:{r['pos']} "
                    f"{r['ref']}>{r['alt']} | {a['location'].get('exon', '')} | {r.get('ad', '')} | {r.get('pop_af', '')} | "
                    f"{(r.get('clnsig') or 'not in ClinVar').replace('_', ' ')} |")
    priv.append("")
    priv.append("Splice lookups: " + "; ".join(f"{s['result']['hgvs_c']} gain {s['result']['gain']:.3f}, loss {s['result']['loss']:.3f}"
                                               for s in splice))
    if stab:
        priv.append(f"Stability map, missense allele ({stab['substitution']}): predicted ΔΔG {stab['ddg_destab']:.2f} "
                    f"(percentile {stab['ddg_pct_ref']:.0f} of {stab['reference_n']} reference substitutions), ESM-2 LLR "
                    f"{stab['esm_llr']:.2f}, pLDDT {stab['plddt']:.0f}, relative accessibility {stab['rsa']:.2f}.")
    t = append_to(t, "B5. *BUB1B* locus deep-dive", "\n".join(b5) + "\n\n" + private("\n".join(priv)))

    # B7: other MVA genes and the SF screen.
    others9 = [g for g in MVA_GENES if g != "BUB1B"]
    import statistics
    lines = ["| Gene | Calls | Median depth of calls | ClinVar P/LP | Rare high or moderate impact |", "|---|---|---|---|---|"]
    weak = []
    for g in others9:
        rs = [r for r in extract if r["gene"] == g]
        dps = [int(r["dp"]) for r in rs if str(r.get("dp", "")).isdigit()]
        med = statistics.median(dps) if dps else None
        if med is not None and med < 15:
            weak.append(g)
        n_plp = sum("athogenic" in r.get("clnsig", "") and "onflicting" not in r.get("clnsig", "") for r in rs)
        n_hm = sum(af(r) <= 0.01 and r.get("impact") in ("HIGH", "MODERATE") for r in rs)
        lines.append(f"| *{g}* | {len(rs)} | {med if med is not None else 'n/a'} | {n_plp} | {n_hm} |")
    b7 = ["**Results in model 2.**", "",
          "- **The other nine MVA genes are negative.** None carries a ClinVar pathogenic or likely pathogenic call, and no "
          "gene has two rare protein-altering alleles or a rare homozygous one.", "", "\n".join(lines), ""]
    if weak:
        b7.append(f"- **Coverage caveat:** {', '.join('*' + g + '*' for g in weak)} has few calls at low depth. A VCF does not "
                  "record depth where nothing is called, so a cause there is not excluded by this data.")
    sf_files = list(d02.rglob("sf_findings.tsv"))
    if sf_files:
        findings = read_tsv(sf_files[0])
        rep = [f for f in findings if f["status"].startswith("reportable")]
        sf_rows = read_tsv(sf_files[0].with_name("sf_plp_rows.tsv")) if sf_files[0].with_name("sf_plp_rows.tsv").is_file() else []
        conflicting = sum("onflicting" in r.get("classification", "") for r in sf_rows)
        b7.append(f"- **ACMG SF v3.3 screen** (84 genes; ClinVar pathogenic or likely pathogenic calls, `pipeline/sf_screen.py`): "
                  f"{len(findings)} gene(s) with such a call; {len(rep)} reportable or for review; "
                  f"{len(findings) - len(rep)} carrier-only or not reportable. "
                  + (f"A further {conflicting} call(s) in these genes carry only conflicting ClinVar classifications and are "
                     "not reportable. " if conflicting else "")
                  + "X-linked findings depend on the sex stated in the phenotype document, which this pipeline never infers.")
        if findings:
            b7.append(private("**SF findings:**\n\n| Gene | Inheritance | P/LP calls | Status |\n|---|---|---|---|\n"
                              + "\n".join(f"| {f['gene']} | {f['inheritance']} | {f['n_plp']} | {f['status']} |" for f in findings)))
    else:
        b7.append("- **ACMG SF v3.3 screen:** not run in this model.")
    t = section(t, "B7. Differential genes and secondary findings",
                "- Coverage and negative results for the other MVA genes (locus dossier §7).\n"
                "- The ACMG SF v3.3 screen, reported as counts only in the public version.\n\n" + "\n".join(b7))

    # B8: ACMG points per allele.
    b8 = ["Points-based classification (Tavtigian 2020, PMID 32720330; `pipeline/acmg_points.py`). Codes were fixed in "
          "`tools/draft_track1_report.py` before this report was drafted. PM2 is applied at Supporting (ClinGen SVI). "
          "PM3 counts only a pathogenic partner, at Supporting because phase is unknown. PP4 is not applied because MVA has "
          "several causal genes. PP3 and PS3 are not applied, because the computational scores are not calibrated for "
          "ACMG strength and no functional data exist for these alleles.", "",
          "| Allele | Evidence | Points | Class |", "|---|---|---|---|"]
    for a, r, codes, res in acmg:
        b8.append(f"| {kind[a['hgvs_c']]} | {', '.join(codes) or 'none'} | {res['points']} | {res['class']} |")
    miss = next((x for x in acmg if kind[x[0]["hgvs_c"]] == "missense"), None)
    if miss:
        pts = miss[3]["points"]
        b8 += ["", f"**What would settle the missense allele.** It stands at {pts} points, a variant of uncertain significance. "
               f"A BUBR1 immunoblot in patient-derived cells showing low protein (PS3, Strong: +4) would bring it to {pts + 4} "
               f"points, {'likely pathogenic' if pts + 4 >= 6 else 'still uncertain'}. Confirming the two alleles in trans by "
               "parental testing would raise PM3 to Moderate (+1). This is the same first experiment that the Track 2 report "
               "puts first (assay 1)."]
    t = section(t, "B8. ACMG/AMP classification", "\n".join(b8))

    t = section(t, "B11. Interpretation",
                "- **The top-ranked pair is the confirmed answer.** Model 1, which submitted the same rows, scored 100 rank points "
                "and F-max 1.000 on the leaderboard: the clinically confirmed variants are the *BUB1B* pair at rank 1.\n"
                "- **The pair fits the classic MVA1 architecture:** a truncating null allele with a missense allele in the "
                "C-terminal region where every tested substitution made an unstable protein (Suijkerbuijk 2010, PMID 20516114).\n"
                "- **Ranking and classification answer different questions.** Under ACMG/AMP the missense allele is still a "
                "variant of uncertain significance (§B8), yet the gene-first, inheritance-aware ranking put the pair first. "
                "A pipeline that kept only ClinVar or ACMG pathogenic calls would have reported the truncating allele alone, "
                "worth half credit; this is why the rule keeps every rare *BUB1B* variant whatever its predicted impact (§B9).\n"
                "- **Mechanism feeds Track 2.** Both alleles are classed as abundance-class, so raising BUBR1 is the coherent "
                "therapeutic goal, and assay 1 (BUBR1 protein in patient cells) is the experiment that settles both the "
                "classification and the drug hypotheses.")
    t = append_to(t, "B12. Limitations",
                  "- The missense allele's mechanism class rests on its position and a computational stability prediction, not on an "
                  "experiment.\n- Low-pass depth and k-mer baiting were not run in this model.")
    return t


def variant_key(chrom, pos, ref, alt):
    return (chrom, str(pos), ref, alt)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--d02", required=True, help="folder collected by tools/kaggle_d02.py")
    ap.add_argument("--model", type=int, required=True, help="submission number (1-6)")
    ap.add_argument("--out", default=str(OUT), help="source to write (outside the repository)")
    ap.add_argument("--d04", help="model 2+: folder with d08_alleles.json, splice_lookup.jsonl and stability.json")
    args = ap.parse_args()
    d02 = Path(args.d02)
    csv_path = next(d02.rglob("*_track1.csv"))
    summary = json.loads(next(d02.rglob("summary.json")).read_text(encoding="utf-8"))
    extract = read_tsv(next(d02.rglob("mva_genes_and_ranked_rows.tsv")))
    baf_files = list(d02.rglob("baf_by_chromosome.tsv"))
    baf = read_tsv(baf_files[0]) if baf_files else []
    with open(csv_path, encoding="utf-8") as fh:
        ranked = list(csv.DictReader(fh))
    by_key = {}
    for r in extract:
        by_key.setdefault(variant_key(r["chrom"], r["pos"], r["ref"], r["alt"]), []).append(r)
    severity = {"HIGH": 3, "MODERATE": 2, "LOW": 1}

    def describe(chrom, pos, ref, alt, gene=None):
        # A variant can lie in two genes' intervals (BUB1B's 3' end is inside PAK6): prefer the
        # ranked gene, then the most severe impact.
        hits = by_key.get(variant_key(chrom, pos, ref, alt), [])
        hits = [h for h in hits if h["gene"] == gene] or hits
        r = max(hits, key=lambda h: severity.get(h.get("impact", ""), 0), default={})
        hgvs = " ".join(x for x in (r.get("dna_change", ""), r.get("aa_change", "")) if x)
        return r, hgvs

    # Private table of the ranked rows.
    lines = ["| Rank | Gene | Variant 1 | Variant 2 | Consequences | ClinVar | Population AF (source) | Genotype / phase | EPCR |",
             "|---|---|---|---|---|---|---|---|---|"]
    top_public = None
    for i, row in enumerate(ranked, 1):
        note_gene = (row.get("notes", "").split() or [None])[0]  # notes start with the ranked gene
        v1, h1 = describe(row["chrom_1"], row["pos_1"], row["ref_1"], row["alt_1"], note_gene)
        cells2 = ("", None, "")
        if row.get("pos_2"):
            v2, h2 = describe(row["chrom_2"], row["pos_2"], row["ref_2"], row["alt_2"], v1.get("gene"))
            cells2 = (f"{row['chrom_2']}:{row['pos_2']} {row['ref_2']}>{row['alt_2']} {h2}".strip(), v2, h2)
        gene = v1.get("gene", "")
        cons = v1.get("consequence", "") + (f" + {cells2[1].get('consequence', '')}" if cells2[1] else "")
        clin = v1.get("clnsig", "") + (f" + {cells2[1].get('clnsig', '')}" if cells2[1] else "")
        af = f"{v1.get('pop_af', '')} ({v1.get('pop_af_source', '')})" + (
            f" + {cells2[1].get('pop_af', '')} ({cells2[1].get('pop_af_source', '')})" if cells2[1] else "")
        gt = f"{v1.get('gt', '')} {v1.get('pgt', '')}".strip() + (
            f" + {cells2[1].get('gt', '')} {cells2[1].get('pgt', '')}".strip() if cells2[1] else "")
        cells = [str(i), gene, f"{row['chrom_1']}:{row['pos_1']} {row['ref_1']}>{row['alt_1']} {h1}".strip(),
                 cells2[0], cons, clin, af, f"{gt}; {row.get('notes', '')}", row["epcr"]]
        lines.append("| " + " | ".join(c.replace("|", "\\|") for c in cells) + " |")  # phased GTs contain '|'
        if i == 1:  # public text: gene, consequence class, ClinVar class and phase only, never the allele
            phase = re.search(r"phase (unknown|trans|cis)|in (trans|cis)\b", row.get("notes", ""))
            top_public = (gene, v1, cells2[1], phase.group(0) if phase else "")
    ranked_table = "\n".join(lines)

    # Public counts for the ten MVA genes.
    def rare(r):
        try:
            return float(r.get("pop_af") or 0) <= 0.01
        except ValueError:
            return True
    gene_counts = []
    for g in MVA_GENES:
        rows = [r for r in extract if r["gene"] == g]
        hi = [r for r in rows if rare(r) and r.get("impact") in ("HIGH", "MODERATE")]
        gene_counts.append(f"| {g} | {len(rows)} | {sum(1 for r in rows if rare(r))} | {len(hi)} |")
    gene_table = "\n".join(["| Gene | Variant rows | Rare (AF ≤ 1%) | Rare, high or moderate impact |",
                            "|---|---|---|---|"] + gene_counts)

    # BAF summary. Public text gives counts only; which chromosomes deviate stays in the private block.
    baf_note = "The per-chromosome B-allele-fraction table was not produced in this run."
    baf_private = ""
    if baf:
        tested, flagged = baf_leads(baf)
        baf_note = (f"B-allele fractions were tested on {len(tested)} chromosomes; "
                    + (f"{len(flagged)} deviate from the autosomal pattern. A deviation is a lead to check "
                       "against read depth, not a diagnosis: the smallest chromosomes, which are GC-rich or acrocentric, "
                       "often show a wider allele balance in short-read data for mapping reasons." if flagged else "none deviates from the autosomal pattern.")
                    + " Variegated mosaicism can sit below this test's detection limit even when a karyotype shows it.")
        if flagged:
            baf_private = private("**Chromosomes flagged by the B-allele-fraction screen:** "
                                  + ", ".join(f"{r['chrom']} (robust z {r['robust_z']})" for r in flagged))

    t = TEMPLATE.read_text(encoding="utf-8")
    today = date.today().isoformat()
    t = t.replace("# Track 1 methods report (skeleton)", "# Track 1 methods report")
    t = re.sub(r"\*\*Rare Disease, Real Kid: The MVA Hackathon 2026 · Track 1 · [^*]+\*\*",
               f"**Rare Disease, Real Kid: The MVA Hackathon 2026 · Track 1 · model {args.model}, {today}**", t)
    t = re.sub(r"> \*\*How this file works\.\*\*.*?\n\n---", "---", t, flags=re.S)
    t = section(t, "A2. Model number", f"{args.model}.")
    t = section(t, "A5. Automated output or manually curated?",
                "Automated. The ranked rows are the direct output of the pipeline (`kaggle/mva_d02/d02.py`, "
                "rule pre-registered in Part B §B9 before data access); no row or EPCR value was edited by hand.")
    t = section(t, "A6. Downstream manual review", "None for this model.")
    if args.model == 1:
        for heading, status in STATUS_MODEL_1.items():
            t = append_to(t, heading, status)
    else:
        for heading, status in STATUS_MODEL_2.items():
            t = append_to(t, heading, status)
        t = section(t, "A6. Downstream manual review",
                    "The ranked rows and EPCR values are model 1's, unedited. The allele-level results in §B5, §B7 and "
                    "§B8 were produced by scripts whose rules were fixed before they ran "
                    "(`reports/track2/d08_worksheet.md`, `reports/track2/stability_map_rules.md`, the evidence codes in "
                    "`tools/draft_track1_report.py`), and were read by the author.")
    t = section(t, "A12. Run time and cost",
                f"- **Gated data:** one private Kaggle CPU notebook (free tier); annotation, ranking and the "
                f"mosaicism statistics took {summary.get('minutes_total', '?')} minutes in total.\n"
                "- **Public references and tests:** one private Kaggle CPU notebook (`kaggle/mva_refs_build`).\n"
                "- **Public splice map:** one Kaggle GPU session.\n- **Cost:** $0.")
    t = section(t, "A1. Team name", "SantanuBanerjee (Hugging Face username, also used as the display name).")
    gene, a1, a2, phase = top_public or ("", {}, None, "")

    def allele_class(a):
        clin = (a.get("clnsig") or "").replace("_", " ")
        return f"a {a.get('consequence') or 'non-coding'} allele" + (f" (ClinVar: {clin})" if clin else " (not in ClinVar)")
    top_sentence = (f"The top-ranked genotype is in *{gene}*: " + (
        f"a compound heterozygous pair, {phase or 'phase as scored'}, of {allele_class(a1)} and {allele_class(a2)}"
        if a2 else f"{allele_class(a1)}") + ".") if gene else "No candidate passed the filters."
    t = section(t, "A13. Method abstract (up to 500 words)",
                "We ranked recessive genotypes from the proband's VCF with a rule fixed before data access. "
                "Variants were normalised and annotated against MANE Select transcripts only (bcftools csq), "
                "with ClinVar classifications and gnomAD v4.1 joint filtering allele frequencies for the ten MVA "
                "genes and the ACMG SF v3.3 genes. Rare variants (AF ≤ 1%) were paired within genes, with "
                "read-backed phase from the PGT/PID tags. Each genotype was scored as gene prior × allele "
                "evidence × phase factor: *BUB1B* 1.0, *CEP57* and *TRIP13* 0.6, genes on a PanelApp list "
                "pre-registered on 2026-10-04 0.1, others 0.01. Scores were normalised over the ten written rows "
                "plus a 0.25 reserve for a cause outside them. In *BUB1B* every rare variant was kept whatever "
                "its predicted impact, because second alleles in MVA1 are often non-coding, and a single "
                "*BUB1B* allele also gets its own row as the hypothesis that its partner was missed. "
                f"{top_sentence} Per-chromosome B-allele fractions were fitted as a screen for mosaic "
                "aneuploidy. All patient-data compute ran in private Kaggle notebooks; the code, the "
                "pre-registered rules and a synthetic end-to-end test are in the repository. "
                "Limitations: no parental data (phase only where reads span both variants), no BAM or RNA, "
                "and VCF-level calling cannot see most transposon insertions or deep structural variants.")
    t = t.replace(private("*Submitted version only: the ranked candidate pairs with HGVS names and evidence weight.*"),
                  private("**Ranked rows as submitted (model %d)**\n\n%s" % (args.model, ranked_table)))
    b10 = (private("**Ranked rows, with HGVS names, ClinVar, population frequency and genotype:**\n\n" + ranked_table)
           + "\nPublic summary: " + top_sentence + "\n\n" + gene_table + "\n\n" + baf_note + "\n\n"
           + (baf_private + "\n" if baf_private else "")
           + f"Input: {summary.get('input', {}).get('records', '?')} VCF records; "
             f"{summary.get('table_rows', '?')} annotated (variant, gene) rows; frequency sources "
             f"{summary.get('pop_af_source_counts', {})}.")
    t = section(t, "B10. Results", b10)
    t = section(t, "B11. Interpretation",
                "Each allele's mechanism class (functional or abundance; Suijkerbuijk 2010, PMID 20516114) is "
                "assigned with the rules fixed in advance in `reports/track2/d08_worksheet.md`, and feeds the "
                "Track 2 report. This model reports the ranking; the per-allele ACMG points follow in the next model.")
    t = section(t, "B12. Limitations",
                "- No parental data: phase is known only where reads span both variants; otherwise it is scored as unknown.\n"
                "- No BAM, CRAM or RNA: splice predictions are not RNA evidence, and read-level checks need the FASTQ files.\n"
                "- VCF-level calling misses most transposon insertions and many structural variants.\n"
                "- gnomAD frequencies cover only the ten MVA genes and the ACMG SF genes; elsewhere ClinVar's "
                "frequencies are used, and a missing frequency counts as rare.\n"
                "- The B-allele-fraction screen cannot see low-level variegated mosaicism.")
    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    t = section(t, "B13. Reproducibility",
                f"- Repository: {REPO_URL} (commit `{head}`).\n"
                "- Public references and tests: `kaggle/mva_refs_build/refs_build.py` (GRCh38 analysis set from NCBI, "
                "md5-checked; ClinVar 2026-09-28; gnomAD v4.1 slices; bcftools/samtools/htslib 1.24).\n"
                "- Patient data: `kaggle/mva_d02/d02.py`, run by `tools/kaggle_d02.py` after "
                "`tools/fetch_vcf_first.py`, `tools/d01_intake.py` and `tools/kaggle_patient_upload.py`.\n"
                "- Validation: `tools/validate_track1_csv.py`.")
    if args.model >= 2:
        if not args.d04:
            raise SystemExit("model 2 and later need --d04")
        t = model2_sections(t, Path(args.d04), d02, extract)
    out = assert_outside_repo(Path(args.out))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(t, encoding="utf-8")
    left = re.findall(r"^\*[^*\n]+\*$", t, re.M)
    print(f"written {out}; italic placeholders left: {len(left)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
