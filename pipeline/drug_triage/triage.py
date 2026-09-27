#!/usr/bin/env python3
"""Rank drugs for a set of checkpoint genes against the six-point screen, from public data only.

For each query gene: Open Targets tractability and drug/clinical candidates, ChEMBL potency against the
target, then the same for its modifier targets, plus a named drug list and negative controls. Every drug is
screened with flags.evaluate() and ranked by tier. Output: a TSV and a Markdown table.

Run (network, cached under D:\\mva-hackathon\\refs\\drug_triage_cache):
  python -m pipeline.drug_triage.triage --panel pipeline/drug_triage/panel.json --out reports/track2/drug_triage
Offline rerun from cache: add --offline.
"""
import argparse
import csv
import datetime as dt
import json
import sys
from pathlib import Path

from pipeline.drug_triage import flags as F
from pipeline.drug_triage import sources as S

TIER_ORDER = {"screen-clear": 0, "caution": 1, "insufficient data": 2, "excluded": 3}
STAGE_ORDER = {"APPROVAL": 0, "PHASE_4": 0, "PHASE_3": 1, "PHASE_2_3": 2, "PHASE_2": 2, "PHASE_1_2": 3,
               "PHASE_1": 3, "EARLY_PHASE_1": 4, "PRECLINICAL": 5, None: 6}
SHORT = {"flag": "FLAG", "note": "note", "no flag found": "-", "not assessed": "n/a", "yes": "yes",
         "not established": "no", "not stated": "?"}


def sort_key(row):
    return (TIER_ORDER[row["tier"]], row["paediatric_approval"]["status"] != "yes",
            sum(row[k]["status"] == "flag" for k in F.FLAG_NAMES), STAGE_ORDER.get(row["stage"], 6), row["drug"])


def screen_drug(http, name, chembl, query_gene, dosage_deficient, via):
    drug = S.ot_drug(http, name=name, chembl=chembl)
    label = S.fda_label(http, drug["name"] if drug else name) or (S.fda_label(http, name) if drug else None)
    fl = F.evaluate(drug, label, query_gene, dosage_deficient)
    return {"query_gene": query_gene, "drug": (drug or {}).get("name", name.upper()), "chembl": (drug or {}).get("chembl"),
            "parent": (drug or {}).get("parent") or name.upper(),
            "stage": (drug or {}).get("stage"), "via": via, "tier": F.tier(fl),
            "label": (label or {}).get("_set_id"), **fl}


def run(panel, http, log=print):
    gene_rows, drug_rows = [], []
    for gene, cfg in panel["query_genes"].items():
        dd = cfg.get("dosage_deficient", True)
        tgt = S.ot_target(http, gene)
        pot = S.chembl_potency(http, gene)
        gene_rows.append({"gene": gene, "disorder": cfg.get("disorder", ""), "ensembl": tgt["ensembl"],
                          "ot_candidates": len(tgt["candidates"]),
                          "tractability": ", ".join(f"{t['modality']}:{t['label']}" for t in tgt["tractability"]) or "none",
                          "chembl_target": pot["target"], "chembl_records": pot["records"], "chembl_potent": pot["potent"],
                          "best_pchembl": pot["best_pchembl"]})
        todo = [(c["name"], c["chembl"], f"OT candidate for {gene}") for c in tgt["candidates"]]
        for mod in cfg.get("modifier_targets", []):
            todo += [(c["name"], c["chembl"], f"OT candidate for modifier {mod}") for c in S.ot_target(http, mod)["candidates"]]
        todo += [(n, None, "named in Track 2 report") for n in panel.get("named_drugs", [])]
        for neg in panel.get("negative_controls", []):
            todo += [(c["name"], c["chembl"], f"negative control ({neg})") for c in S.ot_target(http, neg)["candidates"]]
        seen, rows = set(), []
        for name, chembl, via in todo:
            key = (chembl or name).upper()
            if key in seen:
                continue
            seen.add(key)
            log(f"  {gene}: {name}")
            rows.append(screen_drug(http, name, chembl, gene, dd, via))
        # Collapse salt forms and repeat mentions onto the parent molecule, keeping every route it came by.
        by_parent = {}
        for r in rows:
            if r["parent"] in by_parent:
                keep = by_parent[r["parent"]]
                keep["via"] = "; ".join(dict.fromkeys(keep["via"].split("; ") + [r["via"]]))
                if keep["label"] is None and r["label"]:
                    by_parent[r["parent"]] = {**r, "via": keep["via"]}
            else:
                by_parent[r["parent"]] = r
        drug_rows += by_parent.values()
    drug_rows.sort(key=lambda r: (r["query_gene"],) + sort_key(r))
    return gene_rows, drug_rows


def write_tsv(path, drug_rows):
    cols = ["query_gene", "drug", "chembl", "stage", "via", "tier"] + [f"{k}{s}" for k in F.FLAG_NAMES for s in ("", "_evidence")] + ["label_set_id"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in drug_rows:
            w.writerow([r["query_gene"], r["drug"], r["chembl"] or "", r["stage"] or "", r["via"], r["tier"]]
                       + [x for k in F.FLAG_NAMES for x in (r[k]["status"], r[k]["evidence"])] + [r["label"] or ""])


def render_md(gene_rows, drug_rows, when):
    out = [f"# Drug triage over public data (R06 demo, {when})", "",
           "Generated by `pipeline/drug_triage/triage.py` from Open Targets, ChEMBL and openFDA drug labels. "
           "Public data only: no patient data was used or sent. **Every row is a triage result for expert review, "
           "not a recommendation.** \"-\" means nothing matched in the source read, which is not evidence of safety; "
           "\"n/a\" means the source was missing (usually no US label).", "",
           "## Targets", "",
           "| Gene | Disorder | Open Targets candidates | Positive tractability flags | ChEMBL records (pChEMBL) | pChEMBL ≥ 6 | Best pChEMBL |",
           "|---|---|---|---|---|---|---|"]
    for g in gene_rows:
        out.append(f"| {g['gene']} | {g['disorder']} | {g['ot_candidates']} | {g['tractability']} | {g['chembl_records']} | "
                   f"{g['chembl_potent']} | {g['best_pchembl'] if g['best_pchembl'] is not None else 'none'} |")
    out += ["", "Flags: **Paed** US paediatric approval · **Ckpt** lowers a checkpoint protein or the dosage-deficient gene · "
            "**Aneu** aneugen class or label genotoxicity · **CYP3A4** perpetrator interaction · **Kidney** · **Marrow**. "
            "Tier: *excluded* if Ckpt or Aneu flag; *caution* if a label safety flag; *screen-clear* otherwise.", ""]
    for gene in dict.fromkeys(r["query_gene"] for r in drug_rows):
        rows = [r for r in drug_rows if r["query_gene"] == gene]
        counts = {t: sum(r["tier"] == t for r in rows) for t in TIER_ORDER}
        out += [f"## {gene}", "", " · ".join(f"{t}: {n}" for t, n in counts.items()), "",
                "| # | Drug | OT max stage | Source | Tier | Paed | Ckpt | Aneu | CYP3A4 | Kidney | Marrow | Main reason |",
                "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for i, r in enumerate(rows, 1):
            reason = next((r[k]["evidence"] for k in ("checkpoint_weakening", "aneugen_class", "myelosuppression",
                                                     "nephrotoxicity", "cyp3a4_interaction") if r[k]["status"] == "flag"), "")
            reason = reason.replace("|", "/")[:140]
            out.append(f"| {i} | {r['drug']} | {r['stage'] or ''} | {r['via']} | {r['tier']} | "
                       + " | ".join(SHORT[r[k]["status"]] for k in F.FLAG_NAMES) + f" | {reason} |")
        out.append("")
    out += ["## Limits", "",
            "- Label flags are sentence-level text matches on US prescribing information. They find stated warnings; they do not weigh dose, exposure or frequency.",
            "- \"OT max stage\" is Open Targets' maximumClinicalStage. APPROVAL there can reflect any market listing (for example quercetin and nicotinamide riboside, which are sold as supplements or OTC products). Verify approval status before quoting it.",
            "- Paediatric \"yes\" means at least one US prescription label for the ingredient states established paediatric safety and effectiveness, for some indication and age range. \"no\" means labels state it is not established; \"?\" means they do not say.",
            "- Drugs without a US label (supplements, non-US approvals, research compounds) are `n/a` for label flags; their mechanism flags still apply.",
            "- *screen-clear* only means none of these six flags fired. Concerns outside them (for example immunosuppression, effects on growth, liver toxicity) are not screened; the Track 2 report weighs those by hand.",
            "- The checkpoint and aneugen rules encode the Track 2 report's do-not-use list (§5). They are deliberately conservative.",
            "- Mitotic-duration effects, the report's veto assay (§7), cannot be read from any database and must be measured in patient-derived cells.", ""]
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--panel", default=str(Path(__file__).with_name("panel.json")))
    ap.add_argument("--out", required=True, help="output prefix; writes <out>.tsv and <out>.md")
    ap.add_argument("--cache", default=str(S.DEFAULT_CACHE))
    ap.add_argument("--offline", action="store_true")
    a = ap.parse_args(argv)
    panel = json.loads(Path(a.panel).read_text(encoding="utf-8"))
    http = S.Http(a.cache, offline=a.offline)
    gene_rows, drug_rows = run(panel, http, log=lambda m: print(m, file=sys.stderr))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(out.with_suffix(".tsv"), drug_rows)
    out.with_suffix(".md").write_text(render_md(gene_rows, drug_rows, dt.date.today().isoformat()), encoding="utf-8")
    print(f"{len(drug_rows)} drug rows for {len(gene_rows)} genes -> {out}.tsv / .md")


if __name__ == "__main__":
    main()
