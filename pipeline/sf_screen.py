#!/usr/bin/env python3
"""Secondary-findings screen against the ACMG SF v3.3 gene list (84 genes).

Source: Lee K et al., Genet Med 2025;27:101454 (PMID 40568962), Table 1, read from the PMC
full text on 2026-09-27. Inheritance and reporting rules follow the table:
  AD and SD genes   any P/LP variant is reportable
  AR genes          two P/LP variants (or a homozygote); a single one is a carrier finding
  HFE               only p.Cys282Tyr homozygotes
  XL (GLA, OTC)     hemizygous, heterozygous or homozygous P/LP
  ABCD1             hemizygous, homozygous, or two heterozygous P/LP
  TTN and other DCM genes carry extra conditions in the policy text; they are flagged for review.

X-linked calls depend on the child's sex. It is taken ONLY from the phenotype document
(--sex); it is never inferred. With --sex unknown, X-linked findings are reported as
'depends on stated sex'.

Input is the annotated table used by comphet.py plus a 'classification' column
(Pathogenic, Likely pathogenic, ...) and an optional 'hgvsp' column. It writes a detailed
TSV outside the repository (patient data) and prints only counts, which is what may go
into public reports.

Run: python -m pipeline.sf_screen ANNOTATED.tsv OUT.tsv --sex unknown
"""
import argparse
import csv
from collections import Counter

from pipeline.common import assert_outside_repo

AD = "AD"
SF_V33 = {
    # cancer predisposition
    **{g: AD for g in ("APC", "RET", "BRCA1", "BRCA2", "PALB2", "SDHD", "SDHAF2", "SDHC", "SDHB", "MAX",
                       "TMEM127", "BMPR1A", "SMAD4", "TP53", "MLH1", "MSH2", "MSH6", "PMS2", "MEN1", "NF2",
                       "STK11", "PTEN", "RB1", "TSC1", "TSC2", "VHL", "WT1")},
    "MUTYH": "AR",
    # cardiovascular
    **{g: AD for g in ("FBN1", "TGFBR1", "TGFBR2", "SMAD3", "ACTA2", "MYH11", "PKP2", "DSP", "DSC2", "TMEM43",
                       "DSG2", "RYR2", "TNNT2", "LMNA", "FLNC", "TTN", "BAG3", "DES", "RBM20", "TNNC1", "PLN",
                       "COL3A1", "APOB", "PCSK9", "MYH7", "MYBPC3", "TNNI3", "TPM1", "MYL3", "ACTC1", "PRKAG2",
                       "MYL2", "KCNQ1", "KCNH2", "SCN5A", "CALM1", "CALM2", "CALM3")},
    "CASQ2": "AR", "TRDN": "AR", "LDLR": "SD",
    # inborn errors of metabolism
    "BTD": "AR", "CYP27A1": "AR", "GLA": "XL", "HFE": "AR_HFE", "OTC": "XL", "GAA": "AR", "ABCD1": "XL_ABCD1",
    # other phenotypes
    **{g: AD for g in ("ACVRL1", "ENG", "TTR", "RYR1", "CACNA1S", "HNF1A")},
    "RPE65": "AR", "ATP7B": "AR",
}
POLICY_TEXT_REVIEW = {"TNNT2", "LMNA", "FLNC", "TTN", "BAG3", "DES", "RBM20", "TNNC1", "PLN"}
PLP = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
HOM = {"1/1", "1|1"}


def screen(rows, sex="unknown"):
    if sex not in ("female", "male", "unknown"):
        raise ValueError("sex must be female, male or unknown, as stated in the phenotype document")
    by_gene = {}
    for v in rows:
        if v.get("gene") in SF_V33 and v.get("classification") in PLP:
            by_gene.setdefault(v["gene"], []).append(v)
    findings = []
    for gene, vs in sorted(by_gene.items()):
        rule = SF_V33[gene]
        homs = [v for v in vs if v.get("gt") in HOM]
        if rule in ("AD", "SD"):
            status = "reportable"
        elif rule == "AR":
            status = "reportable (biallelic; check phase)" if homs or len(vs) >= 2 else "carrier (not reportable)"
        elif rule == "AR_HFE":
            c282y = [v for v in homs if "Cys282Tyr" in (v.get("hgvsp") or "")]
            status = "reportable" if c282y else "not reportable (HFE: C282Y homozygotes only)"
        elif sex == "unknown":
            status = "depends on stated sex (X-linked)"
        elif rule == "XL":
            status = "reportable"
        else:  # ABCD1: hemizygous, homozygous or two heterozygous
            status = "reportable" if sex == "male" or homs or len(vs) >= 2 else "carrier (not reportable)"
        if gene in POLICY_TEXT_REVIEW and status == "reportable":
            status = "reportable if policy-text conditions are met (review)"
        findings.append({"gene": gene, "inheritance": rule.split("_")[0], "n_plp": len(vs), "status": status})
    return findings


def public_counts(findings):
    """Only these counts may appear in a public report."""
    c = Counter("reportable" if f["status"].startswith("reportable") else "other" for f in findings)
    return {"sf_genes_screened": len(SF_V33), "genes_with_plp": len(findings),
            "reportable_or_review": c["reportable"], "carrier_or_not_reportable": c["other"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("table")
    ap.add_argument("out")
    ap.add_argument("--sex", default="unknown", help="from the phenotype document only; never inferred")
    a = ap.parse_args()
    assert_outside_repo(a.out)
    with open(a.table, encoding="utf-8", newline="") as fh:
        findings = screen(list(csv.DictReader(fh, delimiter="\t")), a.sex)
    with open(a.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "inheritance", "n_plp", "status"], delimiter="\t")
        w.writeheader()
        w.writerows(findings)
    print(public_counts(findings))


if __name__ == "__main__":
    main()
