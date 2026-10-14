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


def private(block):
    return "<!-- private:start -->\n" + block.rstrip() + "\n<!-- private:end -->\n"


def variant_key(chrom, pos, ref, alt):
    return (chrom, str(pos), ref, alt)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--d02", required=True, help="folder collected by tools/kaggle_d02.py")
    ap.add_argument("--model", type=int, required=True, help="submission number (1-6)")
    ap.add_argument("--out", default=str(OUT), help="source to write (outside the repository)")
    args = ap.parse_args()
    d02 = Path(args.d02)
    csv_path = next(d02.rglob("*_track1.csv"))
    summary = json.loads(next(d02.rglob("summary.json")).read_text(encoding="utf-8"))
    extract = read_tsv(next(d02.rglob("mva_genes_and_ranked_rows.tsv")))
    baf_files = list(d02.rglob("baf_by_chromosome.tsv"))
    baf = read_tsv(baf_files[0]) if baf_files else []
    with open(csv_path, encoding="utf-8") as fh:
        ranked = list(csv.DictReader(fh))
    by_key = {variant_key(r["chrom"], r["pos"], r["ref"], r["alt"]): r for r in extract}

    def describe(chrom, pos, ref, alt):
        r = by_key.get(variant_key(chrom, pos, ref, alt), {})
        hgvs = " ".join(x for x in (r.get("dna_change", ""), r.get("aa_change", "")) if x)
        return r, hgvs

    # Private table of the ranked rows.
    lines = ["| Rank | Gene | Variant 1 | Variant 2 | Consequences | ClinVar | Population AF (source) | Genotype / phase | EPCR |",
             "|---|---|---|---|---|---|---|---|---|"]
    top_public = None
    for i, row in enumerate(ranked, 1):
        v1, h1 = describe(row["chrom_1"], row["pos_1"], row["ref_1"], row["alt_1"])
        cells2 = ("", None, "")
        if row.get("pos_2"):
            v2, h2 = describe(row["chrom_2"], row["pos_2"], row["ref_2"], row["alt_2"])
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
        if i == 1:
            top_public = (gene, v1.get("consequence", ""), cells2[1].get("consequence", "") if cells2[1] else None,
                          row.get("notes", ""))
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

    # Public BAF summary: chromosomes whose robust z-score is notable.
    baf_note = "The per-chromosome B-allele-fraction table was not produced in this run."
    if baf:
        zcol = next((c for c in baf[0] if "z" in c.lower()), None)
        flagged = []
        if zcol:
            for r in baf:
                try:
                    if abs(float(r[zcol])) >= 3:
                        flagged.append(r.get("chrom", "?"))
                except ValueError:
                    pass
        baf_note = (f"B-allele fractions were fitted for {len(baf)} chromosomes. "
                    + (f"{len(flagged)} chromosome(s) deviate from the autosomal pattern (robust |z| ≥ 3): "
                       f"{', '.join(flagged)}. A deviation is a lead to check against read depth, not a diagnosis."
                       if flagged else "None deviates from the autosomal pattern (robust |z| ≥ 3)")
                    + " Variegated mosaicism can sit below this test's detection limit even when a karyotype shows it.")

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
    t = section(t, "A12. Run time and cost",
                f"- **Gated data:** one private Kaggle CPU notebook (free tier); annotation, ranking and the "
                f"mosaicism statistics took {summary.get('minutes_total', '?')} minutes in total.\n"
                "- **Public references and tests:** one private Kaggle CPU notebook (`kaggle/mva_refs_build`).\n"
                "- **Public splice map:** one Kaggle GPU session.\n- **Cost:** $0.")
    gene, c1, c2, note = top_public or ("", "", None, "")
    top_sentence = (f"The top-ranked genotype is in *{gene}*: " + (f"a compound heterozygous pair ({c1} and {c2}; {note})"
                    if c2 else f"a single allele ({c1}; {note})") + ".") if gene else "No candidate passed the filters."
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
    out = assert_outside_repo(Path(args.out))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(t, encoding="utf-8")
    left = re.findall(r"^\*[^*\n]+\*$", t, re.M)
    print(f"written {out}; italic placeholders left: {len(left)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
