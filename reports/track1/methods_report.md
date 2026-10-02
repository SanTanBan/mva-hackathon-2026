# Track 1 methods report (skeleton)

**Rare Disease, Real Kid: The MVA Hackathon 2026 · Track 1 · skeleton v1, 2026-09-30, aligned with the organisers' methods description form**

> **How this file works.** This is the public template. Once analysis starts, the working copy with results lives outside the repository, at `D:\mva-hackathon\work\reports\track1\methods_report.src.md`. From that copy, `python tools/render_reports.py <source> --track track1 --name methods_report` writes:
> - this public file, with every block and inline span marked private withheld (the marker syntax is documented in `tools/render_reports.py`);
> - the full version, into `D:\mva-hackathon\submissions\track1\<date>\`.
>
> Both are produced as Markdown and PDF. The renderer refuses to write a public version that contains variant coordinates, VCF lines or protected identifiers.
>
> **Structure.** Part A answers the organisers' methods description form (`methods_description_form.xlsx`, sheet "Track 1 methods") question by question, in its order. The answers are copied into that form when a submission is staged. Part B gives the detailed methods and results. Italic text marks what is still to be filled in once the analysis has run.

---

## Part A. Methods description form (Track 1)

### A1. Team name
SantanuBanerjee (Hugging Face username). *The display name is confirmed at submission.*

### A2. Model number
*1 for the first submission. The form is filled once per submitted model or approach, up to 6.*

### A3. Model or approach in detail
The approach combines a gene-agnostic pipeline with a mechanism-led search of *BUB1B*:

1. **Start from the VCF.** Quality control and normalisation come first. Consequences are pinned to MANE Select transcripts (Part B §B2–B3).
2. **Rank recessive candidates.** Homozygous and compound heterozygous candidates are ranked against the HPO terms. An independent, panel-based route uses Genomics England PanelApp (Part B §B4).
3. **Search *BUB1B* by mechanism, not by ClinVar status.** This matters because the classic second alleles are VUS or conflicting in ClinVar (locus dossier §3.1). The search covers:
   - splice-region and deep-intronic variants, scored against a precomputed, public splice map of every possible single-nucleotide change in the gene;
   - the promoter and the upstream interval;
   - transposon insertions and structural variants, through k-mer baiting of the locus;
   - read-backed phase (Part B §B5).
4. **Look for mosaic aneuploidy.** The evidence comes from per-chromosome B-allele fractions and low-pass depth (Part B §B6).
5. **Show the differential MVA genes negative.** This includes a coverage check (Part B §B7).
6. **Classify and rank.** Each allele is classified with ACMG/AMP points (Part B §B8), and candidate rows are ranked by EPCR (Part B §B9).

### A4. Generative AI used (required)
Anthropic Claude Max (consumer plan), used through Claude Code, with model-training opt-out enabled and no feedback or ratings given in sessions involving the data. These are Sage's conditions for consumer plans (HF thread #2, 29 September 2026). The wording was confirmed on 2026-09-30 and matches `routine/plan.json` → `gates.llm_terms_confirmed.record`. No patient data was used with an LLM before that date.

### A5. Automated output or manually curated?
*Planned: the ranked rows are the automated output of the pipeline. A documented manual review of the top candidates (ACMG evidence, phase and splice evidence) can change an EPCR value, and any such change is listed in A6.*

### A6. Downstream manual review
*Each manual change: the row, the original and final EPCR, and the reason.*

### A7. Public data only?
Yes. Only public reference data, databases and literature were used, together with the hackathon dataset itself. No proprietary data.

### A8. Public data used
- **Reference:** GRCh38 no-alt analysis set (NCBI), matching the VCF's reference. RefSeq/MANE Select transcripts (NCBI) and Ensembl gene models.
- **Variant and gene knowledge:** ClinVar; gnomAD, for frequency and constraint; Genomics England PanelApp; UniProt features; the ACMG SF v3.3 gene list.
- **Splicing:** Pangolin (Zeng & Li 2022, GPL-3.0), run over every possible single-nucleotide change in *BUB1B* on public sequence only (`kaggle/bub1b_splice_map/`).
- **Literature:** every citation was checked against its PubMed record; the list is in the locus dossier.
- *Versions and download dates are listed in Part B §B10.*

### A9. Proprietary data used
None.

### A10. Compound heterozygous pairs, or single variants only?
Pairs. The pipeline outputs compound heterozygous pairs in the `chrom_1…alt_1` and `chrom_2…alt_2` columns. Phase comes from read-backed phasing tags (PGT/PID) where both variants share reads; otherwise it is scored as unknown, with PM3 at reduced strength. Single variants are reported only for homozygous or hemizygous candidates.

### A11. Secondary or incidental findings
- **How they are found:** a screen of the ACMG SF v3.3 genes (Part B §B7).
- **In the submission file:** an actionable result is a `secondary` row with a short note.
- **In public versions:** only counts appear, never variant details.

### A12. Run time and cost
*Measured after the run.* Planned environment: free tiers only.
- **Gated data:** one laptop (Intel i3-10110U, 4 threads, 7.8 GB RAM), working from the VCF; FASTQ steps run overnight at low priority, one pair at a time.
- **Public splice map:** one Kaggle GPU session.
- **Expected cost:** $0.

### A13. Method abstract (up to 500 words)
*Written after the analysis, covering the strengths and limitations listed in Part B §B12.*

---

## Part B. Detailed methods and results

<!-- private:start -->
*Submitted version only: the ranked candidate pairs with HGVS names and evidence weight.*
<!-- private:end -->

### B1. Data, access and governance
- **Dataset:** whole-genome FASTQ (four pairs), one VCF with its index, and a phenotype document with HPO terms. It is shared under WCG IRB protocol #20252010 through gated access on Hugging Face.
- **Reference build:** GRCh38 (no-alt analysis set with decoys); contigs have no `chr` prefix. Variants were called with Sentieon, then filtered with GATK VariantFiltration.
- **Where the data is processed:** locally, on one laptop. The VCF, its index and the phenotype document are downloaded from the gated repository; FASTQ pairs are fetched one at a time and deleted after use. Gated and derived files stay under `D:\mva-hackathon\`, and nothing gated or derived enters git (a pre-commit guard enforces this).
- **Deletion:** all of it is deleted by 2026-11-23. No patient data went to any other service.

### B2. Compute environment and VCF quality control
- **Compute:**
  - Patient data is processed locally on one laptop (4 threads, 7.8 GB RAM), with bcftools, samtools and htslib 1.24 and minimap2 2.31 under WSL Ubuntu. No patient data leaves the machine.
  - Public-data precomputation (the *BUB1B* splice map) runs on Kaggle; no patient data goes there.
  - Full-genome realignment does not fit this machine, so the analysis works from the VCF, plus targeted read-level checks on FASTQ pairs processed one at a time.
- **VCF quality control:** header, sample and build checks; left-alignment and multi-allelic splitting; quality filters; and the `chr` prefix added on output.

### B3. Annotation and inheritance-aware filtering
- **Consequences** come from `bcftools csq` given only the 19,299 MANE Select transcripts of Ensembl release 116 (`pipeline/mane_gff.py`), so readthrough transcripts cannot give a gene a spurious consequence (locus dossier §1). Indels are left-aligned and multi-allelic records split first (`bcftools norm`). Variants with no coding or splice consequence are placed by gene interval as intronic, upstream or downstream (within 5 kb), so the *BUB1B* upstream interval is not lost (`tools/annotate_vcf.sh`, `pipeline/annotate_table.py`).
- **ClinVar** classification and review status are attached by position and alleles (release of 2026-09-28).
- **Rare-variant filter:** population frequency at most 1%. A missing frequency counts as rare, so the filter can only over-include.
- **Recessive models:** homozygous, and compound heterozygous pairs, with read-backed phase from the PGT/PID tags; pairs without shared phase are scored as phase unknown (`pipeline/comphet.py`). In *BUB1B* every rare variant is kept whatever its predicted impact, because MVA second alleles are often non-coding.
- The whole path was tested end to end on a synthetic VCF built from public ClinVar *BUB1B* alleles (`tools/test_annotate_e2e.py`).

### B4. Phenotype-driven prioritisation
- HPO-based ranking.
- An independent clinical route through Genomics England PanelApp panels, where *BUB1B* is a green, biallelic gene on several phenotype panels (locus dossier §6).

### B5. *BUB1B* locus deep-dive
We seek both alleles by mechanism, not by ClinVar status, because the classic hypomorphic alleles are VUS or conflicting in ClinVar (locus dossier §3.1). The search covers:

- **The promoter and the upstream interval.**
- **Splice-region and deep-intronic variants.**
  - Each is looked up in the public Pangolin map of every possible single-nucleotide change in the gene.
  - Indels are scored individually in the same environment.
- **Transposon insertions and structural variants,** found by k-mer baiting of the locus.
- **Read-backed phase.**
- **Mechanism class and antisense lanes,** applied with the rules fixed in advance in `reports/track2/d08_worksheet.md`.

### B6. Evidence for mosaic aneuploidy
- Per-chromosome B-allele fraction from VCF allele depths.
- Low-pass depth from alignments of reads taken from the start of the FASTQs.
- Stated detection limits, and the caveat that variegated mosaicism is spread across different chromosomes in different cells.

### B7. Differential genes and secondary findings
- Coverage and negative results for the other MVA genes (locus dossier §7).
- The ACMG SF v3.3 screen, reported as counts only in the public version.

### B8. ACMG/AMP classification
Points-based classification shown per allele, with the evidence codes and their strengths (`pipeline/acmg_points.py`).

### B9. Submission file
Candidates and their scores (`epcr`) come from `pipeline/track1_rank.py`, whose rule was fixed on 2026-10-02, before any access to the data:

- **Score = gene prior × allele evidence × phase factor.** Gene prior: *BUB1B* 1.0, *CEP57* and *TRIP13* 0.6, genes on the PanelApp panels for the child's HPO terms 0.1, any other gene 0.01. Allele evidence: ClinVar pathogenic or likely pathogenic 1.0; otherwise by predicted impact (high 0.8, moderate 0.3, low 0.15, modifier 0.05); ClinVar benign 0.01. Phase: trans 1.0, unknown 0.5, cis 0.05.
- **A single *BUB1B* allele whose partner was not found** also gets a row, at 0.3 of its allele score, because MVA second alleles are often non-coding or structural and can be missed by short-read calling.
- **EPCR** is the score divided by the sum of all scores plus 0.25. The 0.25 is held back for a cause the pipeline cannot see, so no row claims certainty.

The file is validated with `tools/validate_track1_csv.py`, which checks the organisers' current field specification: at most 10 rows, `epcr` in (0, 1], `finding_type` either `primary` or `secondary`, and the `chr` prefix.

### B10. Results

<!-- private:start -->
*Submitted version only: variant tables, both alleles with HGVS names, phase evidence, ACMG points, mosaicism figures, differential-gene coverage, and secondary-finding details.*
<!-- private:end -->

Public summary: *what was found, described without coordinates or identifiers.*

### B11. Interpretation
*How the genotype fits MVA1 mechanisms. Each allele is classified as functional class or abundance class (Suijkerbuijk 2010, PMID 20516114), using the rules fixed in advance in the D08 worksheet. Link to the Track 2 report.*

### B12. Limitations
*No BAM or RNA; no parental data; low-pass depth only; what the VCF cannot show (for example transposon insertions); the limits of the mosaicism assay. Splice predictions are not RNA evidence.*

### B13. Reproducibility
*Repository URL, commit hash, environment file, and the exact commands, in order. Public data sources and their versions.*

### B14. Privacy and data deletion
No names, dates, locations more specific than "UK (NHS-diagnosed)", or clinical timeline appear. The public version withholds all patient-level variants. The data must be deleted from every environment by 2026-11-23: the laptop, including the Hugging Face download cache. Deletion is then confirmed by email.

---

## Acknowledgement

*"This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible."*

The organisers require this acknowledgement in any public communication arising from the Hackathon.

References: see `reports/track1/locus_dossier.md` for the full, PubMed-checked reference list.
