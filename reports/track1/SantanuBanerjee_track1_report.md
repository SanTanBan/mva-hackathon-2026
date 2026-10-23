# Track 1 methods report

**Rare Disease, Real Kid: The MVA Hackathon 2026 · Track 1 · model 2, 2026-10-23**

---

## Part A. Methods description form (Track 1)

### A1. Team name
SantanuBanerjee (Hugging Face username, also used as the display name).

### A2. Model number
2.

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

**Status in model 2:** the ranking and the submission file are model 1's, unchanged. Model 2 adds the allele-level work of steps 3, 5 and 6 on the child's data: the splice lookup and phase of the *BUB1B* alleles, their mechanism class, a structure-based stability prediction, the coverage check of the other MVA genes, the ACMG SF v3.3 screen and the per-allele ACMG points. Low-pass depth and k-mer baiting were not run: both need the FASTQ files inside the private compute, and two rare *BUB1B* alleles were already found.

### A4. Generative AI used (required)
Anthropic Claude Max (consumer plan), used through Claude Code, with model-training opt-out enabled and no feedback or ratings given in sessions involving the data. These are Sage's conditions for consumer plans (HF thread #2, 29 September 2026). The wording was confirmed on 2026-09-30 and matches `routine/plan.json` → `gates.llm_terms_confirmed.record`. No patient data was used with an LLM before that date.

### A5. Automated output or manually curated?
Automated. The ranked rows are the direct output of the pipeline (`kaggle/mva_d02/d02.py`, rule pre-registered in Part B §B9 before data access); no row or EPCR value was edited by hand.

### A6. Downstream manual review
The ranked rows and EPCR values are model 1's, unedited. The allele-level results in §B5, §B7 and §B8 were produced by scripts whose rules were fixed before they ran (`reports/track2/d08_worksheet.md`, `reports/track2/stability_map_rules.md`, the evidence codes in `tools/draft_track1_report.py`), and were read by the author.

### A7. Public data only?
Yes. Only public reference data, databases and literature were used, together with the hackathon dataset itself. No proprietary data.

### A8. Public data used
- **Reference:** GRCh38 no-alt analysis set (NCBI), matching the VCF's reference. RefSeq/MANE Select transcripts (NCBI) and Ensembl gene models.
- **Variant and gene knowledge:** ClinVar; gnomAD, for frequency and constraint; Genomics England PanelApp; UniProt features; the ACMG SF v3.3 gene list.
- **Splicing:** Pangolin (Zeng & Li 2022, GPL-3.0), run over every possible single-nucleotide change in *BUB1B* on public sequence only (`kaggle/bub1b_splice_map/`).
- **Literature:** every citation was checked against its PubMed record; the list is in the locus dossier.
- Versions and download dates are listed in Part B §B13.

### A9. Proprietary data used
None.

### A10. Compound heterozygous pairs, or single variants only?
Pairs. The pipeline outputs compound heterozygous pairs in the `chrom_1…alt_1` and `chrom_2…alt_2` columns. Phase comes from read-backed phasing tags (PGT/PID) where both variants share reads; otherwise it is scored as unknown, with PM3 at reduced strength. Single variants are reported only for homozygous or hemizygous candidates.

### A11. Secondary or incidental findings
- **How they are found:** a screen of the ACMG SF v3.3 genes (Part B §B7).
- **In the submission file:** an actionable result is a `secondary` row with a short note.
- **In public versions:** only counts appear, never variant details.

### A12. Run time and cost
- **Gated data:** one private Kaggle CPU notebook (free tier); annotation, ranking and the mosaicism statistics took 8.6 minutes in total.
- **Public references and tests:** one private Kaggle CPU notebook (`kaggle/mva_refs_build`).
- **Public splice map:** one Kaggle GPU session.
- **Cost:** $0.

### A13. Method abstract (up to 500 words)
We ranked recessive genotypes from the proband's VCF with a rule fixed before data access. Variants were normalised and annotated against MANE Select transcripts only (bcftools csq), with ClinVar classifications and gnomAD v4.1 joint filtering allele frequencies for the ten MVA genes and the ACMG SF v3.3 genes. Rare variants (AF ≤ 1%) were paired within genes, with read-backed phase from the PGT/PID tags. Each genotype was scored as gene prior × allele evidence × phase factor: *BUB1B* 1.0, *CEP57* and *TRIP13* 0.6, genes on a PanelApp list pre-registered on 2026-10-04 0.1, others 0.01. Scores were normalised over the ten written rows plus a 0.25 reserve for a cause outside them. In *BUB1B* every rare variant was kept whatever its predicted impact, because second alleles in MVA1 are often non-coding, and a single *BUB1B* allele also gets its own row as the hypothesis that its partner was missed. The top-ranked genotype is in *BUB1B*: a compound heterozygous pair, phase unknown, of a stop_gained allele (ClinVar: Pathogenic/Likely pathogenic) and a missense allele (not in ClinVar). Per-chromosome B-allele fractions were fitted as a screen for mosaic aneuploidy. All patient-data compute ran in private Kaggle notebooks; the code, the pre-registered rules and a synthetic end-to-end test are in the repository. Limitations: no parental data (phase only where reads span both variants), no BAM or RNA, and VCF-level calling cannot see most transposon insertions or deep structural variants.

---

## Part B. Detailed methods and results


*[Withheld from the public version: patient-level detail. It appears in the submitted report.]*


### B1. Data, access and governance
- **Dataset:** whole-genome FASTQ (four pairs), one VCF with its index, and a phenotype document with HPO terms. It is shared under WCG IRB protocol #20252010 through gated access on Hugging Face.
- **Reference build:** GRCh38 (no-alt analysis set with decoys); contigs have no `chr` prefix. Variants were called with Sentieon, then filtered with GATK VariantFiltration.
- **Where the data is processed:**
  - The VCF, its index and the phenotype document are downloaded from the gated repository to one laptop. There, gated and derived files stay under `D:\mva-hackathon\`, and nothing gated or derived enters git (a pre-commit guard enforces this).
  - From 2026-10-14, the compute runs in **private Kaggle notebooks** on the entrant's own account: the whole-VCF annotation, ranking and mosaicism statistics. The VCF and its index sit in a private Kaggle dataset that only the account owner can see. Each notebook keeps only small outputs: the ranked file, the rows for the ten MVA genes and summary statistics.
  - The hackathon rules allow cloud instances and notebooks, provided everything is deleted within 30 days of close.
- **Deletion:** all of it is deleted by 2026-11-23, on the laptop and on Kaggle. No patient data went to any other service.

### B2. Compute environment and VCF quality control
- **Compute:**
  - Patient data is processed in private Kaggle CPU notebooks with bcftools, samtools and htslib from Bioconda (`kaggle/mva_d02/`). Their public inputs come from a separate notebook, `kaggle/mva_refs_build/`:
    - the GRCh38 analysis set, checked against NCBI's md5;
    - ClinVar;
    - the gnomAD slices.

    That notebook also reruns the unit tests and a synthetic end-to-end test before any patient data is touched.
  - The same scripts also run on a laptop (4 threads, 7.8 GB RAM, WSL Ubuntu, the same tool versions), where the pipeline was first built and tested.
  - Public-data precomputation (the *BUB1B* splice map) runs on Kaggle as well.
  - The analysis works from the VCF, plus targeted read-level checks on FASTQ pairs processed one at a time.
- **VCF quality control:** header, sample and build checks; left-alignment and multi-allelic splitting; quality filters; and the `chr` prefix added on output.

### B3. Annotation and inheritance-aware filtering
- **Consequences** come from `bcftools csq` given only the 19,299 MANE Select transcripts of Ensembl release 116 (`pipeline/mane_gff.py`), so readthrough transcripts cannot give a gene a spurious consequence (locus dossier §1). Indels are left-aligned and multi-allelic records split first (`bcftools norm`). Only chromosomes 1–22, X and Y are annotated: MANE Select has no mitochondrial, decoy or unplaced-contig transcripts. Variants with no coding or splice consequence are placed by gene interval as intronic, upstream or downstream (within 5 kb), so the *BUB1B* upstream interval is not lost (`tools/annotate_vcf.sh`, `pipeline/annotate_table.py`).
- **ClinVar** classification and review status are attached by position and alleles (release of 2026-09-28).
- **Population frequency:** gnomAD v4.1 joint (exome + genome) frequencies, fetched from the public gnomAD browser API for the ten MVA genes (gene, 60 kb upstream, 5 kb downstream, all variants) and the 84 ACMG SF v3.3 genes (gnomAD's exonic set), and matched on chromosome, position and alleles after normalisation (`tools/fetch_gnomad_slices.py`). Inside these slices the rarity measure is the maximum filtering allele frequency across genetic ancestry groups (fafmax, 95% CI), and a variant gnomAD does not list counts as 0. Elsewhere it is the highest of ClinVar's ExAC, 1000 Genomes and ESP frequencies. The source is recorded per row (`pop_af_source`).
- **Rare-variant filter:** that frequency at most 1%. A missing frequency counts as rare, so the filter can only over-include.
- **Recessive models:** homozygous, and compound heterozygous pairs, with read-backed phase from the PGT/PID tags; pairs without shared phase are scored as phase unknown (`pipeline/comphet.py`). In *BUB1B* every rare variant is kept whatever its predicted impact, because MVA second alleles are often non-coding.
- The whole path was tested end to end on a synthetic VCF built from public ClinVar *BUB1B* alleles (`tools/test_annotate_e2e.py`).

### B4. Phenotype-driven prioritisation
- HPO-based ranking.
- An independent clinical route through Genomics England PanelApp panels, where *BUB1B* is a green, biallelic gene on several phenotype panels (locus dossier §6).
- **The panel gene list was fixed on 2026-10-04, before any access to the data** (`tools/build_panel_genes.py`). It is the union of green genes on the ten PanelApp panels on which *BUB1B* is green and biallelic (locus dossier §6): childhood solid tumours (two panels), familial rhabdomyosarcoma, severe microcephaly, hydrocephalus, fetal anomalies, DDG2P, intellectual disability, and embryonal tumour and sarcoma of possible germline origin. The panels were chosen by MVA's phenotype classes, not from the child's HPO terms. Genes whose mode of inheritance on those panels is monoallelic only are left out, because the ranker scores recessive genotypes. Symbols are mapped to the Ensembl 116 MANE Select names through HGNC IDs. The list has 1,987 genes. About 1,750 of them come only from the three broad developmental panels (DDG2P, intellectual disability, fetal anomalies), so the 0.1 prior is a weak, broad boost, not a diagnosis.

### B5. *BUB1B* locus deep-dive
We seek both alleles by mechanism, not by ClinVar status, because the classic hypomorphic alleles are VUS or conflicting in ClinVar (locus dossier §3.1). The search covers:

- **The promoter and the upstream interval.**
- **Splice-region and deep-intronic variants.**
  - Each is looked up in the public Pangolin map of every possible single-nucleotide change in the gene.
  - Indels are scored individually in the same environment.
- **Transposon insertions and structural variants,** found by k-mer baiting of the locus.
- **Read-backed phase.**
- **Mechanism class and antisense lanes,** applied with the rules fixed in advance in `reports/track2/d08_worksheet.md`.

**Results in model 2.**

- **Both alleles found.** The VCF carries two rare *BUB1B* variants on MANE Select: a truncating allele and a missense allele, both heterozygous PASS calls. 2 further rare *BUB1B* variant(s) are called, none protein-altering.
- **Phase: unknown.** The two variants are about 11 kb apart, beyond read-backed phasing, and no parental data were provided. Compound heterozygosity is assumed, and PM3 is scored at the reduced strength for unknown phase (§B8).
- **Splicing.** In the public Pangolin map of every possible *BUB1B* single-nucleotide change, neither allele nor the 2 other rare variant(s) looked up changes a splice site (largest |Δ| 0.01; a strong call needs 0.5).
- **Upstream interval.** 2 upstream variants are called, all common in gnomAD v4.1 (lowest allele frequency 0.087).
- **Mechanism class** (rules fixed on 2026-09-29, `reports/track2/d08_worksheet.md`): the truncating allele is in the abundance class (nonsense-mediated decay predicted); the missense allele is in the abundance class (provisional).
- **Stability prediction** (rules fixed before the run, `reports/track2/stability_map_rules.md`): validity gate passed (3 of 5 abundance-class controls above both no-defect controls); call for the missense allele: **inconclusive**. For context only, outside the pre-registered call: its predicted destabilisation is within the range of the five tested abundance-class substitutions and above both no-defect controls.
- **Not run:** k-mer baiting for transposon insertions and structural variants, which needs the FASTQ files.


*[Withheld from the public version: patient-level detail. It appears in the submitted report.]*


### B6. Evidence for mosaic aneuploidy
- Per-chromosome B-allele fraction from VCF allele depths (`pipeline/baf_mosaic.py`): the spread of heterozygous allele fractions around 0.5, scored as a robust z against the other autosomes (median and MAD). Autosomes with at least 1,000 heterozygous sites are tested; a chromosome is a lead when z ≥ 3. X and Y are not tested, because their expected balance depends on sex, and decoy contigs are skipped.
- Low-pass depth from alignments of reads taken from the start of the FASTQs.
- Stated detection limits, and the caveat that variegated mosaicism is spread across different chromosomes in different cells.

**Status in model 2:** B-allele fractions from the VCF only (results in §B10). Low-pass depth was not run.

### B7. Differential genes and secondary findings
- Coverage and negative results for the other MVA genes (locus dossier §7).
- The ACMG SF v3.3 screen, reported as counts only in the public version.

**Results in model 2.**

- **The other nine MVA genes are negative.** None carries a ClinVar pathogenic or likely pathogenic call, and no gene has two rare protein-altering alleles or a rare homozygous one.

| Gene | Calls | Median depth of calls | ClinVar P/LP | Rare high or moderate impact |
|---|---|---|---|---|
| *CEP57* | 7 | 41 | 0 | 0 |
| *TRIP13* | 3 | 7 | 0 | 0 |
| *CENATAC* | 20 | 43.5 | 0 | 0 |
| *SLF2* | 75 | 40 | 0 | 0 |
| *SMC5* | 100 | 39.5 | 0 | 0 |
| *MAD1L1* | 910 | 42.0 | 0 | 1 |
| *MAD2L1BP* | 7 | 33 | 0 | 0 |
| *CEP192* | 192 | 42.0 | 0 | 0 |
| *BUB1* | 15 | 42 | 0 | 0 |

- **Coverage caveat:** *TRIP13* has few calls at low depth. A VCF does not record depth where nothing is called, so a cause there is not excluded by this data.
- **ACMG SF v3.3 screen** (84 genes; ClinVar pathogenic or likely pathogenic calls, `pipeline/sf_screen.py`): 0 gene(s) with such a call; 0 reportable or for review; 0 carrier-only or not reportable. A further 13 call(s) in these genes carry only conflicting ClinVar classifications and are not reportable. X-linked findings depend on the sex stated in the phenotype document, which this pipeline never infers.

### B8. ACMG/AMP classification
Points-based classification (Tavtigian 2020, PMID 32720330; `pipeline/acmg_points.py`). Codes were fixed in `tools/draft_track1_report.py` before this report was drafted. PM2 is applied at Supporting (ClinGen SVI). PM3 counts only a pathogenic partner, at Supporting because phase is unknown. PP4 is not applied because MVA has several causal genes. PP3 and PS3 are not applied, because the computational scores are not calibrated for ACMG strength and no functional data exist for these alleles.

| Allele | Evidence | Points | Class |
|---|---|---|---|
| truncating | PVS1, PM2_Supporting | 9 | Likely pathogenic |
| missense | PM2_Supporting, PM3_Supporting | 2 | Uncertain significance |

**What would settle the missense allele.** It stands at 2 points, a variant of uncertain significance. A BUBR1 immunoblot in patient-derived cells showing low protein (PS3, Strong: +4) would bring it to 6 points, likely pathogenic. Confirming the two alleles in trans by parental testing would raise PM3 to Moderate (+1). This is the same first experiment that the Track 2 report puts first (assay 1).

### B9. Submission file
Candidates and their scores (`epcr`) come from `pipeline/track1_rank.py`, whose rule was fixed on 2026-10-02, before any access to the data:

- **Score = gene prior × allele evidence × phase factor.**
  - Gene prior: *BUB1B* 1.0; *CEP57* and *TRIP13* 0.6; genes on the pre-registered PanelApp gene list (Part B §B4) 0.1; any other gene 0.01.
  - Allele evidence: ClinVar pathogenic or likely pathogenic 1.0. Otherwise it depends on predicted impact: high 0.8, moderate 0.3, low 0.15, modifier 0.05. ClinVar benign scores 0.01. A call that fails a VCF filter scores half.
  - Phase: trans 1.0, unknown 0.5, cis 0.05.
- **A single *BUB1B* allele with allele evidence of at least 0.3** also gets its own row, alongside any pair it belongs to. This row is the hypothesis that its partner was missed, because MVA second alleles are often non-coding or structural and short-read calling can miss them. It scores 0.1 of the allele evidence, which ranks it below a seen moderate-impact partner of unknown phase and above a seen low-impact one.
- **EPCR.** The top 10 genotypes are written, each once. EPCR is the score divided by the sum of the written rows' scores plus 0.25.
  - Normalising over the written rows stops thousands of weak candidates from diluting the strong ones. Intronic variants with no known frequency are the main example.
  - The 0.25 is the probability that the cause is in none of the rows, so no row claims more than 0.8.

The file is validated with `tools/validate_track1_csv.py`, which checks the organisers' current field specification: at most 10 rows, `epcr` in (0, 1], `finding_type` either `primary` or `secondary`, and the `chr` prefix.

### B10. Results

*[Withheld from the public version: patient-level detail. It appears in the submitted report.]*


Public summary: The top-ranked genotype is in *BUB1B*: a compound heterozygous pair, phase unknown, of a stop_gained allele (ClinVar: Pathogenic/Likely pathogenic) and a missense allele (not in ClinVar).

| Gene | Variant rows | Rare (AF ≤ 1%) | Rare, high or moderate impact |
|---|---|---|---|
| BUB1B | 17 | 4 | 2 |
| CEP57 | 7 | 1 | 0 |
| TRIP13 | 3 | 0 | 0 |
| CENATAC | 20 | 1 | 0 |
| SLF2 | 75 | 2 | 0 |
| SMC5 | 100 | 5 | 0 |
| MAD1L1 | 910 | 18 | 1 |
| MAD2L1BP | 7 | 1 | 0 |
| CEP192 | 192 | 1 | 0 |
| BUB1 | 15 | 2 | 0 |

B-allele fractions were tested on 22 chromosomes; 3 deviate from the autosomal pattern. A deviation is a lead to check against read depth, not a diagnosis: the smallest chromosomes, which are GC-rich or acrocentric, often show a wider allele balance in short-read data for mapping reasons. Variegated mosaicism can sit below this test's detection limit even when a karyotype shows it.


*[Withheld from the public version: patient-level detail. It appears in the submitted report.]*


Input: 5012204 VCF records; 5149105 annotated (variant, gene) rows; frequency sources {'none': 5102894, 'clinvar': 43179, 'gnomad_v4.1': 3032}.

### B11. Interpretation
- **The top-ranked pair is the confirmed answer.** Model 1, which submitted the same rows, scored 100 rank points and F-max 1.000 on the leaderboard: the clinically confirmed variants are the *BUB1B* pair at rank 1.
- **The pair fits the classic MVA1 architecture:** a truncating null allele with a missense allele in the C-terminal region where every tested substitution made an unstable protein (Suijkerbuijk 2010, PMID 20516114).
- **Ranking and classification answer different questions.** Under ACMG/AMP the missense allele is still a variant of uncertain significance (§B8), yet the gene-first, inheritance-aware ranking put the pair first. A pipeline that kept only ClinVar or ACMG pathogenic calls would have reported the truncating allele alone, worth half credit; this is why the rule keeps every rare *BUB1B* variant whatever its predicted impact (§B9).
- **Mechanism feeds Track 2.** Both alleles are classed as abundance-class, so raising BUBR1 is the coherent therapeutic goal, and assay 1 (BUBR1 protein in patient cells) is the experiment that settles both the classification and the drug hypotheses.

### B12. Limitations
- No parental data: phase is known only where reads span both variants; otherwise it is scored as unknown.
- No BAM, CRAM or RNA: splice predictions are not RNA evidence, and read-level checks need the FASTQ files.
- VCF-level calling misses most transposon insertions and many structural variants.
- gnomAD frequencies cover only the ten MVA genes and the ACMG SF genes; elsewhere ClinVar's frequencies are used, and a missing frequency counts as rare.
- The B-allele-fraction screen cannot see low-level variegated mosaicism.

- The missense allele's mechanism class rests on its position and a computational stability prediction, not on an experiment.
- Low-pass depth and k-mer baiting were not run in this model.

### B13. Reproducibility
- Repository: https://github.com/SanTanBan/mva-hackathon-2026 (commit `91d0cb9`).
- Public references and tests: `kaggle/mva_refs_build/refs_build.py` (GRCh38 analysis set from NCBI, md5-checked; ClinVar 2026-09-28; gnomAD v4.1 slices; bcftools/samtools/htslib 1.24).
- Patient data: `kaggle/mva_d02/d02.py`, run by `tools/kaggle_d02.py` after `tools/fetch_vcf_first.py`, `tools/d01_intake.py` and `tools/kaggle_patient_upload.py`.
- Validation: `tools/validate_track1_csv.py`.

### B14. Privacy and data deletion
No names, dates, locations more specific than "UK (NHS-diagnosed)", or clinical timeline appear. The public version withholds all patient-level variants. The data must be deleted from every environment by 2026-11-23: the laptop, including the Hugging Face download cache, and the private Kaggle dataset and notebooks that processed it. Deletion is then confirmed by email.

---

## Acknowledgement

*"This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible."*

The organisers require this acknowledgement in any public communication arising from the Hackathon.

References: see `reports/track1/locus_dossier.md` for the full, PubMed-checked reference list.
