# Track 2 methods description

**Rare Disease, Real Kid: The MVA Hackathon 2026 · Track 2 · draft 2026-09-30 · public version**

> This file answers the organisers' methods description form (`methods_description_form.xlsx`, sheet "Track 2 methods") question by question, in its order. The answers are copied into that form when the final Track 2 package is staged. Section numbers (§) refer to `reports/track2/track2_report.md`. Italic text marks what is finalised after package D08 classifies the child's alleles.

## 1. Team name
SantanuBanerjee (Hugging Face username). *The display name is confirmed at submission.*

## 2. Approach: from variant and mechanism to candidate medications
1. **Mechanism first (§2).**
   - MVA1 is a BUBR1 dosage disease: a null allele pairs with a hypomorph.
   - BUBR1 is a pseudokinase, so "activating the kinase" is not a goal.
   - The druggable biology is an acetylation–degradation cycle. PCAF acetylation at K250 protects BUBR1. HDAC2/3 deacetylation starts its APC/C-mediated destruction. SIRT2 deacetylation at K668 stabilises it.
   - Three therapeutic goals follow: raise or protect BUBR1, reduce the burden of senescent cells, and intercept cancer within surveillance.
2. **Landscape check (§3).**
   - We searched five trial registries (ClinicalTrials.gov, ISRCTN, EU CTIS, ANZCTR, WHO ICTRP), PubMed, ChEMBL and Open Targets. None has an MVA trial, and no potent BUB1B ligand exists.
   - Open Targets shows two tractability signals that point the wrong way: a kinase-family heuristic, and degrader tractability, which is harmful for a dosage disease.
3. **Hypotheses (§4).** Five ranked hypotheses, each with its strongest evidence and its main weakness: NAD+ precursors, HDAC2/3 inhibition, metformin inside surveillance, a genotype-matched oligonucleotide, and senolytics. The ranking criterion is how cheaply the first decisive experiment in patient-derived cells could disprove the hypothesis.
4. **Exclusions (§5).** Classes that would weaken the checkpoint further are excluded: MPS1/TTK, Aurora B and PLK1 inhibitors, and pan-HDAC inhibitors at high exposure. Growth hormone and aneugenic benzimidazoles are excluded too.
5. **Safety filters (§6, §7).**
   - A six-point chemotherapy-compatibility screen, filled from hand-read US labels.
   - A five-assay decision rule in patient-derived cells, whose veto is any shortening of mitosis.
6. **Automation (§10).** The drug-triage pipeline (`pipeline/drug_triage/`) applies the landscape and safety screens to any gene list from public data.
7. **Genotype specificity (§8.1).**
   - We fixed in advance how each allele's mechanism class and antisense lane change the ranking (`reports/track2/d08_worksheet.md`, `pipeline/allele_mechanism.py`).
   - A public Pangolin map of every possible single-nucleotide change in *BUB1B* shows where an antisense-correctable cryptic splice site could arise (`kaggle/bub1b_splice_map/`).
   - *The ranking is revisited once the child's alleles are classified.*

## 3. Generative AI used (required)
Anthropic Claude Max (consumer plan), used through Claude Code, with model-training opt-out enabled and no feedback or ratings given in sessions involving the data. These are Sage's conditions for consumer plans (HF thread #2, 29 September 2026). The wording was confirmed on 2026-09-30.

## 4. Automated, or manual review and curation?
Both.
- **Automated:**
  - The drug-triage pipeline queries Open Targets, ChEMBL and openFDA labels, then flags each drug. The flags cover paediatric labelling, checkpoint weakening, aneugen class, CYP3A4 perpetrator interactions, kidney and bone marrow.
  - The allele classifier and the splice map are also automated.
- **Manual:** the hypotheses themselves come from a mechanism-led literature review. Every automated flag behind a ranked hypothesis was checked by hand against the full label text.

## 5. Manual literature review and curation
- **Citations:** each PubMed-indexed citation was checked against its PubMed record for title, authors and journal. The specific claim was then checked against the abstract or the open-access full text.
- **Drug labels:** each US prescribing label was read in full. The openFDA set IDs are listed in the report's references.
- **Registries:** each search was run with positive controls, to show that the search itself worked.
- **Error log:** every error found in our own notes is logged in §12 (12 corrections so far). The hand check of labels also found two pipeline errors, which are now fixed and covered by tests.

## 6. Public data only?
Yes. Only public drug and target databases, registries and literature were used. No proprietary data.

## 7. Public data sources
- **Trial registries:** ClinicalTrials.gov, ISRCTN, EU CTIS, ANZCTR, and the WHO ICTRP portal.
- **Literature:** PubMed, and PubMed Central for open-access full text.
- **Drugs and targets:** ChEMBL web services and the Open Targets Platform.
- **Drug labels:** the openFDA drug label API, with DailyMed set IDs.
- **Genes and proteins:** Ensembl REST, NCBI RefSeq (NM_001211.6), UniProt (O60566), Genomics England PanelApp and gnomAD.
- **Splicing:** the Pangolin splice model (Zeng & Li 2022, GPL-3.0).
- **Dates:** query dates are given with each result in the report.

## 8. Proprietary data sources
None.

## 9. How the variant's mechanism was characterised
- **Loss of function by dosage.** Biallelic *BUB1B* variants lower BUBR1 protein.
- **Two allele classes.** Suijkerbuijk 2010 (PMID 20516114) found that MVA mutations either lower BUBR1 abundance or cause a specific functional defect:
  - **Abundance class:** truncating alleles make no transcript, and kinase-domain substitutions misfold and turn over fast.
  - **Functional class:** a stable protein with a specific defect.
- **Pathway disrupted.** The spindle assembly checkpoint, through the mitotic checkpoint complex and APC/C-mediated degradation, regulated by the acetylation cycle above.
- **Downstream consequences.** Chromosome missegregation and mosaic aneuploidy, premature senescence, and embryonal cancer predisposition.
- **For each of the child's alleles,** the class and the antisense lanes are assigned with rules fixed before the genotype was seen (D08 worksheet). *The result is added here after D08.*

## 10. Time or effort
*Finalised at submission.* So far:
- about ten hours of agent-assisted work in four scheduled sessions (27–29 September 2026);
- one live session on 30 September;
- the earlier literature memos;
- human review time, added at submission.

## 11. Method abstract (up to 500 words)

*Draft; updated after D08.*

Mosaic variegated aneuploidy type 1 (MVA1) is caused by biallelic *BUB1B* variants that lower the amount of BUBR1, a pseudokinase at the centre of the spindle assembly checkpoint. We asked which approved medicines could plausibly raise or protect BUBR1, or blunt the consequences of its loss, without harming a child who may need cancer treatment.

We started from mechanism rather than from drug-target databases. That led to an acetylation–degradation cycle as the druggable biology and to three therapeutic goals. A landscape check of five trial registries, PubMed, ChEMBL and Open Targets found no MVA trial, no potent BUBR1 ligand, and two misleading tractability signals.

We rank five hypotheses by how cheaply the first decisive experiment in patient-derived cells could disprove them: NAD+ precursors, HDAC2/3 inhibition, metformin inside tumour surveillance, a genotype-matched oligonucleotide, and senolytics. Each is stated with its main weakness.

Every candidate must pass two filters:
- a six-point chemotherapy-compatibility screen, filled from hand-read US labels;
- a five-assay decision rule whose veto is any shortening of mitosis, the signature of acetylation-deficient BUBR1.

An open pipeline applies the same screens to any checkpoint gene from public data. On its first run it excluded all 28 checkpoint-kinase control compounds automatically.

The rules that re-rank the hypotheses by the child's own alleles were fixed before the genotype was seen. A public map of every possible single-nucleotide change in *BUB1B* shows where an antisense-correctable cryptic splice site could arise.

**Strengths:** mechanism-first reasoning; stated weaknesses; safety filters designed for a cancer-predisposed child; genotype rules fixed in advance; reproducible public-data tooling.

**Limitations:**
- Most mechanistic evidence comes from mice and cell lines.
- There is no exposure-margin analysis, because by design the report gives no doses.
- No RNA data exist to confirm splice effects.
- Unpublished community drug screens may overlap with our candidates.

All outputs are hypotheses for expert follow-up, not treatment advice.

---

## Acknowledgement

*"This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible."*
