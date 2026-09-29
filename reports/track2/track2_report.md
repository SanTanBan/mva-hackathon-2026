# Repurposing hypotheses for BUB1B-related mosaic variegated aneuploidy (MVA1)

**Rare Disease, Real Kid: The MVA Hackathon 2026 · Track 2 · Draft v1.1, public version · literature and database checks run 2026-09-27, registries and drug labels added 2026-09-29**

> **What this is.** A ranked set of research hypotheses for expert follow-up in patient-derived cells. It is not medical advice and is not addressed to the family. It contains no doses, no schedules and no treatment recommendations. Clinical decisions belong to the child's treating team.
>
> **What this is not yet.** This draft was written before we had access to the hackathon dataset, so it reasons from the published genetic architecture of MVA1, not from the child's own two *BUB1B* alleles. Version 2 will classify those alleles by mechanism and re-rank the hypotheses.
>
> **Privacy.** No patient data were used or sent anywhere. The child is described without sex, name, dates or clinical history; the only location given is "UK (NHS-diagnosed)".

Citations are given as first author, year and PubMed ID (PMID). Every PubMed-indexed source was checked against its PubMed record on 2026-09-27, or on 2026-09-29 for the six added then. Database results are cited by accession and query date. The full reference list, with DOIs, is at the end.

---

## 1. Summary

- **MVA1 is a BUBR1 dosage disease.** Patients pair a truncating allele with a hypomorphic one, total BUBR1 protein is low, and re-expressing BUBR1 restores checkpoint function in patient cells (Suijkerbuijk 2010, PMID 20516114). The coherent therapeutic goal is to raise or protect BUBR1 protein. "Activating the kinase" is not a goal: BUBR1 is a pseudokinase whose catalysis is dispensable (Suijkerbuijk 2012a, PMID 22698286).
- **The druggable biology is an acetylation–degradation cycle.** Acetylation at K250 protects BUBR1 from APC/C-mediated destruction, and HDAC2 and HDAC3 remove that mark to start degradation. Mice that cannot acetylate BubR1 phenocopy MVA (Choi 2009, PMID 19407811; Park 2013, PMID 23878276; Park 2017, PMID 28985013). SIRT2 deacetylation at a different lysine, K668, stabilises BubR1, and the NAD+ precursor NMN raised BubR1 levels in mice (North 2014, PMID 24825348).
- **Nothing exists yet (checked 2026-09-27 and 2026-09-29).** ClinicalTrials.gov, ISRCTN, EU CTIS, ANZCTR and the WHO ICTRP portal list no studies for MVA, premature chromatid separation or *BUB1B*. PubMed contains no report of disease-modifying drug therapy in MVA. ChEMBL's most potent BUB1B measurement is about 10 µM. Open Targets lists no clinical-stage agent, and its two positive tractability signals point the wrong way (§3.2).
- **Ranked hypotheses**, by priority for the first decisive experiment in patient cells:
  - (H1) NAD+ precursors, to stabilise BUBR1
  - (H2) HDAC2/3 inhibition, to protect acetylated BUBR1
  - (H3) metformin as chemoprevention embedded in tumour surveillance, following the design of the UK MILI trial and a German trial that includes adolescents
  - (H4) genotype-matched oligonucleotides, conditional on the child's alleles
  - (H5) senolytics

  Each comes with its main weakness stated plainly.
- **Every candidate must pass two safety filters.** The first is a six-point chemotherapy-compatibility screen (§6). The second is a five-assay decision rule in patient cells that vetoes any compound that shortens mitosis (§7). Shortened mitosis is the signature of the acetylation-deficient BubR1 that phenocopies MVA.
- **Do-not-use list:** MPS1/TTK inhibitors, Aurora B inhibitors, PLK1 inhibitors, pan-HDAC inhibitors at high exposure, and growth hormone (§5).

---

## 2. Mechanism: what BUBR1 insufficiency breaks

### 2.1 A dosage disease with a null-plus-hypomorph architecture

- In patients with mutations on both copies of *BUB1B*, a missense mutation pairs with a truncating mutation. Patient cell lines show an impaired mitotic checkpoint, chromosome-alignment defects and low overall BUBR1 levels. The mutations fall into two classes (Suijkerbuijk 2010, PMID 20516114):
  - those that impose specific defects in checkpoint activity or microtubule attachment;
  - those that lower BUBR1 abundance. Truncating alleles produce no transcript, and missense proteins, whose substitutions lie in or near the kinase-like domain, turn over fast.
- **The second allele is often non-coding.** Published examples include:
  - an intergenic variant 44 kb upstream of a *BUB1B* transcription start site, whose causality was proven by TALEN single-base editing (Ochiai 2013, PMID 24344301);
  - a homozygous intronic variant, c.2386-11A>G, that creates a new splice site favoured over the authentic one (Rio Frio 2010, PMID 21190457);
  - a 3.0 kb intronic transposon insertion that genome sequencing missed and RNA-seq plus targeted long-read sequencing found (Kawakami 2023, PMID 38102195).
- **Dose drives phenotype in both directions in mice.** Mice with low BubR1 levels develop progressive aneuploidy and premature-ageing (progeroid) features (Baker 2004, PMID 15208629). Sustained high BubR1 expression reduces tumour formation and extends healthy lifespan (Baker 2013a, PMID 23242215). An allelic series modelling human MVA (*BubR1^H/L1002P*) reproduces cancer predisposition and progeroid features, and shows that allelic effects beyond protein level drive variation between patients (Sieben 2020, PMID 31738183).

**Implication.** Raising total BUBR1, or slowing its destruction, is the most direct disease-modifying goal. Re-expressing BUBR1 rescues checkpoint function in patient cells, so the defect is reversible in principle (Suijkerbuijk 2010).

### 2.2 A pseudokinase: "activate the kinase" is not a goal

Human BUBR1 keeps a kinase-like domain, but its putative catalytic activity is dispensable for error-free chromosome segregation. The residues that would bind ATP in a working kinase instead hold the protein's structure stable (Suijkerbuijk 2012a, PMID 22698286). This explains why missense alleles in this domain act by destabilising the protein. It also makes any "BUBR1 agonist" or "kinase activator" proposal incoherent.

### 2.3 The acetylation–degradation cycle: the druggable biology

| Step | What happens | Evidence |
|---|---|---|
| Upstream signal | When a chromosome is not yet attached, Aurora B phosphorylates BubR1 at S39 (and S16). This step is indispensable for K250 acetylation, for maintaining the fibrous corona and for mitotic checkpoint complex (MCC) stability. | Choi 2026, PMID 41533580 |
| Scaffold | BRCA2 brings the acetyltransferase PCAF to BubR1. Without BRCA2, BubR1 is not acetylated and the checkpoint weakens. | Choi 2012, PMID 22340495 |
| Protective mark | PCAF acetylates BubR1 at K250 in prometaphase, which blocks its ubiquitin-dependent degradation. A mutant that mimics acetylation (K250Q) is stable and delays chromosome segregation. A mutant that cannot be acetylated (K250R) is unstable and speeds up mitosis. | Choi 2009, PMID 19407811 |
| Mark removed, protein destroyed | HDAC2 and HDAC3 bind acetylated BubR1 only in mitosis. This leads to its polyubiquitination, degradation and MCC disassembly. PLK1 is needed for BRCA2 to reach the kinetochore in prometaphase. | Park 2017, PMID 28985013 |
| Genetic phenocopy | Mice carrying one copy of acetylation-deficient BubR1 (K243R, the mouse equivalent of K250R) show early BubR1 degradation, a weaker checkpoint with shorter mitoses, micronuclei and spontaneous tumours. | Park 2013, PMID 23878276 |
| A second, stabilising deacetylation | NAD+-dependent SIRT2 keeps BubR1 deacetylated at K668. Overexpressing SIRT2, or giving mice the NAD+ precursor NMN, raises BubR1 levels in vivo. SIRT2 overexpression extends median lifespan in BubR1 hypomorphic mice. | North 2014, PMID 24825348 |
| Context-dependence | In cells, SIRT2 also deacetylates K250, although total BubR1 did not measurably change (Suematsu 2014, PMID 25285631). In mouse oocytes, the non-acetylatable K243R mutant partly rescued SIRT2-depletion defects (Qiu 2018, PMID 29067790). In post-mitotic neurons, HDAC11 deacetylates BubR1 (Watanabe 2014, PMID 24726361). | as cited |

**The lesson for Track 2.** Removing acetyl groups from BUBR1 has opposite effects depending on the site: deacetylation at K250 starts its destruction, while deacetylation at K668 stabilises it. A hypothesis that touches this cycle has to say which site dominates in the child's cells. The decisive experiment (§7) measures this rather than assuming it.

### 2.4 Beyond the checkpoint

- **Senescence, not aneuploidy alone, drives the progeroid phenotype.**
  - Mice with half doses of both Bub3 and Rae1 carry more aneuploidy than BubR1 hypomorphs yet age more slowly (Baker 2006, PMID 16476774).
  - In BubR1-insufficient mice, deleting p16^Ink4a reduces senescence and premature ageing, while deleting p19^Arf makes both worse (Baker 2008, PMID 18516091).
  - Clearing p16-positive senescent cells delays age-related disorders, and clearing them late in life slows disorders already under way (Baker 2011, PMID 22048312).
  - p53 protects muscle and fat through p21, but p21 promotes lens senescence and cataract (Baker 2013b, PMID 23602569).
- **mTORC1.** In BubR1-mutant mice, susceptibility to muscle wasting correlates with mTORC1 hyperactivity (Sieben 2020, PMID 31738183).
- **PLK1 and cilia.**
  - BUBR1 restrains PLK1 at centrosomes between divisions. Cells from patients with premature chromatid separation (PCS)/MVA show excess PLK1 activity and extra centrosomes, and knocking PLK1 down suppresses the extra centrosomes (Izumi 2009, PMID 19503101).
  - BUBR1 insufficiency impairs formation of primary cilia (Miyamoto 2011, PMID 21389084). In patient cells, a constantly active PLK1–KIF2A pathway causes the ciliary defect (Miyamoto 2015, PMID 25660017).
- **Inflammation.** In MVA caused by *MAD1L1*, single-cell data show interferon and NF-κB activation in both aneuploid and normal-karyotype cells (Villarroya-Beltri 2022, PMID 36322655). The field's leading review names a chronic immune response as a feature of MVA (Malumbres 2024, PMID 39169218). This has not been shown for *BUB1B*.

### 2.5 Three therapeutic goals follow

| Goal | Target | Hypotheses |
|---|---|---|
| A: raise or protect BUBR1 | the cause | H1, H2, H4 |
| B: reduce the burden of senescent cells | the progeroid consequences | H5 |
| C: intercept cancer within surveillance | the embryonal-tumour risk | H3; KIF18A as a research lane (§8.2) |

---

## 3. Landscape check (run 2026-09-27)

### 3.1 Registries, literature and drug-target databases

| Resource | Query | Result | What it means |
|---|---|---|---|
| ClinicalTrials.gov | condition: "mosaic variegated aneuploidy" OR "premature chromatid separation" OR BUB1B | **0 studies** of any type | No trial has ever been registered for MVA. |
| ClinicalTrials.gov | Li-Fraumeni syndrome × metformin | 1 study: NCT01981525 (NCI pilot, phase 1, 26 participants, completed) | The closest chemoprevention precedent on the US registry |
| ClinicalTrials.gov | KIF18A inhibitors (ATX-295, VLS-1488, sovilnesib, AMG 650, MEN2501, GenSci122) | 6 trials, **all in adults with advanced cancer** (for example NCT06799065 and NCT05902988) | Aneuploidy-selective agents exist in the clinic, but only in adult oncology. |
| ISRCTN | "mosaic variegated aneuploidy"; "premature chromatid separation"; "BUB1B" | **0 results** for each | The UK registry is empty too. |
| ISRCTN | MILI, ISRCTN16699730 | Recruiting and ongoing. Sponsor: University of Oxford. Funders: NIHR, the Pantziarka TP53 Trust, Cancer Research UK. Participants aged 16 or over. | A live UK precision-prevention trial built on surveillance in an inherited cancer-predisposition syndrome (Dixon-Zegeye 2024, PMID 38308321) |
| EU CTIS (run 2026-09-29) | "mosaic variegated aneuploidy"; "premature chromatid separation"; "BUB1B" (positive controls "Duchenne" and "Wilms" returned 32 and 8 trials) | **0 results** for each | The EU register is empty too. |
| EU CTIS | "Li-Fraumeni" | 1 trial: 2025-520714-63-00, the LFS-Metformin-Trial. It is a phase IIIb randomised open-label trial of metformin for cancer prevention in **adolescents and adults** with Li-Fraumeni syndrome. Authorised in Germany in August 2025, with recruitment from January 2026. Sponsor: Medizinische Hochschule Hannover (non-commercial). Funder: German Cancer Aid. Primary end point: cancer-free survival. | A metformin chemoprevention trial that names adolescents in its title; the closest age precedent for H3 that we found |
| ANZCTR (run 2026-09-29) | the same three MVA terms (positive controls "Li-Fraumeni" and "Wilms" returned 2 and 30 records) | **0 results** for each | The Australian and New Zealand register is empty too. |
| WHO ICTRP (run 2026-09-29) | "mosaic AND variegated AND aneuploidy"; "premature AND chromatid AND separation"; "BUB1B" | **0 results** for each. The positive control "Li-Fraumeni AND metformin" returned the three trials above (NCT01981525, ISRCTN16699730, CTIS 2025-520714-63-00). | No MVA trial is registered in any registry that feeds the WHO portal. |
| PubMed | "mosaic variegated aneuploidy"[Title/Abstract] AND (treatment OR therapy OR therapeutic) | 11 records. **None reports disease-modifying therapy.** They cover cancer care, transplantation and three MVA2 case reports of poor growth-hormone response. | The repurposing question is open. |
| ChEMBL | target CHEMBL4295998 (BUB1B, UniProt O60566) | 24 activity records. Most potent IC50 values are 9.67 µM and 10.36 µM (pChEMBL 5.01 and 4.99), then ceritinib at 13.08 µM. Dasatinib leaves more than 75% of BUB1B activity at 1 µM. Thermal-shift assays (assay CHEMBL6125447) show some ligands destabilising BUB1B at 20 µM and others stabilising it. | No potent BUB1B ligand exists, and binding BUBR1 can destabilise it. |
| Open Targets | target ENSG00000156970, tractability | No approved or clinical-stage agent in any modality. Small-molecule "High-Quality Ligand" and "Druggable Family" flags are true. All antibody flags are false. Protein-degrader "Literature", "UniProt Ubiquitination" and "Database Ubiquitination" flags are true. | See the two traps below. |

**All five registries have now been searched** (ClinicalTrials.gov, ISRCTN, EU CTIS, ANZCTR and the WHO ICTRP portal). None lists a trial for MVA.

### 3.2 Two tractability traps in Open Targets

1. **"Druggable Family: true" is a false positive for this biology.** The flag is inherited from the kinase annotation. But BUBR1's catalytic activity is dispensable, and its ATP-site residues keep the protein structurally stable (Suijkerbuijk 2012a, PMID 22698286). Inhibiting that pocket is not a therapeutic route, and the ChEMBL thermal-shift data show that ligands binding BUBR1 can destabilise it.
2. **The only other positive signals point in the harmful direction.** Outside the kinase-family heuristic, the only positive flags are for targeted protein degradation, based on ubiquitination evidence. This fits the biology: HDAC2/3 deacetylation leads to ubiquitination and APC/C-dependent destruction (Choi 2009, PMID 19407811; Park 2017, PMID 28985013). But targeted degradation is exactly what a dosage-deficiency disease must avoid, so the strategy has to be stabilising the protein, not engaging it. Antibodies are irrelevant because BUBR1 is inside the cell.

---

## 4. Ranked hypotheses

**Ranking criterion.** Priority goes to what should be tested first in patient-derived cells. Hypotheses that act on BUBR1 levels, the causal variable, and that can be disproved cheaply rank highest. Our earlier research notes ranked metformin second; it moves to third here because it acts on a consequence rather than the cause.

Each item is a hypothesis for expert follow-up, not a proposal to treat.

| Rank | Hypothesis | Acts on | Strongest evidence | Main weakness |
|---|---|---|---|---|
| H1 | NAD+ precursors (NR, NMN) raise BUBR1 through SIRT2 at K668 | the cause (goal A) | NMN raised BubR1 in mice, and SIRT2 overexpression extended hypomorph lifespan (North 2014) | SIRT2 also removes the protective K250 mark, and mouse biochemistry is not patient benefit. |
| H2 | HDAC2/3 inhibition protects acetylated BUBR1 from degradation | the cause (goal A) | HDAC2/3 start BubR1 destruction (Park 2017), and acetylation-deficient mice phenocopy MVA (Park 2013) | HDAC inhibitors can themselves cause aneuploidy and suppress blood counts. |
| H3 | Metformin as chemoprevention embedded in surveillance | a consequence (goal C) | The energy-stress drug AICAR selectively kills aneuploid cells (Tang 2011), and the MILI trial provides a ready design (Dixon-Zegeye 2024) | In a child whose normal tissues are mosaic, selectivity for aneuploid cells may mean selectivity against the child's own tissue. |
| H4 | Genotype-matched oligonucleotide (splice correction or expression up-regulation) | the cause (goal A) | A bespoke antisense drug was built for one patient (Kim 2019), and up-regulation from hypomorphic alleles has been demonstrated (Lim 2020) | Entirely conditional on the child's alleles, and delivery to dividing tissues is unsolved |
| H5 | Senolytics (the dasatinib + quercetin class) | a consequence (goal B) | Genetic clearance of p16-positive cells rescued BubR1 progeroid mice (Baker 2011) | No senolytic drug has been tested in BubR1 models, and senescence may be a brake on tumours in this child. |

### H1: NAD+ precursors to stabilise BUBR1

- **Mechanism.** As NAD+ declines, SIRT2 can no longer keep BubR1 deacetylated at K668. Overexpressing SIRT2, or treating mice with NMN, raises BubR1 levels in vivo, and SIRT2 overexpression increases median lifespan in BubR1^H/H mice (North 2014, PMID 24825348).
- **The benefit looks specific to low BubR1, which is supported rather than assumed.** SIRT2 overexpression does not extend lifespan in normal mice (Wu 2023, PMID 38009412). That fits a benefit conditional on low BubR1, which is this child's situation.
- **Precedent in children.** Nicotinamide riboside has been given open-label to children with ataxia-telangiectasia, another genome-instability disorder, with improved ataxia scores (Panvino 2026, PMID 42387837).
- **Weaknesses.**
  1. SIRT2 also deacetylates K250, the mark that protects BubR1 from APC/C (Suematsu 2014, PMID 25285631). In oocytes, the non-acetylatable K243R state was the protective one (Qiu 2018, PMID 29067790). Which site dominates in MVA1 cells is unknown.
  2. NMN produced the biochemical change, but the lifespan effect came from the SIRT2 transgene. These are different experiments.
  3. In ataxia-telangiectasia, promising open-label results have repeatedly been followed by negative phase 3 trials (Panvino 2026).
  4. NAD+ supports dividing cells, so any use in a cancer-predisposed child would need explicit tumour surveillance. This is a reasoned concern that has not been tested.
- **Decisive experiment.** In patient-derived cells, does NR or NMN raise BUBR1, and how does it change K250 versus K668 acetylation, without shortening mitosis? (§7)

### H2: HDAC2/3 inhibition to protect acetylated BUBR1

- **Mechanism.** HDAC2/3 bind acetylated BubR1 only in mitosis and start its polyubiquitination and destruction (Park 2017, PMID 28985013). Losing BubR1 acetylation weakens the checkpoint, shortens mitosis and causes tumours in mice (Park 2013, PMID 23878276). Blocking deacetylation is the most direct way to keep more BUBR1 through mitosis.
- **Precedent in children.**
  - Givinostat, a class I/II HDAC inhibitor, is approved for Duchenne muscular dystrophy in patients aged 6 and over (Lamb 2024, PMID 38967716). Long-term open-label use in DMD has been reported (McDonald 2025, PMID 40830818).
  - Valproate, an established antiseizure medicine, inhibits HDACs at therapeutically relevant concentrations (Michaelis 2004, PMID 14978230). Seizures were reported in 45.1% of published MVA1 cases (Pavone 2022, PMID 35804254).
- **Weaknesses.**
  1. The HDAC inhibitor trichostatin A both breaks chromosomes and causes aneuploidy in human lymphoblasts, through missegregation (Olaharski 2006, PMID 16857700). A pan-HDAC inhibitor could therefore worsen the defect it is meant to treat. Isoform selectivity and exposure are the substance of this hypothesis, not refinements.
  2. The US label for givinostat reports thrombocytopenia in 33% of treated patients in the pivotal study and in none on placebo, with anaemia and neutropenia also listed (DUVYZAT prescribing information). A published commentary gives about 20% (Anjum 2025, PMID 40330090). Either figure conflicts directly with the bone-marrow criterion in §6.
  3. No approved agent is selective for HDAC2/3.
- **Decisive experiment.** The same five assays (§7). The micronucleus and PCS readout is the likely point of failure.

### H3: Metformin as chemoprevention embedded in surveillance

- **Mechanism.** In the key screen for aneuploidy-selective compounds, the energy-stress agent AICAR caused p53-mediated death of mouse cells carrying an extra chromosome (Tang 2011, PMID 21315436). Metformin is the clinically available drug that induces energy stress. Aneuploidy causes protein-folding stress, shown in yeast (Oromendia 2012, PMID 23222101).
- **A trial design already exists.** MILI randomises adults with Li-Fraumeni syndrome to metformin plus annual MRI surveillance, or to surveillance alone. It is built into NHS standard-of-care surveillance, and its primary end point is cancer-free survival. Its no-placebo design was chosen after consulting the patient community (Dixon-Zegeye 2024, PMID 38308321; ISRCTN16699730, recruiting). A German phase IIIb trial now tests metformin for cancer prevention in adolescents as well as adults with Li-Fraumeni syndrome, with the same primary end point (LFS-Metformin-Trial, EU CTIS 2025-520714-63-00; §3.1). An MVA1 version would sit on top of the renal ultrasound surveillance already recommended for MVA (§9).
- **Weaknesses.**
  1. **The mosaic objection.** AICAR selects against aneuploid cells, not tumour cells. A child with MVA1 has normal tissues that are constitutionally mosaic-aneuploid, so selectivity for aneuploid cells may mean selectivity against the child's own tissue.
  2. The selectivity depended on p53 (Tang 2011), and p53 mutations were frequent in tumours from acetylation-deficient BubR1 mice (Park 2013, PMID 23878276).
  3. AICAR is not metformin.
  4. Metformin is mainly cleared by the kidneys, and chronic kidney disease is a contraindication (Dixon-Zegeye 2024). Renal function therefore gates it. That matters because the kidney is the surveillance organ, and surgery can remove kidney tissue (Ochiai 2019, PMID 31081598).
  5. Gastrointestinal side effects are common (Dixon-Zegeye 2024) and unwelcome in a child with growth failure.
- **Why it is still ranked.** It is the lowest-risk and most deliverable hypothesis, and the only one with a ready-made UK trial design.

### H4: Genotype-matched oligonucleotide (conditional)

Details are in §8.1. This ranks fourth only because eligibility depends on alleles we have not yet seen. If the child carries a splice-altering or otherwise antisense-amenable allele, H4 becomes the highest-potential hypothesis in this report.

### H5: Senolytics

- **Mechanism and evidence.**
  - Clearing senescent cells is the best-validated rescue in this exact model. The INK-ATTAC clearance system was demonstrated in BubR1 progeroid mice, where clearance delayed and reduced age-related disease (Baker 2011, PMID 22048312). Deleting p16 rescued the mice, while deleting p19^Arf made them worse (Baker 2008, PMID 18516091).
  - Dasatinib plus quercetin reduced senescent cells in aged, irradiated and Ercc1-progeroid mice (Zhu 2015, PMID 25754370). It was feasible in adult pilot studies in idiopathic pulmonary fibrosis (Justice 2019, PMID 30616998; Nambiar 2023, PMID 36857968).
  - In ChEMBL, dasatinib does not inhibit BUB1B at 1 µM (§3.1).
- **Weaknesses.**
  1. The BubR1 evidence is genetic, not pharmacological. No senolytic drug has been tested in BubR1 mice.
  2. All human data come from adults.
  3. **The cancer objection.** In a child whose cells keep missegregating chromosomes, p16-driven senescence may be a main brake on aneuploid clones, and clearing it could plausibly raise tumour risk. This is a reasoned concern.
  4. The senolytic navitoclax causes on-target, dose-limiting thrombocytopenia (Mahadevan 2025, PMID 41225624), which rules it out here.
- **Best use.** The experiment the field has not done: testing a senolytic drug in the BubR1^H/L1002P mouse (Sieben 2020, PMID 31738183), rather than any use in a child.

### Considered and not ranked

- **Rapamycin analogues.** mTORC1 hyperactivity gives a rationale (Sieben 2020), but the class is immunosuppressive and a poor fit for a child with growth failure. Deferred until there is patient-cell evidence of mTORC1 activation in MVA1.
- **HSP90 inhibitors.** 17-AAG was one of the compounds selective for aneuploid cells (Tang 2011). Our research notes said no HSP90 inhibitor is approved; that is wrong. Pimitespib is approved in Japan for advanced gastrointestinal stromal tumours (Muraoka 2025, PMID 40167031). It is an oncology drug with frequent gastrointestinal and eye-related adverse reactions (Kurokawa 2025, PMID 41264161), so it is excluded.
- **Hydroxychloroquine.** Chloroquine was the third compound selective for aneuploid cells (Tang 2011), but there are no in vivo chemoprevention data. Excluded in this version.
- **Direct APC/C or proteasome inhibition.** Mechanistically this is on target, since APC/C destroys BubR1 (Choi 2009), but there is no clinical path for long-term use in children. Excluded.

---

## 5. Do-not-use list (exclusion criteria for this hypothesis set)

These classes are excluded because they would further weaken a checkpoint that is already weak, or because they carry a consensus contraindication. This is a research exclusion list, not clinical guidance. The treating team's decisions, for example about cancer care, are outside its scope.

| Class | Why it is excluded | Evidence |
|---|---|---|
| **MPS1/TTK inhibitors** (clinical candidates exist, such as BOS172722) | MPS1 is one of the main components of the spindle assembly checkpoint, and inhibiting it disables the checkpoint. | Woodward 2018, PMID 30199249 |
| **Aurora B inhibitors** | Aurora B phosphorylation of BubR1 is indispensable for the protective K250 acetylation, for corona maintenance and for MCC stability. | Choi 2026, PMID 41533580 |
| **PLK1 inhibitors** | PLK1 is needed for BRCA2 to reach the kinetochore, and so for BubR1 acetylation (Park 2017, PMID 28985013). It also phosphorylates the BUBR1 KARD region that recruits PP2A-B56 to stabilise attachments (Suijkerbuijk 2012b, PMID 23079597). PLK1 inhibition does restore cilia in PCS patient cells (Miyamoto 2015, PMID 25660017). That points research toward the downstream KIF2A step, not toward PLK1 inhibitor drugs. | as cited |
| **Pan-HDAC inhibitors at high exposure** | Trichostatin A causes aneuploidy in human lymphoblasts through missegregation. H2 survives only if an isoform-selective agent passes the aneuploidy veto at exposures that raise BUBR1. | Olaharski 2006, PMID 16857700 |
| **Growth hormone** | A 2022 consensus statement says GH is generally contraindicated in cancer-predisposition syndromes, and considered only cautiously in select patients (Boguszewski 2022, PMID 35319491). Tumour growth on GH has been reported in Noonan syndrome (McWilliams 2016, PMID 26377682). Three MVA2 (*CEP57*) case reports describe poor response to GH (Feng 2022, PMID 35434947; Hübner 2022, PMID 36588761; Viudes 2026, PMID 41700350). Growth failure in MVA1 starts before birth and is plausibly intrinsic to the cells (Hanks 2004, PMID 15475955; Baker 2004, PMID 15208629). | as cited |

**Relative caution (a class inference, not yet a finding).** Drugs that bind tubulin and poison the spindle cause aneuploidy in human lymphocytes (Guccini 2012, PMID 22952150), and MVA fibroblasts fail to arrest in colcemid (Matsuura 2000, PMID 10877982). How such drugs behave in BUBR1-insufficient cells is an untested research question.

**Benzimidazole anthelmintics** (for example mebendazole, a popular oncology repurposing candidate) now have specific evidence:

- Albendazole and its active metabolite caused chromosome non-disjunction in cultured human lymphocytes, at lower concentrations than those that caused micronuclei (Ramírez 2006, PMID 17055335).
- The mebendazole US label states that it is aneugenic in mammalian cells in vitro, and that its in vivo mouse micronucleus result suggested aneugenicity.
- In a small before-and-after study of 14 children treated with albendazole, lymphocyte micronucleus and sister-chromatid-exchange frequencies rose (Oztas 2007, PMID 17331781). The albendazole US label, by contrast, reports a negative in vivo mouse micronucleus test.
- Non-disjunction caused by spindle inhibitors, mebendazole among them, appears to have a concentration threshold (Elhajouji 1997, PMID 9175637), so exposure matters. Whether that threshold is lower in BUBR1-insufficient cells is untested.

This class is therefore not a candidate here, and the R06 pipeline excludes it automatically. Whether an anthelmintic is needed to treat an infection is a clinical decision outside this report's scope.

---

## 6. Six-point chemotherapy-compatibility screen

A child with MVA1 may need cancer treatment at short notice. Genomic-instability disorders often show heightened sensitivity to chemotherapy and ionising radiation (Nakano 2024, PMID 39264246). In one published case of biallelic *BUB1B* mutations, a severe reaction to dactinomycin stopped all further chemotherapy (Ochiai 2019, PMID 31081598). Any hypothesis therefore has to stay compatible with future cancer care. It must meet six criteria:

1. **It must not weaken the checkpoint further.** Cancer risk tracks the severity of checkpoint failure (Yost 2017, PMID 28553959), and losing BubR1 acetylation alone is enough to cause tumours in mice (Park 2013).
2. **It must not have a major CYP3A4 interaction with vincristine-based regimens.** A CYP3A4 inhibitor has been reported to worsen vincristine neurotoxicity in a child (Eiden 2009, PMID 19346885).
3. **It must not harm the kidneys and must remain usable with reduced kidney function.** The kidney is the surveillance target, and nephrectomy has been needed (Ochiai 2019).
4. **It must not suppress the bone marrow.** Cytopenia and myelodysplasia occur in MVA1 (Laberko 2019, PMID 31053147).
5. **It must show no signal of causing aneuploidy or chromosome breaks** (Olaharski 2006; Guccini 2012).
6. **It must be quick to stop:** within days if cancer treatment begins.

| Hypothesis | 1 Checkpoint | 2 CYP3A4 | 3 Kidney | 4 Bone marrow | 5 Genotoxicity | 6 Quick to stop |
|---|---|---|---|---|---|---|
| H1 NAD+ precursors | no flag found | not assessed: NR and NMN have no US prescription label | note: kidney was among the target organs at high doses in a 90-day rat study (Conze 2016) | note: nicotinamide, a related NAD+ precursor, raised the risk of thrombocytopenia in haemodialysis patients (Ketteler 2023); we found no such report for NR itself | no flag found: NR was negative in the Ames, chromosome-aberration and micronucleus tests (Conze 2016) | yes |
| H2 HDAC2/3 inhibition | flag: class causes aneuploidy at high exposure | note: givinostat is a weak intestinal CYP3A4 inhibitor; the valproate label names no CYP3A4 effect | no flag found: givinostat is not cleared mainly by the kidney, and valproate clearance falls only slightly in renal failure. See the OCT2 note below. | flag: dose-related thrombocytopenia, with anaemia and neutropenia (givinostat label); dose-related thrombocytopenia (valproate label) | flag: class effect (Olaharski 2006); givinostat was positive in the Ames test and negative for chromosome aberrations (label) | yes |
| H3 Metformin | no flag found | no flag found: not metabolised by the liver; excreted unchanged in urine (label) | flag: renal gating (Dixon-Zegeye 2024); contraindicated in severe renal impairment (label) | note: can lower vitamin B12, possibly with anaemia, reversibly (label) | no flag found: negative Ames, mouse-lymphoma, human-lymphocyte chromosome-aberration and mouse micronucleus tests (label) | yes |
| H4 Oligonucleotide (nusinersen label as the class example) | passes: raises BUBR1 by design | passes: nusinersen is not a substrate, inhibitor or inducer of CYP450 enzymes | flag: class warning of renal toxicity, including glomerulonephritis; urine protein testing is required before every administration | flag: thrombocytopenia and coagulation abnormalities | no flag found in the Ames, chromosome-aberration and micronucleus tests. Note: vascular tumours increased in mice at the highest dose in a 2-year study. | flag: long-acting by design |
| H5 Senolytics | passes: dasatinib spares BUB1B (ChEMBL) | note: dasatinib is a time-dependent CYP3A4 inhibitor; it has been combined with multi-agent chemotherapy in 81 children with Ph+ ALL (label) | note: renal failure is a rare adverse reaction to dasatinib, under 1% (label) | flag: severe thrombocytopenia, neutropenia and anaemia on continuous dasatinib (label); navitoclax excluded (Mahadevan 2025) | flag: dasatinib was clastogenic in hamster cells in vitro but negative in the rat micronucleus test (label). Quercetin is mutagenic in vitro; a review found no genotoxicity in vivo (Harwood 2007). | yes: adult pilots used short courses; dasatinib's half-life is 3 to 5 hours (label) |

"No flag found" means nothing surfaced in the sources we reviewed; it is not proof of safety. "Label" means the current US prescribing information, read in full by hand on 2026-09-29 (openFDA set IDs are listed under References). For H4 we use nusinersen as the class example, because the bespoke precedent, milasen (Kim 2019), has no label. Oligonucleotide safety depends on chemistry and route, so a genotype-matched drug would need its own assessment.

**An interaction between hypotheses.** Givinostat weakly inhibits OCT2, a kidney transporter; its label says the effect was seen in creatinine measurements in clinical trials. The metformin label warns that OCT2/MATE inhibitors can raise metformin exposure. So H2 and H3 should not be tested together without planning for this. A creatinine rise from givinostat could also blur the kidney monitoring that H3 relies on.

**How the screen was run.** The drug-triage pipeline in this repository (work package R06, `pipeline/drug_triage/`) screens these drugs, and every Open Targets candidate for a gene, against checkpoint, aneugen, CYP3A4, kidney, bone-marrow and paediatric-label criteria from public data. Its output is in `reports/track2/drug_triage/drug_triage_demo.md`. The table above started from those text-mined flags, and each one was then checked against the label text. The check found two pipeline errors, which are now fixed and covered by tests. Givinostat had been recorded as a CYP3A4 substrate when its label calls it a weak inhibitor. Dasatinib's time-dependent inhibition had been missed.

---

## 7. Five-assay decision rule, with the mitotic-duration veto

The assays run in patient-derived lymphoblastoid cells or fibroblasts, the systems the MVA literature already uses (Matsuura 2000, PMID 10877982; Suijkerbuijk 2010, PMID 20516114; Miyamoto 2015, PMID 25660017). Each candidate is tested at clinically achievable free concentrations. BUBR1 re-expression serves as the positive control (Suijkerbuijk 2010).

| # | Assay | Pass | Rationale |
|---|---|---|---|
| 1 | BUBR1 immunoblot: total protein, plus K250 and K668 acetylation where antibodies allow | ↑ | BUBR1 level is the causal variable (Suijkerbuijk 2010; Baker 2013a). This assay also answers the H1 question of which lysine matters. |
| 2 | Live-cell mitotic duration, from nuclear envelope breakdown to anaphase | Not shortened. **Shortening is a veto.** | Shorter mitosis is the signature of acetylation-deficient BubR1, which phenocopies MVA (Choi 2009, PMID 19407811; Park 2013, PMID 23878276). |
| 3 | Checkpoint maintenance when the spindle is poisoned | ↑ or unchanged | MVA cells fail to arrest in colcemid (Matsuura 2000), and checkpoint maintenance fails with acetylation-deficient BubR1 (Park 2013). |
| 4 | Micronuclei, plus premature chromatid separation and aneuploidy | ↓ or unchanged | Catches compounds that cause aneuploidy (Olaharski 2006, PMID 16857700) |
| 5 | Senescence markers (SA-β-gal, p16) | ↓ or unchanged | Readout for goal B (Baker 2008, PMID 18516091; Baker 2011, PMID 22048312) |

**Decision rule.** Advance a compound only if assay 1 rises, assay 2 does not shorten and assay 4 does not worsen. **Any shortening of mitosis is a veto, whatever the other results.** The rule does not depend on the gene and transfers unchanged to other checkpoint-gene disorders (§10).

---

## 8. Two lanes beyond small-molecule repurposing

### 8.1 The oligonucleotide lane

- **Why it fits MVA1.** The disease is defined by dosage, and second alleles are often splice-altering or non-coding (Rio Frio 2010, PMID 21190457; Ochiai 2013, PMID 24344301; Kawakami 2023, PMID 38102195). One published MVA1 genotype pairs a frameshift with an intronic splice-region variant, c.1288+5G>A (Laberko 2019, PMID 31053147). That is exactly the kind of variant splice-correcting oligonucleotides are designed for.
- **Base rate.** One paediatric hospital applied the N=1 Collaborative's VAEAOT eligibility guidelines to all diagnostic variants found in 532 patients. Of those patients, 26 (4.9%), carrying 25 distinct variants, were eligible or likely eligible for an antisense oligonucleotide (ASO). After disease and delivery factors were considered, 11 diagnoses remained candidates (Cheerie 2025, PMID 41090344).
- **Two sub-types.**
  - *Splice correction* for alleles that create cryptic splice sites. The precedent is milasen, where tests in the patient's own cells supported starting an N-of-1 study within one year of first contact (Kim 2019, PMID 31597037).
  - *Expression up-regulation*: blocking naturally occurring non-productive splicing can raise output from normal or hypomorphic alleles (Lim 2020, PMID 32647108). This is the relevant case in MVA1, where one allele is null and the other is intact but produces too little.
- **Annotation check (Ensembl, 2026-09-29).** Besides the MANE Select transcript, Ensembl annotates 12 *BUB1B* transcripts as nonsense-mediated-decay (NMD) targets and 5 as retained-intron transcripts. We compared each with MANE Select, using the exon numbering of our locus dossier:
  - The most recurrent non-productive event is **skipping of exon 5** (197 nt), which shifts the reading frame; three NMD transcripts carry it.
  - Other frameshifting events are single skips of exon 9 (230 nt) and exon 20 (143 nt), a joint skip of exons 7 and 8, and an extra 124-nt exon between exons 2 and 3 that occurs only in an NMD transcript, a candidate poison exon.
  - The published splice-region variant c.1288+5G>A lies at the exon 9 donor site. It would be expected to promote the same exon 9 skip, but that is an inference, not a measurement.

  These are candidate targets for up-regulation, not validated ones. An annotation shows that an event exists, not how often it happens in the child's cells. The decisive measurement is RNA sequencing of patient-derived cells with and without NMD inhibition, to show how much *BUB1B* output each event removes (the approach of Lim 2020). An ASO that suppresses one of these events would raise output from both alleles; for MVA1, that means more protein from the hypomorphic allele. Raising BubR1 above normal levels was protective, not harmful, in mice (Baker 2013a, PMID 23242215).
- **Version 2 deliverable (package D08).** Classify each of the child's alleles by mechanism, as functional-class or abundance-class (Suijkerbuijk 2010), and by VAEAOT category. Then re-rank H4.
- **Weaknesses.** This is bespoke drug development, not repurposing. Oligonucleotides are long-acting by design, a poor fit for criterion 6. MVA1 affects many dividing tissues, and delivery beyond a few organs is unsolved.

### 8.2 KIF18A: an aneuploidy-selective target outside the checkpoint

- **Evidence.**
  - Across about 1,000 human cancer cell lines, aneuploid cells were more sensitive to disrupting core checkpoint components. They became increasingly sensitive to checkpoint inhibition over time, and were particularly vulnerable to losing KIF18A (Cohen-Sharir 2021, PMID 33505028).
  - A selective KIF18A inhibitor causes chromosome-alignment defects, mitotic accumulation and cell death, specifically in chromosomally unstable cancer cells (Phillips 2025, PMID 39747049).
  - Six KIF18A inhibitors are in adult cancer trials (§3.1).
- **Why it is a research lane, not a candidate.**
  1. The mosaic objection applies in full.
  2. The inhibitor kills by holding cells in mitosis, and MVA cells arrest poorly (Matsuura 2000, PMID 10877982). How it behaves in BUBR1-insufficient cells is unpredictable.
  3. There are no data in children.
- **The useful question.** Do BUBR1-insufficient patient cells share the KIF18A dependence of chromosomally unstable cancer cells? A tool-compound experiment with the §7 readouts would answer it. The answer matters for understanding MVA-associated tumours.

### 8.3 Measure the inflammation before targeting it

MVA7 cells show interferon and NF-κB activation (Villarroya-Beltri 2022, PMID 36322655). The obvious explanation, that micronuclei activate the cGAS–STING pathway, is disputed: live-cell reporter experiments found micronuclei are not a potent trigger of that pathway (Sato 2024, PMID 38307626). The right first step is to measure an interferon-stimulated gene signature in *BUB1B* patient cells, not to propose an immune-modulating drug.

---

## 9. Impact: fitting the care the child already receives

- **Surveillance is the standard of care.**
  - The AACR guideline update (Nakano 2024, PMID 39264246) reports a recommendation of renal ultrasound every 3 months from birth until age 7 for all MVA conditions, including those with an unknown genetic cause. It also advises avoiding radiation exposure and recommends regular clinical review for rhabdomyosarcoma and other cancers.
  - The same interval applies in North America to every syndrome with at least a 1% Wilms tumour risk (Kalish 2017, PMID 28674120).
  - SIOP-Europe gives the corresponding European recommendations (Hol 2021, PMID 34134020).
- **Cancer risk, stated honestly.**
  - In a pooled early series, 7 of 10 infants developed proven or probable Wilms tumour, diagnosed at 2–16 months (Kajii 2001, PMID 11746029).
  - Malignancy was reported in 38.7% of 31 published MVA1 cases (Pavone 2022, PMID 35804254).
  - Both figures are case-series proportions subject to strong ascertainment bias. The guideline authors state that the rarity of MVA prevents accurate risk estimates (Nakano 2024). We quote no other risk figure.
- **A realistic path.**
  1. Patient-derived cells and the §7 rule come first.
  2. If a licensed drug passes, a formal N-of-1 design would follow, reported to CENT standards (Vohra 2015, PMID 25976398). A cytogenetic biomarker (micronuclei, aneuploidy, PCS rate) would be the pharmacodynamic end point.
  3. MILI shows the UK research system will run a drug-prevention trial inside routine surveillance for an inherited cancer-predisposition syndrome (Dixon-Zegeye 2024, PMID 38308321).
- **Benefit to the wider community.** SIOP-Europe asks that every patient with a cancer-predisposition syndrome be entered in a national or international registry, so that future risk estimates improve (Hol 2021). The same patient-cell rule applies to every MVA subtype.

---

## 10. Scalability

The mechanism map (§2), the six-point screen (§6) and the five-assay rule (§7) do not depend on the gene. They apply unchanged to:

- *TRIP13* (MVA3), where patient cells show severe checkpoint failure (Yost 2017, PMID 28553959);
- *MAD1L1* (MVA7) (Villarroya-Beltri 2022);
- *CEP57* (MVA2), where the checkpoint defect is minimal (Yost 2017), so the ranking would change.

Work package R06 automates the landscape and screening steps over public data (ChEMBL, Open Targets, openFDA drug labels). It runs on *BUB1B*, *TRIP13* and *CEP57* from one configuration file, and adding a gene is one line. In the *BUB1B* demo, all 28 negative-control compounds (TTK, Aurora B and PLK1 inhibitors) were excluded automatically, as were 11 of the 12 HDAC2/3 candidates. That matches the do-not-use list (§5) without hand curation. The twelfth, mocetinostat, is held as "insufficient data": its recorded mechanism does not meet the broad-HDAC rule and it has no US label.

---

## 11. Limitations

- **Genotype-agnostic.** This version was written before we had dataset access. The child's alleles may change the ranking, H4 especially.
- **Registries.** Five registries were searched by keyword (§3.1). A trial registered under a different disease name, for example a generic chromosomal-instability or cancer-predisposition label, could still be missed.
- **Model systems.** Most mechanistic evidence comes from mice, cancer cell lines or yeast. Direct data from human MVA1 cells are few.
- **No exposure-margin analysis.** By design this report gives no doses. Comparing active concentrations in cells with achievable exposure in patients is left to expert follow-up.
- **Possible overlap with community research.** Public sources do not list compounds already screened by MVA community drug-repurposing projects, so duplication cannot be ruled out.
- **"No flag found" in §6** reflects the sources reviewed, not proof of safety.

---

## 12. Corrections to our own research notes (audit trail)

1. The *CEP57* discovery paper is Snape et al. 2011 (PMID 21552266), not "Snijders Blok".
2. We removed an unsourced figure of "about a one-in-three risk of embryonal cancer before age 5" and an unsourced surveillance interval. They are replaced by Nakano 2024 and the case-series figures in §9.
3. PMID 29067790 is by Qiu et al., not "Zhang".
4. The antisense-eligibility base rate is 4.9% of 532 **patients**, not of variants (Cheerie 2025).
5. According to the ISRCTN record, MILI's funders are NIHR, the Pantziarka TP53 Trust and Cancer Research UK.
6. Kajii 2001's figure of 7 in 10 infants with Wilms tumour pools five new infants with five previously reported ones.
7. An HSP90 inhibitor is approved (pimitespib, in Japan). Our notes said otherwise.
8. Oromendia 2012 is a yeast study and is now described as one.
9. The mouse "~10% BubR1" figure is not in the Baker 2004 abstract, so this report says "low levels".
10. The language is sex-neutral throughout, because none of the public sources we rely on states the child's sex.
11. We had quoted a givinostat thrombocytopenia rate of about 20% from a commentary (Anjum 2025). The US label reports 33% in the pivotal study, against none on placebo, so the label figure now leads (§4, H2).

---

## 13. Governance, privacy and LLM use

- **The data.** The hackathon dataset is a real child's genome, shared under WCG IRB protocol #20252010. This public draft uses no patient data (no variants, genotypes, reads or phenotype text) and nothing that could identify the child.
- **The sources.** All sources are public. Each PubMed-indexed citation was checked on 2026-09-27 against its PubMed record: title, authors and journal, with the specific claim checked against the abstract or open-access full text. The six added on 2026-09-29 (Conze 2016, Elhajouji 1997, Harwood 2007, Ketteler 2023, Oztas 2007, Ramírez 2006) were checked the same way. Drug-label statements were read from the current US prescribing information. Registry and database results are cited by accession and query date.
- **LLM use.** An LLM agent (Claude, Anthropic) helped with literature synthesis and drafting, using public sources only. *The provider, plan or tier, relevant setting and data-handling terms still need to be recorded here; the organisers require this, and the confirmed wording is pending.*

---

## References

- Anjum AF et al. 2025. Unleashing the potential of givinostat: a novel therapy for Duchenne muscular dystrophy. *Curr Ther Res Clin Exp* 102:100787. PMID 40330090. https://doi.org/10.1016/j.curtheres.2025.100787
- Baker DJ et al. 2004. BubR1 insufficiency causes early onset of aging-associated phenotypes and infertility in mice. *Nat Genet* 36:744-9. PMID 15208629. https://doi.org/10.1038/ng1382
- Baker DJ et al. 2006. Early aging-associated phenotypes in Bub3/Rae1 haploinsufficient mice. *J Cell Biol* 172:529-40. PMID 16476774. https://doi.org/10.1083/jcb.200507081
- Baker DJ et al. 2008. Opposing roles for p16Ink4a and p19Arf in senescence and ageing caused by BubR1 insufficiency. *Nat Cell Biol* 10:825-36. PMID 18516091. https://doi.org/10.1038/ncb1744
- Baker DJ et al. 2011. Clearance of p16Ink4a-positive senescent cells delays ageing-associated disorders. *Nature* 479:232-6. PMID 22048312. https://doi.org/10.1038/nature10600
- Baker DJ et al. 2013a. Increased expression of BubR1 protects against aneuploidy and cancer and extends healthy lifespan. *Nat Cell Biol* 15:96-102. PMID 23242215. https://doi.org/10.1038/ncb2643
- Baker DJ et al. 2013b. p21 both attenuates and drives senescence and aging in BubR1 progeroid mice. *Cell Rep* 3:1164-74. PMID 23602569. https://doi.org/10.1016/j.celrep.2013.03.028
- Boguszewski MCS et al. 2022. Safety of growth hormone replacement in survivors of cancer and intracranial and pituitary tumours: a consensus statement. *Eur J Endocrinol* 186:P35-P52. PMID 35319491. https://doi.org/10.1530/EJE-21-1186
- Cheerie D et al. 2025. Screening rare genetic diagnoses for amenability to bespoke antisense oligonucleotide therapy development: A retrospective cohort study. *Genet Med* 28:101597. PMID 41090344. https://doi.org/10.1016/j.gim.2025.101597
- Choi E et al. 2009. BubR1 acetylation at prometaphase is required for modulating APC/C activity and timing of mitosis. *EMBO J* 28:2077-89. PMID 19407811. https://doi.org/10.1038/emboj.2009.123
- Choi E et al. 2012. BRCA2 fine-tunes the spindle assembly checkpoint through reinforcement of BubR1 acetylation. *Dev Cell* 22:295-308. PMID 22340495. https://doi.org/10.1016/j.devcel.2012.01.009
- Choi SY et al. 2026. Linking kinetochore attachment to checkpoint control: the role of Aurora B in BubR1 acetylation. *Nucleic Acids Res* 54. PMID 41533580. https://doi.org/10.1093/nar/gkaf1517
- Cohen-Sharir Y et al. 2021. Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition. *Nature* 590:486-491. PMID 33505028. https://doi.org/10.1038/s41586-020-03114-6
- Conze DB, Crespo-Barreto J, Kruger CL. 2016. Safety assessment of nicotinamide riboside, a form of vitamin B3. *Hum Exp Toxicol* 35:1149-60. PMID 26791540. https://doi.org/10.1177/0960327115626254
- Dixon-Zegeye M et al. 2024. Cancer Precision-Prevention trial of Metformin in adults with Li Fraumeni syndrome (MILI) undergoing yearly MRI surveillance: a randomised controlled trial protocol. *Trials* 25:103. PMID 38308321. https://doi.org/10.1186/s13063-024-07929-w
- Eiden C et al. 2009. Posaconazole-increased vincristine neurotoxicity in a child: a case report. *J Pediatr Hematol Oncol* 31:292-5. PMID 19346885. https://doi.org/10.1097/MPH.0b013e31819b9d01
- Elhajouji A, Tibaldi F, Kirsch-Volders M. 1997. Indication for thresholds of chromosome non-disjunction versus chromosome lagging induced by spindle inhibitors in vitro in human lymphocytes. *Mutagenesis* 12:133-40. PMID 9175637. https://doi.org/10.1093/mutage/12.3.133
- Feng B et al. 2022. A novel CEP57 variant associated with mosaic variegated aneuploidy syndrome in a Chinese female presenting with short stature, microcephaly, brachydactyly, and small teeth. *Mol Genet Genomic Med* 10:e1951. PMID 35434947. https://doi.org/10.1002/mgg3.1951
- Guccini S et al. 2012. Effects of spindle poisons in peripheral human lymphocytes by the in vitro cytokinesis-block micronucleus assay. *Mutagenesis* 27:749-58. PMID 22952150. https://doi.org/10.1093/mutage/ges044
- Hanks S et al. 2004. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. *Nat Genet* 36:1159-61. PMID 15475955. https://doi.org/10.1038/ng1449
- Harwood M et al. 2007. A critical review of the data related to the safety of quercetin and lack of evidence of in vivo toxicity, including lack of genotoxic/carcinogenic properties. *Food Chem Toxicol* 45:2179-205. PMID 17698276. https://doi.org/10.1016/j.fct.2007.05.015
- Hol JA et al. 2021. Wilms tumour surveillance in at-risk children: Literature review and recommendations from the SIOP-Europe Host Genome Working Group and SIOP Renal Tumour Study Group. *Eur J Cancer* 153:51-63. PMID 34134020. https://doi.org/10.1016/j.ejca.2021.05.014
- Hübner CT et al. 2022. Mosaic variegated aneuploidy syndrome and Noonan syndrome in the same family. *Mol Syndromol* 13:402-408. PMID 36588761. https://doi.org/10.1159/000523877
- Izumi H et al. 2009. BubR1 localizes to centrosomes and suppresses centrosome amplification via regulating Plk1 activity in interphase cells. *Oncogene* 28:2806-20. PMID 19503101. https://doi.org/10.1038/onc.2009.141
- Justice JN et al. 2019. Senolytics in idiopathic pulmonary fibrosis: Results from a first-in-human, open-label, pilot study. *EBioMedicine* 40:554-563. PMID 30616998. https://doi.org/10.1016/j.ebiom.2018.12.052
- Kajii T et al. 2001. Cancer-prone syndrome of mosaic variegated aneuploidy and total premature chromatid separation: report of five infants. *Am J Med Genet* 104:57-64. PMID 11746029. https://doi.org/10.1002/ajmg.1580
- Kalish JM et al. 2017. Surveillance Recommendations for Children with Overgrowth Syndromes and Predisposition to Wilms Tumors and Hepatoblastoma. *Clin Cancer Res* 23:e115-e122. PMID 28674120. https://doi.org/10.1158/1078-0432.CCR-17-0710
- Kawakami R et al. 2023. RNA sequencing and target long-read sequencing reveal an intronic transposon insertion causing aberrant splicing. *J Hum Genet* 69:91-99. PMID 38102195. https://doi.org/10.1038/s10038-023-01211-8
- Ketteler M et al. 2023. Modified-release nicotinamide for the treatment of hyperphosphataemia in haemodialysis patients: 52-week efficacy and safety results of the phase 3 randomized controlled NOPHOS trial. *Nephrol Dial Transplant* 38:982-991. PMID 35751625. https://doi.org/10.1093/ndt/gfac206
- Kim J et al. 2019. Patient-Customized Oligonucleotide Therapy for a Rare Genetic Disease. *N Engl J Med* 381:1644-1652. PMID 31597037. https://doi.org/10.1056/NEJMoa1813279
- Kurokawa Y et al. 2025. Pooled safety analysis of pimitespib for the treatment of patients with advanced gastrointestinal stromal tumors. *Int J Clin Oncol* 31:118-127. PMID 41264161. https://doi.org/10.1007/s10147-025-02917-9
- Laberko A et al. 2019. Hematopoietic stem cell transplantation in a patient with type 1 mosaic variegated aneuploidy syndrome. *Orphanet J Rare Dis* 14:97. PMID 31053147. https://doi.org/10.1186/s13023-019-1073-x
- Lamb YN. 2024. Givinostat: First Approval. *Drugs* 84:849-856. PMID 38967716. https://doi.org/10.1007/s40265-024-02052-1
- Lim KH et al. 2020. Antisense oligonucleotide modulation of non-productive alternative splicing upregulates gene expression. *Nat Commun* 11:3501. PMID 32647108. https://doi.org/10.1038/s41467-020-17093-9
- Mahadevan D et al. 2025. First in human phase 1 study of DT2216, a selective BCL-xL degrader, in patients with relapsed/refractory solid malignancies. *J Hematol Oncol* 18:98. PMID 41225624. https://doi.org/10.1186/s13045-025-01753-8
- Malumbres M, Villarroya-Beltri C. 2024. Mosaic variegated aneuploidy in development, ageing and cancer. *Nat Rev Genet* 25:864-878. PMID 39169218. https://doi.org/10.1038/s41576-024-00762-6
- Matsuura S et al. 2000. Chromosomal instability syndrome of total premature chromatid separation with mosaic variegated aneuploidy is defective in mitotic-spindle checkpoint. *Am J Hum Genet* 67:483-6. PMID 10877982. https://doi.org/10.1086/303022
- McDonald CM et al. 2025. Long-term evaluation of givinostat in Duchenne muscular dystrophy, and natural history comparisons. *Ann Clin Transl Neurol* 12:2335-2348. PMID 40830818. https://doi.org/10.1002/acn3.70165
- McWilliams GD et al. 2016. Occurrence of DNET and other brain tumors in Noonan syndrome warrants caution with growth hormone therapy. *Am J Med Genet A* 170A:195-201. PMID 26377682. https://doi.org/10.1002/ajmg.a.37379
- Michaelis M et al. 2004. Valproic acid inhibits angiogenesis in vitro and in vivo. *Mol Pharmacol* 65:520-7. PMID 14978230. https://doi.org/10.1124/mol.65.3.520
- Miyamoto T et al. 2011. Insufficiency of BUBR1, a mitotic spindle checkpoint regulator, causes impaired ciliogenesis in vertebrates. *Hum Mol Genet* 20:2058-70. PMID 21389084. https://doi.org/10.1093/hmg/ddr090
- Miyamoto T et al. 2015. The Microtubule-Depolymerizing Activity of a Mitotic Kinesin Protein KIF2A Drives Primary Cilia Disassembly Coupled with Cell Proliferation. *Cell Rep* 10:664-673. PMID 25660017. https://doi.org/10.1016/j.celrep.2015.01.003
- Muraoka H et al. 2025. Pimitespib, an HSP90 inhibitor, enhances the efficacy of PARP inhibitors in PARP inhibitor-insensitive breast cancer cells. *Cancer Sci* 116:1745-1757. PMID 40167031. https://doi.org/10.1111/cas.70058
- Nakano Y et al. 2024. Update on Recommendations for Cancer Screening and Surveillance in Children with Genomic Instability Disorders. *Clin Cancer Res* 30:5009-5020. PMID 39264246. https://doi.org/10.1158/1078-0432.CCR-24-1098
- Nambiar A et al. 2023. Senolytics dasatinib and quercetin in idiopathic pulmonary fibrosis: results of a phase I, single-blind, single-center, randomized, placebo-controlled pilot trial on feasibility and tolerability. *EBioMedicine* 90:104481. PMID 36857968. https://doi.org/10.1016/j.ebiom.2023.104481
- North BJ et al. 2014. SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. *EMBO J* 33:1438-53. PMID 24825348. https://doi.org/10.15252/embj.201386907
- Ochiai H et al. 2013. TALEN-mediated single-base-pair editing identification of an intergenic mutation upstream of BUB1B as causative of PCS (MVA) syndrome. *Proc Natl Acad Sci U S A* 111:1461-6. PMID 24344301. https://doi.org/10.1073/pnas.1317008111
- Ochiai K et al. 2019. Long-term remission of bilateral Wilms tumors that developed from premature separation of chromatids/mosaic variegated aneuploidy syndrome due to bilateral nephrectomy and peritoneal dialysis. *Pediatr Blood Cancer* 66:e27804. PMID 31081598. https://doi.org/10.1002/pbc.27804
- Olaharski AJ et al. 2006. The histone deacetylase inhibitor trichostatin a has genotoxic effects in human lymphoblasts in vitro. *Toxicol Sci* 93:341-7. PMID 16857700. https://doi.org/10.1093/toxsci/kfl068
- Oromendia AB et al. 2012. Aneuploidy causes proteotoxic stress in yeast. *Genes Dev* 26:2696-708. PMID 23222101. https://doi.org/10.1101/gad.207407.112
- Oztas S et al. 2007. Genotoxic effect of albendazole in pediatric patients with hepatic hydatid disease. *Int J Infect Dis* 11:446-9. PMID 17331781. https://doi.org/10.1016/j.ijid.2007.01.002
- Panvino F et al. 2026. Pharmacological Treatment of Cerebellar Ataxia in Pediatric Ataxia-Telangiectasia: A Systematic Review. *Eur J Neurol* 33:e70683. PMID 42387837. https://doi.org/10.1111/ene.70683
- Park I et al. 2013. Loss of BubR1 acetylation causes defects in spindle assembly checkpoint signaling and promotes tumor formation. *J Cell Biol* 202:295-309. PMID 23878276. https://doi.org/10.1083/jcb.201210099
- Park I et al. 2017. HDAC2/3 binding and deacetylation of BubR1 initiates spindle assembly checkpoint silencing. *FEBS J* 284:4035-4050. PMID 28985013. https://doi.org/10.1111/febs.14286
- Pavone P et al. 2022. Pathogenic correlation between mosaic variegated aneuploidy 1 (MVA1) and a novel BUB1B variant: a reappraisal of a severe syndrome. *Neurol Sci* 43:6529-6538. PMID 35804254. https://doi.org/10.1007/s10072-022-06247-w
- Phillips AF et al. 2025. Targeting chromosomally unstable tumors with a selective KIF18A inhibitor. *Nat Commun* 16:307. PMID 39747049. https://doi.org/10.1038/s41467-024-55300-z
- Qiu D et al. 2018. Sirt2-BubR1 acetylation pathway mediates the effects of advanced maternal age on oocyte quality. *Aging Cell* 17. PMID 29067790. https://doi.org/10.1111/acel.12698
- Ramírez T, Eastmond DA, Herrera LA. 2006. Non-disjunction events induced by albendazole in human cells. *Mutat Res* 626:191-5. PMID 17055335. https://doi.org/10.1016/j.mrgentox.2006.09.004
- Rio Frio T et al. 2010. Homozygous BUB1B mutation and susceptibility to gastrointestinal neoplasia. *N Engl J Med* 363:2628-37. PMID 21190457. https://doi.org/10.1056/NEJMoa1006565
- Sato Y, Hayashi MT. 2024. Micronucleus is not a potent inducer of the cGAS/STING pathway. *Life Sci Alliance* 7:e202302424. PMID 38307626. https://doi.org/10.26508/lsa.202302424
- Sieben CJ et al. 2020. BubR1 allelic effects drive phenotypic heterogeneity in mosaic-variegated aneuploidy progeria syndrome. *J Clin Invest* 130:171-188. PMID 31738183. https://doi.org/10.1172/JCI126863
- Snape K et al. 2011. Mutations in CEP57 cause mosaic variegated aneuploidy syndrome. *Nat Genet* 43:527-9. PMID 21552266. https://doi.org/10.1038/ng.822
- Suematsu T et al. 2014. Deacetylation of the mitotic checkpoint protein BubR1 at lysine 250 by SIRT2 and subsequent effects on BubR1 degradation during the prometaphase/anaphase transition. *Biochem Biophys Res Commun* 453:588-94. PMID 25285631. https://doi.org/10.1016/j.bbrc.2014.09.128
- Suijkerbuijk SJ et al. 2010. Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. *Cancer Res* 70:4891-900. PMID 20516114. https://doi.org/10.1158/0008-5472.CAN-09-4319
- Suijkerbuijk SJ et al. 2012a. The vertebrate mitotic checkpoint protein BUBR1 is an unusual pseudokinase. *Dev Cell* 22:1321-9. PMID 22698286. https://doi.org/10.1016/j.devcel.2012.03.009
- Suijkerbuijk SJ et al. 2012b. Integration of kinase and phosphatase activities by BUBR1 ensures formation of stable kinetochore-microtubule attachments. *Dev Cell* 23:745-55. PMID 23079597. https://doi.org/10.1016/j.devcel.2012.09.005
- Tang YC et al. 2011. Identification of aneuploidy-selective antiproliferation compounds. *Cell* 144:499-512. PMID 21315436. https://doi.org/10.1016/j.cell.2011.01.017
- Villarroya-Beltri C et al. 2022. Biallelic germline mutations in MAD1L1 induce a syndrome of aneuploidy with high tumor susceptibility. *Sci Adv* 8:eabq5914. PMID 36322655. https://doi.org/10.1126/sciadv.abq5914
- Viudes CP et al. 2026. A novel CEP57 gene mutation in mosaic variegated aneuploidy syndrome 2: case report. *J Pediatr Endocrinol Metab* 39:381-387. PMID 41700350. https://doi.org/10.1515/jpem-2025-0643
- Vohra S et al. 2015. CONSORT extension for reporting N-of-1 trials (CENT) 2015 Statement. *BMJ* 350:h1738. PMID 25976398. https://doi.org/10.1136/bmj.h1738
- Watanabe Y et al. 2014. Dendrite development regulated by the schizophrenia-associated gene FEZ1 involves the ubiquitin proteasome system. *Cell Rep* 7:552-564. PMID 24726361. https://doi.org/10.1016/j.celrep.2014.03.022
- Woodward HL et al. 2018. Introduction of a Methyl Group Curbs Metabolism of Pyrido[3,4- d]pyrimidine Monopolar Spindle 1 (MPS1) Inhibitors and Enables the Discovery of the Phase 1 Clinical Candidate N-(2-Ethoxy-4-(4-methyl-4 H-1,2,4-triazol-3-yl)phenyl)-6-methyl- N-neopentylpyrido[3,4- d]pyrimidine-2,8-diamine (BOS172722). *J Med Chem* 61:8226-8240. PMID 30199249. https://doi.org/10.1021/acs.jmedchem.8b00690
- Wu LE et al. 2023. SIRT2 transgenic over-expression does not impact lifespan in mice. *Aging Cell* 22:e14027. PMID 38009412. https://doi.org/10.1111/acel.14027
- Yost S et al. 2017. Biallelic TRIP13 mutations predispose to Wilms tumor and chromosome missegregation. *Nat Genet* 49:1148-1151. PMID 28553959. https://doi.org/10.1038/ng.3883
- Zhu Y et al. 2015. The Achilles' heel of senescent cells: from transcriptome to senolytic drugs. *Aging Cell* 14:644-58. PMID 25754370. https://doi.org/10.1111/acel.12344

Database and registry sources (queried 2026-09-27): ClinicalTrials.gov (NCT01981525, NCT06799065, NCT05902988, NCT06084416, NCT04293094, NCT07226427, NCT06772415); ISRCTN (ISRCTN16699730; searches for mosaic variegated aneuploidy, premature chromatid separation, BUB1B); ChEMBL target CHEMBL4295998 (assays incl. CHEMBL6125423, CHEMBL6125447, CHEMBL4807760); Open Targets Platform target ENSG00000156970.

Queried 2026-09-29: EU CTIS public search (2025-520714-63-00; searches for mosaic variegated aneuploidy, premature chromatid separation, BUB1B, Li-Fraumeni); ANZCTR and the WHO ICTRP Search Portal (the same MVA searches, plus Li-Fraumeni controls); Ensembl REST API, gene ENSG00000156970 with all transcript models (GRCh38; MANE Select ENST00000287598.11; NMD transcripts include ENST00000558715, ENST00000559733 and ENST00001082076–ENST00001082085).

US prescribing information (openFDA drug label API, read 2026-09-29; set IDs resolve on DailyMed): givinostat, DUVYZAT (1a28a3db-7880-4ac2-a2f2-2b4e362aafb6); valproic acid (1edd220d-d151-4177-9a29-9fee36489fb7); metformin hydrochloride extended-release (05999192-ebc6-4198-bd1e-f46abbfb4f8a); dasatinib (58ce42c6-10bd-4bb4-bf89-eae8ed82c339); nusinersen, SPINRAZA (dd70cd5f-b0fc-4ba4-a5ea-89a34778bd94); mebendazole (a720c424-517a-461b-8e13-6c37caf8aa3f); albendazole (8de337eb-92e1-15b5-e053-2a95a90ac4f1).
