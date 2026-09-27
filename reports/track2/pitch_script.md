# Track 2 pitch: script and slide outline (3 minutes)

**For:** the 3-minute video (package U6, record by 2026-10-21). **Pace:** about 140 words a minute. The spoken text below is about 360 words, roughly 2:35 when read aloud, which leaves about 25 seconds for slide transitions and pauses. Timings are cumulative.

**Rules for the video.** Hypotheses only: no doses, no schedules, no advice to the family. Say "the child", never a name, sex, age, date or place beyond "UK". Show no genotypes or variant coordinates. The on-screen footer on every slide reads: *Research hypotheses for expert follow-up. Not medical advice.*

**Before recording:** if report v2 (package D08) has re-ranked the hypotheses using the child's alleles, update slide 4 and the line marked [v2].

---

## Slide 1 · 0:00–0:25 · The veto rule

**On screen:** one sentence, large: *"Any drug that shortens mitosis is out — whatever else it does."* Underneath, smaller: *the five-assay decision rule, §7.*

**Say:**
> Every drug in this pitch has to survive one rule first. If it shortens mitosis in the patient's own cells, it is out, however good its other results look. That rule comes straight from the biology of this disease, and it is why several obvious drug classes are excluded before we start.

## Slide 2 · 0:25–0:55 · A dosage disease

**On screen:** the acetylation–degradation cycle as a loop (PCAF acetylates K250 → protected; HDAC2/3 remove the mark → APC/C destroys BUBR1). A red tag: *acetylation-deficient BubR1 → shorter mitosis → MVA-like mice* (Park 2013).

**Say:**
> MVA type 1 is a dosage disease. One copy of *BUB1B* makes no protein, the other makes too little, and adding BUBR1 back rescues patient cells. So the goal is more BUBR1, not a stronger kinase: BUBR1 is a pseudokinase. Its level is controlled by an acetylation cycle, and mice that cannot acetylate it get shorter mitoses and look like MVA. Shorter mitosis is the fingerprint of the disease. That is the veto.

## Slide 3 · 0:55–1:20 · Nothing exists, and the databases point the wrong way

**On screen:** three zeros: *0 trials (ClinicalTrials.gov, ISRCTN) · 0 disease-modifying reports (PubMed) · best BUB1B ligand ≈ 10 µM (ChEMBL)*. Two warning icons: *"Druggable family" = kinase heuristic · "Degrader-tractable" = the harmful direction.*

**Say:**
> We checked the registries, PubMed, ChEMBL and Open Targets. There has never been a trial for MVA. The only positive tractability signals are traps: one comes from the kinase family this protein no longer behaves like, and the other says it can be degraded, which is the opposite of what a dosage disease needs.

## Slide 4 · 1:20–2:00 · Five hypotheses, each with its weakness

**On screen:** a five-row table: hypothesis · acts on · main weakness.

| | Hypothesis | Main weakness |
|---|---|---|
| H1 | NAD+ precursors via SIRT2 | SIRT2 also removes the protective mark |
| H2 | HDAC2/3 inhibition | the class can cause aneuploidy and low platelets |
| H3 | Metformin inside surveillance (MILI design) | may select against the child's own mosaic tissue |
| H4 | Genotype-matched oligonucleotide | depends on the child's alleles |
| H5 | Senolytics | senescence may be a brake on tumours |

**Say:**
> We rank five hypotheses by how cheaply the first experiment could disprove them. NAD+ precursors and HDAC2/3 inhibition both aim to protect BUBR1 itself. Metformin could sit inside the kidney surveillance already recommended in MVA, following the UK MILI trial design. A genotype-matched oligonucleotide depends on the child's own alleles. Senolytics come last. Every row carries its main weakness, because a hypothesis without one isn't ready for an expert. [v2: once the alleles are classified, say here whether H4 moves up.]

## Slide 5 · 2:00–2:35 · Two filters, and a pipeline that applies them

**On screen:** left, the six-point chemotherapy-compatibility screen as six icons (checkpoint · CYP3A4 · kidney · marrow · genotoxicity · quick to stop). Right, a bar from the R06 demo: *28 of 28 checkpoint-kinase controls excluded automatically · 11 of 12 HDAC2/3 candidates excluded*.

**Say:**
> A child with MVA may need cancer treatment at short notice, so every idea must also pass a six-point compatibility screen: no further weakening of the checkpoint, no CYP3A4 clash with vincristine, and no harm to kidney or bone marrow. We built an open pipeline that runs these checks over public data. On its first run it excluded every checkpoint-kinase inhibitor we fed it as a control, with no hand curation.

## Slide 6 · 2:35–2:55 · How it scales and where it goes next

**On screen:** three genes, one config file: *BUB1B · TRIP13 · CEP57*. A path: *patient cells → five-assay rule → N-of-1 design (CENT) → registry*.

**Say:**
> Nothing in the rule depends on the gene. The same pipeline runs for *TRIP13* and *CEP57* from one line of configuration. The path is patient-derived cells first, then the five-assay rule, and only then a formal N-of-1 design.

## Slide 7 · 2:55–3:00 · Close

**On screen:** *Veto first. Hypotheses, not prescriptions.* Repository link.

**Say:**
> Veto first, then hypotheses for the experts. Thank you.

---

## Checks before upload
- [ ] Runs at 3:00 or less when read aloud at a natural pace.
- [ ] No name, sex, age, date or place beyond "UK"; no genotype or coordinate on any slide.
- [ ] No dose, schedule or "should take" wording anywhere, spoken or on screen.
- [ ] Slide 4 matches the latest Track 2 report version.
- [ ] Numbers on slide 5 match `reports/track2/drug_triage/drug_triage_demo.md` from the latest run.
- [ ] Uploading to YouTube or Vimeo, and the visibility setting, need Santanu's live yes.
