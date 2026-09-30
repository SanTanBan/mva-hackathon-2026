# D08 worksheet: how the child's *BUB1B* alleles will re-rank the Track 2 hypotheses

**Rare Disease, Real Kid: The MVA Hackathon 2026 · Track 2 companion to the report · rules fixed 2026-09-29, before dataset access · public version**

> **Why this exists.** Package D08 has to classify each of the child's two *BUB1B* alleles by mechanism, assess whether an antisense oligonucleotide (ASO) could target it, and revisit the Track 2 ranking. We wrote the rules down, and coded and tested them, **before seeing any patient data**. The re-ranking therefore cannot be tuned to the answer. The code is `pipeline/allele_mechanism.py`, with 18 tests in `pipeline/tests/test_allele_mechanism.py`.
>
> **Privacy.** This page contains no patient data. The worked examples are published ClinVar or literature alleles. The public version of the D08 result will report classes and tiers only. It will not give the child's variant names or coordinates, and it will not say whether any worked example matches the child.

All variant names are on NM_001211.6 (MANE Select). Exon numbers and domain positions follow the locus dossier (`reports/track1/locus_dossier.md`, §2). Sequence facts (exon frames, junction codons, KEN-box positions) were computed from the RefSeq record NM_001211.6 (NCBI, retrieved 2026-09-29).

---

## Step 1. Mechanism class

Suijkerbuijk 2010 (PMID 20516114) replaced endogenous BUBR1 in human cells with each MVA mutant and studied patient cell lines. Every MVA mutation fell into one of two classes:

- **Abundance class**: the allele lowers BUBR1 protein. Truncating alleles make no transcript. Substitutions in or near the kinase-like domain make a misfolding-prone protein that turns over fast and needs HSP90 to fold.
- **Functional class**: the protein is stable but carries a specific defect.

Two further results from that paper shape every rule below:

- When the abundance-class substitutions I909T and L1012P were overexpressed to wild-type levels, they fully restored the checkpoint. **For this class, raising the level is enough.**
- In HeLa cells, cutting BUBR1 to about 6% of normal caused missegregation in most cells, while about 13% had virtually no effect. **Chromosome segregation may have a threshold, so modest increases could matter.**

**The substitutions that paper tested:**

| Protein change | Exon | Class | What the experiment showed |
|---|---|---|---|
| p.Tyr155Cys | 5 (TPR domain) | functional | Stable. It restored chromosome alignment but not checkpoint activity (separation of function). |
| p.Arg550Gln | 14 | no defect detected | Stable; rescued both functions |
| p.Arg727Cys | 17 | abundance | Unstable, with faster turnover in patient cells |
| p.Arg814His | 19 | abundance | Unstable in cells |
| p.Leu844Phe | 19 | abundance | Poorly expressed; failed to rescue |
| p.Ile909Thr | 21 | abundance | Unstable and HSP90-dependent; fully functional when overexpressed |
| p.Gln921His | 21 | no defect detected | Indistinguishable from wild type; a passenger on the L844F allele |
| p.Leu1012Pro | 23 | abundance | Unstable and HSP90-dependent; fully functional when overexpressed |

**Rules for alleles not in that table:**

- **Truncating allele, with NMD predicted:** abundance class.
- **Truncating allele that escapes NMD:** unclassified. The BUBR1 immunoblot (report §7, assay 1) decides.
- **Substitution at residue 727 or later:** provisionally abundance class. Every tested substitution in that region was unstable.
- **Substitution at residue 550 or earlier:** provisionally functional class or no defect.
- **Substitution at residues 551–726:** unclassified, because nothing there has been tested.
- **Intronic, splice-region or upstream allele:** provisionally abundance class, pending RNA.

## Step 2. Nonsense-mediated decay (NMD)

These rules are from Lindeboom 2016 (PMID 27618451), based on matched exome and transcriptome data from 9,769 tumours.

| Premature stop codon position | Prediction |
|---|---|
| Last exon (exon 23, from c.2958) | escapes NMD |
| Last 50 nt of exon 22 (c.2908–2957) | escapes NMD |
| First 100 coding nucleotides | rarely targeted |
| Coding nucleotides 101–200 | only partly targeted |
| Anywhere else | NMD predicted |

A frameshift is scored at its new stop codon, not at the variant itself. When the annotation does not give the stop position, the code uses the variant position and flags the call as approximate.

## Step 3. Antisense lanes

**Tiers.** We use the five tiers of the N=1 Collaborative's VARIANT guidelines (N1C VARIANT, "variant assessments toward eligibility for antisense oligonucleotide treatment"; Cheerie 2025, PMID 40139194): eligible, likely eligible, unlikely eligible, not eligible, and unable to assess. We apply the splice-switching and exon-skipping criteria of the Dutch Center for RNA Therapeutics (Zardetto 2024, PMID 40225926).

- **"Eligible" needs functional evidence** that an ASO works for that variant or exon.
- **None exists for *BUB1B*.** A PubMed search on 2026-09-29 ("BUB1B AND (antisense OR oligonucleotide)", 16 records) found none. The best tier available here is therefore "likely eligible".
- **The tiers are about the variant only.** The guideline authors stress that an eligible variant does not make an eligible person (Cheerie 2025).

### 3a. Splice correction: blocking a cryptic splice site

| Situation | Tier |
|---|---|
| Canonical ±1/2 site destroyed | not eligible: an ASO cannot rebuild a splice site |
| RNA shows the variant only weakens the canonical site (exon skipping or intron retention) | not eligible |
| RNA-confirmed cryptic site under 5 nt from the canonical site | not eligible (hard minimum) |
| RNA-confirmed cryptic site 5–14 nt away | unlikely eligible (under the recommended 15 nt) |
| Acceptor side, 18–40 nt upstream (branch-point window) | unlikely eligible |
| RNA-confirmed cryptic site at least 15 nt away | likely eligible |
| RNA-confirmed cryptic site deep intronic (more than 100 nt) | likely eligible (the ideal case) |
| Any intronic or splice-region variant without RNA data | **unable to assess**. The code also reports the tier the variant would get if RNA confirmed a cryptic site. |

A synonymous exonic variant with an RNA-confirmed cryptic site at least 15 nt from the exon end is likely eligible. A missense variant with one is unlikely eligible, because the amino-acid change may be pathogenic on its own.

**The dataset has no RNA** (dossier §4). Any splice-type allele will therefore stay "unable to assess" until RNA from patient cells is available. This makes RNA sequencing the first experiment for such an allele.

#### What the BUB1B splice map shows (added 2026-10-01)

We scored all 180,165 possible single-nucleotide changes across the *BUB1B* MANE Select transcript with Pangolin (Zeng & Li 2022, PMID 35449021; GPL-3.0). The run used public sequence only, on Kaggle (`kaggle/bub1b_splice_map/`, Pangolin commit 5cf94b8).

**Checks.**
- Our batched scorer matched Pangolin's own scoring function to within 5×10⁻⁵ (2×10⁻⁴ in half precision).
- All seven published positive controls behaved as expected:
  - The canonical splice variants c.1402-1G>T, c.2386-2A>G and c.581+1G>T each lose about 0.85.
  - c.1401G>T, on the last base of exon 10, loses 0.84 at the donor. This supports the dossier's inference that this ClinVar "missense" acts through splicing.
  - c.2386-11A>G gains 0.76 at a new acceptor 10 nt upstream while the canonical acceptor loses 0.71, as Rio Frio 2010 reported.

Each variant's largest gain and loss within ±50 nt, after masking MANE Select splice sites:

| Where the change is | Possible changes | Strong cryptic gain (≥ 0.5) | Strong loss (≤ −0.5) |
|---|---|---|---|
| Exonic | 11,007 | 51 | 68 |
| Canonical site (±1–2) | 264 | 86 | 256 |
| Intronic, 3–4 nt from the exon | 264 | 12 | 55 |
| Intronic, 5–14 nt | 1,320 | 30 | 102 |
| Branch-point window (acceptor −18 to −40) | 1,518 | 1 | 7 |
| Intronic, 15–100 nt | 9,516 | 11 | 7 |
| Deep intronic (more than 100 nt) | 156,276 | 64 | 0 |
| **Total** | **180,165** | **255** | **495** |

**What it means for the splice-correction lane:**
- **The target class is rare.** Only 75 of the 169,158 possible intronic changes (0.04%) create a strong cryptic site at least 15 nt from a canonical site and outside the branch-point window. 64 of these are deep intronic, at 58 positions across 13 introns; introns 4 and 8 have the most, with 10 each.
- **The two published cryptic-acceptor alleles fall outside it.** c.2386-11A>G and c.1402-5A>G both create strong new acceptors, but only 10 and 4 nt upstream of the canonical acceptor. That is inside the zone where splice correction is "unlikely eligible". For the *BUB1B* second alleles published so far, splice correction is therefore unlikely to apply. Up-regulation remains the more plausible antisense route.
- **At the looser threshold of 0.2,** 398 intronic changes at least 15 nt from a canonical site show some gain.

The map covers single-nucleotide changes only. Indels and transposon insertions need their own Pangolin run, and a prediction is not RNA evidence. D04 looks up the child's single-nucleotide variants locally with `python -m pipeline.splice_map`, so they never leave this machine.

### 3b. Exon skipping: removing the exon that carries a truncating variant

An exon can be skipped only if it is in-frame and is neither the first nor the last coding exon. It must also not encode an essential region. Exon 1 and exon 23 are excluded. Of the 21 internal exons, only four are in-frame, and none of the four is a good candidate:

| Exon | Residues removed | Verdict | Why |
|---|---|---|---|
| 2 | 12–60 (4.6%) | not eligible | Contains the first KEN box (residues 26–28), which BUBR1 needs to form the mitotic checkpoint complex (Di Fiore 2015, PMID 25669885). The exon 2 acceptor variant c.36-1G>A is likely pathogenic in ClinVar. |
| 3 | 60–80 (1.9%) | unlikely eligible | Removes the start of the TPR domain (62–226), which checkpoint signalling needs (Suijkerbuijk 2010). The exon 3 donor variant c.239+2T>C is likely pathogenic in ClinVar. |
| 17 | 715–762 (4.5%) | unlikely eligible | For skipping, it has no annotated domain or KEN box and no pathogenic splice-site variant at its boundaries. Residues 1–715 alone supported checkpoint and attachment in cell lines (the 731X construct; Suijkerbuijk 2010). **Against skipping**, the MVA substitution R727C in this exon destabilises BUBR1, so the segment probably helps the kinase-like domain fold. Zardetto 2024 treats a destabilising or inactivating missense in an exon as evidence that skipping it would fail too. |
| 19 | 796–845 (4.8%) | not eligible | Lies inside the kinase-like domain, whose integrity keeps BUBR1 stable (Suijkerbuijk 2010). R814H and L844F in this exon lower BUBR1 levels. The exon 19 acceptor variant c.2386-2A>G is likely pathogenic in ClinVar. |

Skipping any of these four creates no new stop codon at the junction (checked on the reference sequence). Adjacent out-of-frame pairs that together keep the frame are exons 4+5, 12+13, 14+15, 20+21 and 21+22. Each pair would remove 37 to 127 residues from a functional region and would need two ASOs, so they are not pursued.

**Exon 17 is still worth one cheap experiment.** It is the only in-frame exon clear of the KEN boxes and of the TPR and kinase-like domains. Four of the 86 named ClinVar pathogenic or likely pathogenic *BUB1B* variants lie in it (dossier appendix). If the child's truncating allele is one of them or another exon 17 variant, the first experiment is this: express BUBR1 lacking exon 17 (or treat patient cells with an exon-17-skipping ASO), then run assay 1 (protein level) and assay 3 (checkpoint).

### 3c. Knockdown

Not eligible. MVA1 is caused by too little BUBR1, and the N1C guidelines list dosage-sensitive genes as not eligible for knockdown (Cheerie 2025).

### 3d. Up-regulation (reported, but not given a tier)

The N1C guidelines give considerations for up-regulation but define no tiers for it. The idea here is to block non-productive splicing so that more productive *BUB1B* mRNA is made (Lim 2020; candidate events are in report §8.1). Such an ASO acts on both alleles, but only an allele with a usable product benefits:

| Allele | Up-regulation |
|---|---|
| Abundance-class substitution | **candidate**: more of a protein that works when abundant (I909T and L1012P) |
| Regulatory or leaky-splice allele that still makes some normal transcript | candidate |
| Truncating allele | not for this allele: its transcripts carry the stop codon |
| Functional-class substitution | unlikely to help: more of a defective protein |

## Step 4. Pre-registered changes to the Track 2 ranking

The current order (report §4) is H1 NAD+ precursors, H2 HDAC2/3 inhibition, H3 metformin inside surveillance, H4 genotype-matched oligonucleotide, H5 senolytics.

| If D04 finds… | Then |
|---|---|
| An allele with an RNA-confirmed cryptic splice site at least 15 nt from the canonical site | **H4 moves to rank 1**: splice correction of that allele. MVA1 is recessive (biallelic in every PanelApp panel, dossier §6), so one corrected allele is expected to be enough. |
| A truncating allele in exon 17 | H4 stays at rank 4. The first experiment is the exon-17-skip stability test above. H4 rises only if the shortened protein is stable (assay 1) and supports the checkpoint (assay 3). |
| An intronic, splice-region or transposon allele without RNA data | H4 stays at rank 4. RNA sequencing of patient cells, with and without NMD inhibition, becomes the first experiment, because it decides both the splice-correction and the up-regulation lanes. |
| The non-null allele is an abundance-class substitution | H1 and H2 keep their place with stronger support, because raising the level of these mutants restores function. HSP90 inhibitors get a second reason for exclusion, because HSP90 inhibition sharply lowers these mutants (Suijkerbuijk 2010). The up-regulation lane becomes a candidate. |
| The non-null allele is a stable, functional-class substitution (like Y155C) | Goal A loses its direct rationale for that allele. The order becomes **H3, H1, H2, H5**, and H4 applies only if a splice or skip lane does. |
| A regulatory (upstream or promoter) second allele | Up-regulation is a candidate lane. H4 stays at rank 4, and H1 and H2 are unchanged. |
| Two apparently null alleles | Re-check NMD escape (Step 2) and leaky splicing before interpreting. In the biallelic families studied by Suijkerbuijk 2010, a truncation always paired with a substitution. |
| No second allele found | H4 stays conditional at rank 4, and the gap is reported as package D04 requires. |

## Worked examples on published alleles

Tiers from `python -m pipeline.allele_mechanism`, run on alleles already in the literature or in ClinVar:

| Allele | Mechanism | Splice correction | Exon skipping | Up-regulation |
|---|---|---|---|---|
| c.2210T>G, p.(Leu737Ter), exon 17 | abundance (NMD) | not applicable | unlikely eligible (R727C) | not for this allele |
| c.1833del, p.(Phe611fs), exon 15 | abundance (NMD) | not applicable | not eligible (out of frame) | not for this allele |
| c.3035T>C, p.(Leu1012Pro) | abundance (tested) | not applicable | not eligible (last exon) | candidate |
| c.2386-11A>G, a new splice site favoured over the authentic one (Rio Frio 2010); scored with that published RNA evidence | abundance (provisional) | unlikely eligible (11 nt from the acceptor) | not applicable | candidate if some normal transcript remains |
| c.1288+5G>A (Laberko 2019), no RNA data | abundance (provisional) | unable to assess | not applicable | candidate if some normal transcript remains |

**The pattern matters for Track 2.** For the published MVA1 alleles, the realistic antisense routes are:

- splice correction, when a second allele creates a cryptic site far enough from the canonical one;
- up-regulation, when a second allele is an unstable substitution.

Classic exon skipping is not a realistic route for this gene.

## Limits

- The tiers describe molecular eligibility only. Delivery is the larger barrier. Local delivery works best for disorders mainly affecting the brain, spinal cord or eye, and bespoke ASO development takes about two years (Zardetto 2024). MVA1 affects many dividing tissues (report §8.1).
- Mechanism classes for untested substitutions are provisional until assay 1 is run in patient-derived cells.
- The splice rules use distance cut-offs. Real splice outcomes need RNA, which this dataset does not include.

## References

- Cheerie D et al. 2025. Consensus guidelines for assessing eligibility of pathogenic DNA variants for antisense oligonucleotide treatments. *Am J Hum Genet* 112:975–983. PMID 40139194. https://doi.org/10.1016/j.ajhg.2025.02.017
- Di Fiore B et al. 2015. The ABBA motif binds APC/C activators and is shared by APC/C substrates and regulators. *Dev Cell* 32:358–372. PMID 25669885. https://doi.org/10.1016/j.devcel.2015.01.003
- Laberko A et al. 2019. PMID 31053147 (full entry in the Track 2 report).
- Lim KH et al. 2020. PMID 32647108 (full entry in the Track 2 report).
- Lindeboom RGH et al. 2016. The rules and impact of nonsense-mediated mRNA decay in human cancers. *Nat Genet* 48:1112–1118. PMID 27618451. https://doi.org/10.1038/ng.3664
- Rio Frio T et al. 2010. PMID 21190457 (full entry in the Track 2 report).
- Suijkerbuijk SJE et al. 2010. Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. *Cancer Res* 70:4891–4900. PMID 20516114. https://doi.org/10.1158/0008-5472.CAN-09-4319
- Zardetto B et al. 2024. Practical recommendations for the selection of patients for individualized splice-switching ASO-based treatments. *Hum Mutat* 2024:9920230. PMID 40225926. https://doi.org/10.1155/2024/9920230
- Zeng T, Li YI. 2022. Predicting RNA splicing from DNA sequence using Pangolin. *Genome Biol* 23:103. PMID 35449021. https://doi.org/10.1186/s13059-022-02664-4
