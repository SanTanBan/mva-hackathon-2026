# BUB1B locus dossier (public data only)

**Rare Disease, Real Kid: The MVA Hackathon 2026 · Track 1 reference document · compiled 2026-09-27**

> **Scope.** Everything here comes from public databases and published literature. The document contains no patient data and no coordinates of the child's variants. Published variants are cited by ClinVar accession and HGVS name. The only genomic coordinates given are for the public reference gene model (GRCh38), which our analysis code uses to define the regions it inspects.
>
> **Why it exists.** In MVA1 a "solved" case needs **both** *BUB1B* alleles accounted for, and the second allele is often non-coding. This dossier fixes, before we look at the data:
> - where to look;
> - what known pathogenic variation looks like;
> - which traps to avoid;
> - which other genes must be shown negative.
>
> Packages D04 (both alleles, phase, structural variants and mobile-element insertions, upstream region, splicing) and D05 (differential genes, ACMG) use it as their checklist.

Sources and query dates: Ensembl REST (GRCh38), NCBI Datasets and ClinVar E-utilities, UniProt, the gnomAD API, and the Genomics England PanelApp API, all queried 2026-09-27. Literature is cited by PMID; every PMID was checked against its PubMed record.

---

## 1. Gene, transcript and protein

| Item | Value |
|---|---|
| Gene | *BUB1B* (HGNC:1149), NCBI Gene 701, Ensembl ENSG00000156970 |
| Location (GRCh38) | chromosome 15, 40,160,984–40,221,137, **plus strand** (15q15.1). The dataset VCF names this contig `15` with no `chr` prefix; Track 1 CSVs need `chr15`. |
| **MANE Select transcript** | **NM_001211.6 = ENST00000287598.11**: 23 exons, 3,669 nt |
| Protein | NP_001202.5 / UniProt O60566 (BUBR1), 1,050 aa |
| Coding sequence | c.1–3153, including the stop codon; exon 1 carries only 35 coding nucleotides |
| gnomAD constraint | pLI ≈ 0; observed/expected loss-of-function (LoF) 0.62 (90% CI 0.52–0.75, so LOEUF 0.75), from 86 observed vs 137.7 expected LoF variants; missense Z 1.30 (o/e 0.91); synonymous Z −0.06 |

**Reading the constraint.** *BUB1B* tolerates heterozygous loss of function, which is what a recessive disease predicts. A heterozygous LoF variant in the child is therefore expected in carriers and is not evidence on its own. It also means gnomAD's LoF observations cannot argue against pathogenicity in a biallelic setting.

**Always report on NM_001211.6.** A readthrough gene, *BUB1B-PAK6* (ENSG00000259288, HGNC:52276; GRCh38 15:40,218,500–40,264,890, plus strand), overlaps the 3′ end of *BUB1B*. Annotation tools can attach consequences to readthrough transcripts. Pin every consequence call to MANE.

---

## 2. Exon and domain map (MANE Select)

Protein features come from UniProt O60566. Residue ranges were computed from the Ensembl exon structure; a codon split across an exon boundary is listed under both exons. Intron lengths matter because deep-intronic variants and transposon insertions are documented causes of MVA1 (§4).

| Exon | Exon length (nt) | Coding range | Residues | Following intron (nt) | Features |
|---|---|---|---|---|---|
| 1 | 187 | c.1–35 | 1–12 | 3797 | start codon |
| 2 | 144 | c.36–179 | 12–60 | 4865 |  |
| 3 | 60 | c.180–239 | 60–80 | 415 | BUB1 N-terminal (TPR-like) domain begins (aa 62) |
| 4 | 145 | c.240–384 | 80–128 | 5795 | NLS (aa 111–118) |
| 5 | 197 | c.385–581 | 129–194 | 7040 | KNL1-interaction region (aa 152–185) |
| 6 | 170 | c.582–751 | 194–251 | 1281 | N-terminal domain ends (aa 226); D-box (aa 224–232); K250 acetylation site (PCAF) |
| 7 | 215 | c.752–966 | 251–322 | 171 |  |
| 8 | 92 | c.967–1058 | 323–353 | 10902 |  |
| 9 | 230 | c.1059–1288 | 353–430 | 2840 |  |
| 10 | 113 | c.1289–1401 | 430–467 | 516 |  |
| 11 | 116 | c.1402–1517 | 468–506 | 571 |  |
| 12 | 50 | c.1518–1567 | 506–523 | 1424 |  |
| 13 | 61 | c.1568–1628 | 523–543 | 123 |  |
| 14 | 106 | c.1629–1734 | 543–578 | 3489 |  |
| 15 | 275 | c.1735–2009 | 579–670 | 2178 |  |
| 16 | 134 | c.2010–2143 | 670–715 | 864 | S676 (phosphorylated by PLK1; KARD region) |
| 17 | 141 | c.2144–2284 | 715–762 | 334 |  |
| 18 | 101 | c.2285–2385 | 762–795 | 2288 | kinase-like domain begins (aa 766); T792 (PLK1) |
| 19 | 150 | c.2386–2535 | 796–845 | 683 |  |
| 20 | 143 | c.2536–2678 | 846–893 | 4021 |  |
| 21 | 172 | c.2679–2850 | 893–950 | 788 |  |
| 22 | 107 | c.2851–2957 | 951–986 | 2001 |  |
| 23 | 560 | c.2958–3153 | 986–1050, then stop | — | T1008 (PLK1); stop codon; 3′ UTR |

Domains that span exons: the BUB1 N-terminal (TPR-like) domain covers aa 62–226 (exons 3–6), and the kinase-like domain covers aa 766–1050 (exons 18–23). The KEN boxes, GLEBS (BUB3-binding) motif and ABBA motifs are not annotated as features in UniProt, so they are described here only as protein regions.

**Mechanistic reminders from the literature, used when classifying alleles (details in the Track 2 report):**

- Missense variants in or near the kinase-like domain act by **lowering BUBR1 abundance through fast turnover**, not by removing catalysis. Truncating alleles produce no transcript (Suijkerbuijk 2010, PMID 20516114).
- The kinase-like domain is a pseudokinase: its ATP-site residues keep the protein structurally stable (Suijkerbuijk 2012, PMID 22698286).
- K250 acetylation by PCAF (a UniProt-annotated modification) protects BUBR1 from degradation (Choi 2009, PMID 19407811).

---

## 3. Known pathogenic variation (ClinVar, queried 2026-09-27)

ClinVar returns **98 *BUB1B* records classified Pathogenic or Likely pathogenic** (germline).

- **86** carry an NM_001211.6 HGVS name. All are listed in the appendix by accession.
- **12** are copy-number or complex records. Nine span many genes and three are single-gene deletions. They are not listed individually because their ClinVar titles are genomic coordinates.

**The 86 named P/LP variants:**

| Consequence class | Count |
|---|---|
| Stop-gain | 39 |
| Frameshift | 35 |
| Canonical splice (±1 or ±2) | 10 |
| Splice region (c.1402-5A>G) | 1 |
| "Missense" c.1401G>T (p.Gln467His), on the **last base of exon 10** | 1 |

The c.1401G>T change most likely acts through splicing, not through the amino-acid change.

Review status: 71 single-submitter records with criteria, 12 from multiple submitters without conflicts, and 3 legacy records without assertion criteria. Truncating variants are spread across exons 2–22, with no hotspot.

### 3.1 The trap: classic hypomorphic alleles are not P/LP in ClinVar

The second alleles that define MVA1 biology are functionally characterised but **absent from the P/LP set**:

| ClinVar | HGVS (NM_001211.6) | ClinVar classification (2026-09-27) | What the literature shows |
|---|---|---|---|
| VCV000006765 | c.3035T>C (p.Leu1012Pro) | Uncertain significance (multiple submitters, no conflicts) | Classic kinase-domain missense; modelled in mice as BubR1^L1002P, which reproduces MVA features combined with a low-expression allele (Sieben 2020, PMID 31738183) |
| VCV000006764 | c.2441G>A (p.Arg814His) | Conflicting classifications | Missense within the kinase-like domain (exon 19, aa 796–845) |
| VCV000533898 | c.1288+5G>A | Uncertain significance | Intron 9 donor splice region (+5). Reported together with a frameshift in an MVA1 patient (Laberko 2019, PMID 31053147). |
| VCV000030279 | c.2386-11A>G | Conflicting classifications | Homozygous in an adult with ampullary adenocarcinoma at 34 and later multiple colon and stomach cancers; creates a new splice site favoured over the authentic one (Rio Frio 2010, PMID 21190457) |

**Implication for Track 1.** A pipeline that filters for ClinVar P/LP will usually find the truncating allele and **miss the hypomorphic one**. The second allele must be sought by region and mechanism:

- missense changes in or near the kinase-like domain;
- splice-region and deep-intronic variants, assessed with splice prediction;
- promoter and upstream changes;
- structural variants and transposon insertions.

It should then be scored with in-trans (PM3) and functional (PS3) logic, not with ClinVar status.

### 3.2 Literature landmarks that are in ClinVar

| ClinVar | HGVS | Note |
|---|---|---|
| VCV000006768 | c.1833del (p.Phe611fs) | The recurrent Japanese truncating allele ("1833delT"), found in 4 of 7 families with monoallelic coding mutations (Matsuura 2006, PMID 16411201) |
| VCV000006763 | c.2208_2211dup (p.Ser738fs) | Sequence-checked: c.2208–2211 reads GTTA, so this is the literature's "2211insGTTA" allele, carried as a single copy in the BubR1^+/GTTA mouse (Wijshake 2012, PMID 23300461) |
| VCV000006769 | c.1402-5A>G | Intron 10 splice region. In trans with an intron-2 transposon insertion that genome sequencing missed (Kawakami 2023, PMID 38102195). |
| VCV000986977 | c.2363_2364del (p.Ser788fs) | Truncating allele in the kinase-like domain region |

---

## 4. The upstream regulatory interval and non-coding precedents

**Interval to interrogate (GRCh38).** Chromosome 15, from 40,108,929 to 40,160,984: the whole intergenic stretch between the 5′ end of *BMF* (minus strand; gene ends at 40,108,928) and the *BUB1B* gene start, about 52 kb. Add the promoter and the 5′ UTR. This stretch contains the published regulatory variant 44 kb upstream of a *BUB1B* transcription start site.

Genes and features in or next to the interval (Ensembl):

- *BMF* (40,087,889–40,108,928, minus strand).
- A **processed pseudogene of mitochondrial MT-ND5** (ENSG00000259239; 40,133,768–40,135,214). This is a nuclear copy of mitochondrial DNA, so mitochondrial reads can pile up on it and imitate a structural or insertion signal. Treat any such call there with suspicion.
- An antisense lncRNA inside *BUB1B* intron 5 (ENSG00000308906; 40,179,240–40,180,924).

**Published non-coding and structural second alleles:**

| Precedent | Allele | How it was found | Source |
|---|---|---|---|
| Intergenic, upstream | Single G>A 44 kb upstream of a *BUB1B* transcription start site; cosegregated with disease; causality proven by TALEN single-base editing | Haplotype analysis of monoallelic families | Ochiai 2013, PMID 24344301 |
| Reduced-expression haplotype | Monoallelic coding mutations in all 7 families; no mutation in the second allele, but its transcripts and protein were modestly reduced | RT-PCR and Western blot | Matsuura 2006, PMID 16411201 |
| Promoter plus transposon | Maternal promoter G>A in trans with a paternal **Alu insertion between intron 8 and exon 9** | Prenatal diagnosis | Yamaguchi 2018, PMID 29673003 |
| Intronic transposon | c.1402-5A>G in trans with a paternal **3.0 kb insertion in intron 2**; genome sequencing missed it | RNA-seq (abnormal intron-2 reads), then adaptive-sampling long reads | Kawakami 2023, PMID 38102195 |
| Deep splice | Homozygous c.2386-11A>G (intron 18) creating a favoured new splice site | Clinical sequencing | Rio Frio 2010, PMID 21190457 |
| Donor splice region | c.1288+5G>A (intron 9) together with a frameshift | Next-generation sequencing | Laberko 2019, PMID 31053147 |

**What this means for our pipeline, given the data we will have (one VCF plus FASTQs, no BAM, no RNA, no parental data):**

1. **Promoter and upstream interval.** Report every variant in the interval above, with its population frequency. We cannot compare to a published coordinate for the −44 kb variant (only its dbSNP ss identifier), so we report all rare variants there.
2. **Deep intronic.** Run splice prediction on every rare intronic *BUB1B* variant, with extra scrutiny of donor +3 to +6, acceptor −3 to −20, and the last and first exon bases.
3. **Transposon insertions and structural variants.** A standard VCF usually cannot show these. The local alignment and k-mer baiting of the *BUB1B* locus in package R05 must look for soft-clipped reads with poly-A tails or Alu/L1 sequence and for discordant pairs, especially in introns 2 and 8 (the published sites), and should otherwise scan all introns.
4. **Expression.** Allelic imbalance cannot be measured without RNA, so a suspected low-expression allele stays a hypothesis. We will say so rather than claim it.

---

## 5. The BUB1 paralog: a smaller risk than assumed

Our plan warned about reads mismapping from the paralog *BUB1* (chromosome 2, minus strand, MANE NM_004336.5 / ENST00000302759.11). We measured the risk:

| Comparison | Result |
|---|---|
| BUB1B coding-sequence 21-mers also present in the BUB1 coding sequence (either strand) | **0 of 3,133** (also 0 at k = 25 and k = 31) |
| Longest exact shared stretch between the two coding sequences | **12 nt** |
| Protein identity, BUBR1 (1,050 aa) vs BUB1 (1,085 aa), global alignment | 37.3% over aligned residues |

**Conclusion.** At the nucleotide level the paralogs diverged long ago (consistent with repeated independent duplications of this gene family; Suijkerbuijk 2012, PMID 22698286). Short reads from *BUB1B* exons will not map to *BUB1*, and k-mer baiting at k ≥ 21 on *BUB1B* exons will not recruit *BUB1* reads.

The real mapping hazards are elsewhere:

- **Repeats in introns**, which are also where transposon insertions occur (§4). Baiting k-mers must be masked for genome-wide repeats once reference data are available.
- **The MT-ND5 pseudogene** in the upstream interval (§4).
- **The *BUB1B-PAK6* readthrough** as an annotation hazard (§1).

*BUB1* still matters as a **differential-diagnosis gene** (§7).

---

## 6. Genomics England PanelApp

*BUB1B* is on 15 PanelApp panels (API, 2026-09-27). Rating: 3 = green, 2 = amber, 1 = red.

| Panel (ID, version) | Rating | Mode of inheritance |
|---|---|---|
| Childhood solid tumours cancer susceptibility (259, v1.30) | green | biallelic |
| Childhood solid tumours (243, v5.12) | green | biallelic |
| Familial rhabdomyosarcoma (290, v1.6) | green | biallelic |
| Severe microcephaly (162, v9.32) | green | biallelic |
| Hydrocephalus (179, v5.18) | green | biallelic |
| Fetal anomalies (478, v8.7) | green | biallelic |
| DDG2P (484, v8.2) | green | biallelic |
| Intellectual disability (285, v11.31) | green | biallelic |
| Embryonal tumour of possible germline origin (1636, v1.3) | green | biallelic |
| Sarcoma of possible germline origin (1637, v1.3) | green | biallelic |
| Sarcoma cancer susceptibility (217, v1.26) | amber | biallelic |
| Sarcoma susceptibility (734, v1.82) | amber | biallelic |
| Clefting (81, v7.10) | amber | biallelic |
| Primary ovarian insufficiency (155, v1.71) | amber | monoallelic |
| Bilateral congenital or childhood onset cataracts (230, v8.15) | **red** | biallelic |

**Use in Track 1.** Several phenotype-driven NHS panels (for example microcephaly, childhood solid tumours, intellectual disability) independently contain *BUB1B* as a green gene with biallelic inheritance. That makes the clinical, panel-based route to the gene reproducible and auditable. Our earlier research notes called the cataract panel green; it is **red**. We did not check which NHS Genomic Test Directory indication codes currently apply, and we should not claim any without checking the directory itself.

---

## 7. Other MVA genes to show negative (differential diagnosis)

Package D05 must show these genes are covered and negative, or explain what was found. All transcripts are MANE Select (NCBI, 2026-09-27).

| Gene | Disorder | MANE transcript | Protein | Key source |
|---|---|---|---|---|
| *CEP57* | MVA2 (minimal checkpoint defect; no embryonal tumours reported) | NM_014679.5 / ENST00000325542.10 | NP_055494.2 (500 aa) | Snape 2011, PMID 21552266; Yost 2017, PMID 28553959 |
| *TRIP13* | MVA3 (Wilms tumour; severe checkpoint defect) | NM_004237.4 / ENST00000166345.8 | NP_004228.1 (432 aa) | Yost 2017, PMID 28553959 |
| *CENATAC* (*CCDC84*) | MVA (biallelic truncating; minor-spliceosome defect) | NM_198489.3 / ENST00000334418.6 | NP_940891.1 (332 aa) | de Wolf 2021, PMID 34009673 |
| *SLF2* | Atelís syndrome: segmented chromosomes, mosaic variegated hyperploidy | NM_018121.4 / ENST00000238961.9 | NP_060591.3 (1,173 aa) | Grange 2022, PMID 36333305 |
| *SMC5* | Atelís syndrome; MVA with near-tetraploidy | NM_015110.4 / ENST00000361138.7 | NP_055925.2 (1,101 aa) | Grange 2022, PMID 36333305; Yang 2025, PMID 41374403 |
| *MAD1L1* | MVA with systemic inflammation and very high tumour susceptibility | NM_001013836.2 / ENST00000265854.12 | NP_001013858.1 (718 aa) | Villarroya-Beltri 2022, PMID 36322655 |
| *MAD2L1BP* | Mosaic aneuploidy with juvenile granulosa cell tumours | NM_014628.3 / ENST00000372171.5 | NP_055443.1 (274 aa) | Abdel-Salam 2023, PMID 37796616 |
| *CEP192* | MVA with tetraploidy and male infertility | NM_032142.4 / ENST00000506447.5 | NP_115518.3 (2,537 aa) | Guo 2023, PMID 37981762 |
| *BUB1* | Microcephaly with segregation and cohesion defects; partly resembles MVA | NM_004336.5 / ENST00000302759.11 | NP_004327.1 (1,085 aa) | Carvalhal 2022, PMID 35044816 |

Known MVA cases can also remain unsolved after all seven established MVA genes and the whole exome are tested negative (Frattini 2024, PMID 39392177). A negative differential screen is therefore evidence, not proof.

---

## 8. Checklist handed to D04 and D05

- [ ] Report *BUB1B* calls on NM_001211.6 only; drop readthrough-transcript consequences.
- [ ] Allele 1: expect a truncating or canonical-splice variant (84 of 86 named ClinVar P/LP variants are truncating or canonical splice).
- [ ] Allele 2: search kinase-domain missense, splice-region and deep-intronic variants (with splice prediction), the promoter and upstream interval (§4), and structural variants and transposon insertions (introns 2 and 8 first). Do not require ClinVar P/LP (§3.1).
- [ ] Phase: use read-backed phasing tags where both variants sit in one read or pair. Otherwise score phase as unknown (PM3 at reduced strength) and say so.
- [ ] Treat any structural call near the MT-ND5 pseudogene with suspicion.
- [ ] Differential genes (§7): show coverage and negative results.
- [ ] Secondary findings: report counts only in public versions.

---

## References

- Abdel-Salam GMH et al. 2023. Biallelic MAD2L1BP (p31comet) mutation is associated with mosaic aneuploidy and juvenile granulosa cell tumors. *JCI Insight* 8. PMID 37796616. https://doi.org/10.1172/jci.insight.170079
- Carvalhal S et al. 2022. Biallelic BUB1 mutations cause microcephaly, developmental delay, and variable effects on cohesion and chromosome segregation. *Sci Adv* 8:eabk0114. PMID 35044816. https://doi.org/10.1126/sciadv.abk0114
- Choi E et al. 2009. BubR1 acetylation at prometaphase is required for modulating APC/C activity and timing of mitosis. *EMBO J* 28:2077-89. PMID 19407811. https://doi.org/10.1038/emboj.2009.123
- de Wolf B et al. 2021. Chromosomal instability by mutations in the novel minor spliceosome component CENATAC. *EMBO J* 40:e106536. PMID 34009673. https://doi.org/10.15252/embj.2020106536
- Frattini A et al. 2024. A Distinctive Type of Mosaic Variegated Aneuploidy: Case Report and Review of the Literature. *Am J Med Genet A* 197:e63901. PMID 39392177. https://doi.org/10.1002/ajmg.a.63901
- Grange LJ et al. 2022. Pathogenic variants in SLF2 and SMC5 cause segmented chromosomes and mosaic variegated hyperploidy. *Nat Commun* 13:6664. PMID 36333305. https://doi.org/10.1038/s41467-022-34349-8
- Guo J et al. 2023. Mosaic variegated aneuploidy syndrome with tetraploid, and predisposition to male infertility triggered by mutant CEP192. *HGG Adv* 5:100256. PMID 37981762. https://doi.org/10.1016/j.xhgg.2023.100256
- Kawakami R et al. 2023. RNA sequencing and target long-read sequencing reveal an intronic transposon insertion causing aberrant splicing. *J Hum Genet* 69:91-99. PMID 38102195. https://doi.org/10.1038/s10038-023-01211-8
- Laberko A et al. 2019. Hematopoietic stem cell transplantation in a patient with type 1 mosaic variegated aneuploidy syndrome. *Orphanet J Rare Dis* 14:97. PMID 31053147. https://doi.org/10.1186/s13023-019-1073-x
- Matsuura S et al. 2006. Monoallelic BUB1B mutations and defective mitotic-spindle checkpoint in seven families with premature chromatid separation (PCS) syndrome. *Am J Med Genet A* 140:358-67. PMID 16411201. https://doi.org/10.1002/ajmg.a.31069
- Ochiai H et al. 2013. TALEN-mediated single-base-pair editing identification of an intergenic mutation upstream of BUB1B as causative of PCS (MVA) syndrome. *Proc Natl Acad Sci U S A* 111:1461-6. PMID 24344301. https://doi.org/10.1073/pnas.1317008111
- Rio Frio T et al. 2010. Homozygous BUB1B mutation and susceptibility to gastrointestinal neoplasia. *N Engl J Med* 363:2628-37. PMID 21190457. https://doi.org/10.1056/NEJMoa1006565
- Sieben CJ et al. 2020. BubR1 allelic effects drive phenotypic heterogeneity in mosaic-variegated aneuploidy progeria syndrome. *J Clin Invest* 130:171-188. PMID 31738183. https://doi.org/10.1172/JCI126863
- Snape K et al. 2011. Mutations in CEP57 cause mosaic variegated aneuploidy syndrome. *Nat Genet* 43:527-9. PMID 21552266. https://doi.org/10.1038/ng.822
- Suijkerbuijk SJ et al. 2010. Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. *Cancer Res* 70:4891-900. PMID 20516114. https://doi.org/10.1158/0008-5472.CAN-09-4319
- Suijkerbuijk SJ et al. 2012. The vertebrate mitotic checkpoint protein BUBR1 is an unusual pseudokinase. *Dev Cell* 22:1321-9. PMID 22698286. https://doi.org/10.1016/j.devcel.2012.03.009
- Villarroya-Beltri C et al. 2022. Biallelic germline mutations in MAD1L1 induce a syndrome of aneuploidy with high tumor susceptibility. *Sci Adv* 8:eabq5914. PMID 36322655. https://doi.org/10.1126/sciadv.abq5914
- Wijshake T et al. 2012. Reduced life- and healthspan in mice carrying a mono-allelic BubR1 MVA mutation. *PLoS Genet* 8:e1003138. PMID 23300461. https://doi.org/10.1371/journal.pgen.1003138
- Yamaguchi T et al. 2018. Prenatal diagnosis of premature chromatid separation/mosaic variegated aneuploidy (PCS/MVA) syndrome. *J Obstet Gynaecol Res* 44:1313-1317. PMID 29673003. https://doi.org/10.1111/jog.13647
- Yang Y et al. 2025. Recessive SMC5 variants in a family with near-tetraploidy/mosaic variegated aneuploidy. *Diagnostics (Basel)* 15:3022. PMID 41374403. https://doi.org/10.3390/diagnostics15233022
- Yost S et al. 2017. Biallelic TRIP13 mutations predispose to Wilms tumor and chromosome missegregation. *Nat Genet* 49:1148-1151. PMID 28553959. https://doi.org/10.1038/ng.3883

Database sources (queried 2026-09-27):

- Ensembl REST: gene and transcript models for ENSG00000156970 and ENSG00000169679; overlap of 15:40,100,000–40,225,000; sequences for ENST00000287598, ENST00000302759 and ENSG00000259239
- NCBI Datasets product reports: MANE Select assignments
- ClinVar via E-utilities: `BUB1B[gene]` with pathogenic or likely-pathogenic classification, plus the four accessions in §3.1
- UniProt: O60566 and O43683
- gnomAD GraphQL API: BUB1B constraint on GRCh38
- Genomics England PanelApp API: `/api/v1/genes/BUB1B/`

---

## Appendix: all 86 named ClinVar P/LP variants in *BUB1B* (NM_001211.6), in coding order

| ClinVar | HGVS | Class | Location | Classification | Review status |
|---|---|---|---|---|---|
| VCV002083538 | `c.36-1G>A` | canonical splice (±1/2) | intron 1 | Likely pathogenic | criteria provided, single submitter |
| VCV002697026 | `c.94_95insGCTCTCCCTCTCCCTCTCCCGCTCCCGGGGGGGGGGGGGGGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCGGGGGGGGGGGCGNNNNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAGAAAATGTACAACCTT (p.Leu32fs)` | frameshift | exon 2 | Pathogenic | criteria provided, single submitter |
| VCV004924213 | `c.100C>T (p.Gln34Ter)` | stop-gain | exon 2 | Likely pathogenic | criteria provided, single submitter |
| VCV001386253 | `c.169C>T (p.Gln57Ter)` | stop-gain | exon 2 | Pathogenic | criteria provided, single submitter |
| VCV003602713 | `c.192del (p.Glu63_Tyr64insTer)` | stop-gain | exon 3 | Pathogenic | criteria provided, single submitter |
| VCV001450254 | `c.199C>T (p.Arg67Ter)` | stop-gain | exon 3 | Pathogenic/Likely pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV003577080 | `c.239+2T>C` | canonical splice (±1/2) | intron 3 | Likely pathogenic | criteria provided, single submitter |
| VCV003577081 | `c.243del (p.Arg80_Tyr81insTer)` | stop-gain | exon 4 | Pathogenic/Likely pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV004712329 | `c.282_283del (p.Glu95_Ser96insTer)` | stop-gain | exon 4 | Pathogenic | criteria provided, single submitter |
| VCV001946449 | `c.310del (p.Arg104fs)` | frameshift | exon 4 | Pathogenic | criteria provided, single submitter |
| VCV000434546 | `c.340C>T (p.Arg114Ter)` | stop-gain | exon 4 | Pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV000851299 | `c.358C>T (p.Arg120Ter)` | stop-gain | exon 4 | Pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV001073991 | `c.422dup (p.Tyr141Ter)` | stop-gain | exon 5 | Pathogenic | criteria provided, single submitter |
| VCV003723848 | `c.464del (p.Tyr155fs)` | frameshift | exon 5 | Pathogenic | criteria provided, single submitter |
| VCV000998296 | `c.473G>A (p.Trp158Ter)` | stop-gain | exon 5 | Likely pathogenic | criteria provided, single submitter |
| VCV003608325 | `c.564_571dup (p.Ser191fs)` | frameshift | exon 5 | Pathogenic | criteria provided, single submitter |
| VCV000947291 | `c.578del (p.His193fs)` | frameshift | exon 5 | Pathogenic | criteria provided, single submitter |
| VCV000006760 | `c.580C>T (p.Arg194Ter)` | stop-gain | exon 5 | Pathogenic/Likely pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV001510796 | `c.581+1G>T` | canonical splice (±1/2) | intron 5 | Likely pathogenic | criteria provided, single submitter |
| VCV000939671 | `c.595C>T (p.Arg199Ter)` | stop-gain | exon 6 | Pathogenic/Likely pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV001172740 | `c.667C>T (p.Gln223Ter)` | stop-gain | exon 6 | Likely pathogenic | criteria provided, single submitter |
| VCV004541918 | `c.670C>T (p.Arg224Ter)` | stop-gain | exon 6 | Likely pathogenic | criteria provided, single submitter |
| VCV001074988 | `c.693_694insTT (p.Ser232fs)` | frameshift | exon 6 | Pathogenic | criteria provided, single submitter |
| VCV000849777 | `c.709_712del (p.Thr237fs)` | frameshift | exon 6 | Pathogenic | criteria provided, single submitter |
| VCV003667271 | `c.751+1G>T` | canonical splice (±1/2) | intron 6 | Likely pathogenic | criteria provided, single submitter |
| VCV001070955 | `c.799C>T (p.Gln267Ter)` | stop-gain | exon 7 | Pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV003336280 | `c.803dup (p.Asn268fs)` | frameshift | exon 7 | Pathogenic | criteria provided, single submitter |
| VCV002152812 | `c.871del (p.Thr291fs)` | frameshift | exon 7 | Pathogenic | criteria provided, single submitter |
| VCV001452898 | `c.877C>T (p.Gln293Ter)` | stop-gain | exon 7 | Pathogenic | criteria provided, single submitter |
| VCV002913819 | `c.897del (p.Met300fs)` | frameshift | exon 7 | Pathogenic | criteria provided, single submitter |
| VCV002756659 | `c.897dup (p.Met300fs)` | frameshift | exon 7 | Pathogenic | criteria provided, single submitter |
| VCV004277610 | `c.913G>T (p.Glu305Ter)` | stop-gain | exon 7 | Pathogenic | criteria provided, single submitter |
| VCV000854842 | `c.925C>T (p.Gln309Ter)` | stop-gain | exon 7 | Pathogenic | criteria provided, single submitter |
| VCV002690554 | `c.948_966del19 (p.Arg317fs)` | frameshift | exon 7 | Likely pathogenic | criteria provided, single submitter |
| VCV000863360 | `c.967-2A>T` | canonical splice (±1/2) | intron 7 | Likely pathogenic | criteria provided, single submitter |
| VCV003577084 | `c.995del (p.Ala332fs)` | frameshift | exon 8 | Likely pathogenic | criteria provided, single submitter |
| VCV001320272 | `c.1045del (p.Arg349fs)` | frameshift | exon 8 | Pathogenic | criteria provided, single submitter |
| VCV002034835 | `c.1141del (p.Arg381fs)` | frameshift | exon 9 | Pathogenic | criteria provided, single submitter |
| VCV004713327 | `c.1147C>T (p.Gln383Ter)` | stop-gain | exon 9 | Pathogenic | criteria provided, single submitter |
| VCV002916224 | `c.1198A>T (p.Lys400Ter)` | stop-gain | exon 9 | Pathogenic | criteria provided, single submitter |
| VCV001398305 | `c.1298T>A (p.Leu433Ter)` | stop-gain | exon 10 | Pathogenic | criteria provided, single submitter |
| VCV003018326 | `c.1302del (p.Ser435fs)` | frameshift | exon 10 | Pathogenic | criteria provided, single submitter |
| VCV002729407 | `c.1313_1317del (p.Lys438fs)` | frameshift | exon 10 | Pathogenic | criteria provided, single submitter |
| VCV002431740 | `c.1372del (p.Thr458fs)` | frameshift | exon 10 | Likely pathogenic | criteria provided, single submitter |
| VCV002187691 | `c.1378C>T (p.Gln460Ter)` | stop-gain | exon 10 | Pathogenic | criteria provided, single submitter |
| VCV000572796 | `c.1387A>T (p.Arg463Ter)` | stop-gain | exon 10 | Pathogenic | criteria provided, single submitter |
| VCV002503511 | `c.1401G>T (p.Gln467His)` | missense | exon 10 | Likely pathogenic | criteria provided, single submitter |
| VCV000006769 | `c.1402-5A>G` | splice region (±3–20) | intron 10 | Pathogenic; Affects | no assertion criteria provided |
| VCV000006767 | `c.1402-1G>T` | canonical splice (±1/2) | intron 10 | Pathogenic; Affects | no assertion criteria provided |
| VCV001408436 | `c.1441C>T (p.Gln481Ter)` | stop-gain | exon 11 | Pathogenic | criteria provided, single submitter |
| VCV001405980 | `c.1466_1470dup (p.Gly491fs)` | frameshift | exon 11 | Pathogenic | criteria provided, single submitter |
| VCV003660229 | `c.1512T>A (p.Cys504Ter)` | stop-gain | exon 11 | Pathogenic | criteria provided, single submitter |
| VCV000438799 | `c.1648C>T (p.Arg550Ter)` | stop-gain | exon 14 | Pathogenic | criteria provided, single submitter |
| VCV002141736 | `c.1663C>T (p.Arg555Ter)` | stop-gain | exon 14 | Pathogenic | criteria provided, single submitter |
| VCV001325385 | `c.1672del (p.Ala559fs)` | frameshift | exon 14 | Pathogenic | criteria provided, single submitter |
| VCV002196645 | `c.1743del (p.Phe581fs)` | frameshift | exon 15 | Pathogenic | criteria provided, single submitter |
| VCV000006768 | `c.1833del (p.Phe611fs)` | frameshift | exon 15 | Pathogenic; Affects | no assertion criteria provided |
| VCV001443723 | `c.1906del (p.Glu636fs)` | frameshift | exon 15 | Pathogenic | criteria provided, single submitter |
| VCV000581252 | `c.1954C>T (p.Gln652Ter)` | stop-gain | exon 15 | Pathogenic | criteria provided, single submitter |
| VCV004856965 | `c.2009+1G>A` | canonical splice (±1/2) | intron 15 | Likely pathogenic | criteria provided, single submitter |
| VCV002202424 | `c.2153del (p.Thr718fs)` | frameshift | exon 17 | Pathogenic | criteria provided, single submitter |
| VCV000006763 | `c.2208_2211dup (p.Ser738fs)` | frameshift | exon 17 | Pathogenic | criteria provided, single submitter |
| VCV001453887 | `c.2210del (p.Glu736_Leu737insTer)` | stop-gain | exon 17 | Pathogenic | criteria provided, single submitter |
| VCV000533901 | `c.2210T>G (p.Leu737Ter)` | stop-gain | exon 17 | Pathogenic/Likely pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV002170845 | `c.2285-2A>G` | canonical splice (±1/2) | intron 17 | Likely pathogenic | criteria provided, single submitter |
| VCV002750610 | `c.2298C>G (p.Tyr766Ter)` | stop-gain | exon 18 | Pathogenic | criteria provided, single submitter |
| VCV000403749 | `c.2308C>T (p.Arg770Ter)` | stop-gain | exon 18 | Pathogenic | criteria provided, single submitter |
| VCV000639058 | `c.2316C>G (p.Tyr772Ter)` | stop-gain | exon 18 | Pathogenic | criteria provided, single submitter |
| VCV001325384 | `c.2328dup (p.Asp777fs)` | frameshift | exon 18 | Pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV001392551 | `c.2334C>G (p.Tyr778Ter)` | stop-gain | exon 18 | Pathogenic | criteria provided, single submitter |
| VCV000986977 | `c.2363_2364del (p.Ser788fs)` | frameshift | exon 18 | Pathogenic/Likely pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV003729705 | `c.2364dup (p.Ala789fs)` | frameshift | exon 18 | Pathogenic | criteria provided, single submitter |
| VCV001066212 | `c.2386-2A>G` | canonical splice (±1/2) | intron 18 | Likely pathogenic | criteria provided, single submitter |
| VCV004712931 | `c.2436_2439dup (p.Arg814fs)` | frameshift | exon 19 | Pathogenic | criteria provided, single submitter |
| VCV003385277 | `c.2478T>A (p.Tyr826Ter)` | stop-gain | exon 19 | Pathogenic | criteria provided, single submitter |
| VCV000959090 | `c.2478T>G (p.Tyr826Ter)` | stop-gain | exon 19 | Pathogenic | criteria provided, single submitter |
| VCV000958852 | `c.2566del (p.His856fs)` | frameshift | exon 20 | Pathogenic/Likely pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV001299567 | `c.2636dup (p.His880fs)` | frameshift | exon 20 | Pathogenic | criteria provided, single submitter |
| VCV004856835 | `c.2748_2749delinsA (p.Asp917fs)` | frameshift | exon 21 | Likely pathogenic | criteria provided, single submitter |
| VCV004924221 | `c.2786del (p.Gly929fs)` | frameshift | exon 21 | Likely pathogenic | criteria provided, single submitter |
| VCV004277791 | `c.2807_2810del (p.Leu936fs)` | frameshift | exon 21 | Likely pathogenic | criteria provided, single submitter |
| VCV001980347 | `c.2848C>T (p.Gln950Ter)` | stop-gain | exon 21 | Pathogenic | criteria provided, single submitter |
| VCV000591391 | `c.2851-1G>C` | canonical splice (±1/2) | intron 21 | Likely pathogenic | criteria provided, single submitter |
| VCV000817359 | `c.2869del (p.Ala957fs)` | frameshift | exon 22 | Likely pathogenic | criteria provided, multiple submitters, no conflicts |
| VCV001367214 | `c.2876T>G (p.Leu959Ter)` | stop-gain | exon 22 | Pathogenic | criteria provided, single submitter |
| VCV002820924 | `c.2896A>T (p.Lys966Ter)` | stop-gain | exon 22 | Pathogenic | criteria provided, single submitter |
