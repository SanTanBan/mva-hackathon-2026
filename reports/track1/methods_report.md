# Track 1 methods report (skeleton)

**Rare Disease, Real Kid: The MVA Hackathon 2026 · Track 1 · skeleton v0, 2026-09-27**

> **How this file works.** This is the public template. Once analysis starts, the working copy with results lives outside the repository, at `D:\mva-hackathon\work\reports\track1\methods_report.src.md`. From that copy, `python tools/render_reports.py <source> --track track1 --name methods_report` writes:
> - this public file, with every block and inline span marked private withheld (the marker syntax is documented in `tools/render_reports.py`);
> - the full version, into `D:\mva-hackathon\submissions\track1\<date>\`.
>
> Both are produced as Markdown and PDF. The renderer refuses to write a public version that contains variant coordinates, VCF lines or protected identifiers.
>
> **Section order is provisional.** It will be aligned with the organisers' methods template, including its required LLM-usage question, once that template has been downloaded (setup package U3).

---

## 1. Summary

*Two or three sentences: the diagnosis we support, the evidence for each* BUB1B *allele, and what is still uncertain.*

<!-- private:start -->
*Submitted version only: the ranked candidate pairs with HGVS names and evidence weight.*
<!-- private:end -->

## 2. Team and roles

*Name, Hugging Face username (SantanuBanerjee), and who did what. State where an LLM agent was used (see §9).*

## 3. Data, access and governance

- Dataset: whole-genome FASTQ (four pairs), one VCF with its index, and a phenotype document with HPO terms, shared under WCG IRB protocol #20252010 through gated access on Hugging Face.
- Reference build: GRCh38 (no-alt analysis set with decoys); contigs have no `chr` prefix. Variants were called with Sentieon, then filtered with GATK VariantFiltration.
- Handling: gated and derived data stay under `D:\mva-hackathon\`, never in git, and are deleted by 2026-11-23. No patient data was sent to any external service.

## 4. Methods

### 4.1 Compute environment
One Windows 11 laptop (4 threads, 7.8 GB RAM) with Linux tools in WSL Ubuntu. Tool names, versions and commands will be listed here and in the repository. Full-genome realignment was out of reach, so the analysis works from the VCF plus targeted read-level checks.

### 4.2 VCF quality control and normalisation
Header, sample and build checks; left-alignment and multi-allelic splitting; quality filters; `chr` prefix added on output.

### 4.3 Annotation and inheritance-aware filtering
- Rare-variant filters using population frequency.
- Recessive models: homozygous, and compound heterozygous pairs.
- Consequence calls pinned to MANE Select transcripts.

### 4.4 Phenotype-driven prioritisation
- HPO-based ranking.
- An independent clinical route through Genomics England PanelApp panels, where *BUB1B* is a green, biallelic gene on several phenotype panels (locus dossier §6).

### 4.5 *BUB1B* locus deep-dive
Both alleles are sought by mechanism, not by ClinVar status: the classic hypomorphic alleles are VUS or conflicting in ClinVar (locus dossier §3.1). The search covers:
- the promoter and the upstream interval;
- splice-region and deep-intronic variants, with splice prediction;
- structural variants and transposon insertions, via k-mer baiting of the locus;
- read-backed phase.

### 4.6 Evidence for mosaic aneuploidy
- Per-chromosome B-allele fraction from VCF allele depths.
- Low-pass depth from alignments of reads taken from the start of the FASTQs.
- Stated detection limits, and the caveat that variegated mosaicism is spread across different chromosomes in different cells.

### 4.7 Differential genes and secondary findings
- Coverage and negative results for the other MVA genes (locus dossier §7).
- The ACMG SF v3.3 screen, reported as counts only in the public version.

### 4.8 ACMG/AMP classification
Points-based classification shown per allele, with the evidence codes and their strengths.

### 4.9 Submission file
How candidate pairs and their scores (`epcr`) were chosen. Validated with `tools/validate_track1_csv.py`.

## 5. Results

<!-- private:start -->
*Submitted version only: variant tables, both alleles with HGVS names, phase evidence, ACMG points, mosaicism figures, differential-gene coverage, and secondary-finding details.*
<!-- private:end -->

Public summary: *what was found, described without coordinates or identifiers.*

## 6. Interpretation

*How the genotype fits MVA1 mechanisms. Classify each allele as functional class or abundance class (Suijkerbuijk 2010, PMID 20516114). Link to the Track 2 report.*

## 7. Limitations

*No BAM or RNA; no parental data; low-pass depth only; what the VCF cannot show (for example transposon insertions); the limits of the mosaicism assay.*

## 8. Reproducibility

*Repository URL, commit hash, environment file, and the exact commands, in order. Public data sources and their versions.*

## 9. LLM usage (required)

*The organisers require the provider, plan or tier, relevant setting and data-handling terms. The confirmed wording is pending: it will be copied from `routine/plan.json` → `gates.llm_terms_confirmed.record` once that gate is set. No patient data was used with an LLM before that gate was confirmed.*

## 10. Privacy and data deletion

*No names, dates, locations more specific than "UK (NHS-diagnosed)", or clinical timeline appear. The public version withholds all patient-level variants. Data deletion is due by 2026-11-23, then confirmed by email.*

---

References: see `reports/track1/locus_dossier.md` for the full, PubMed-checked reference list.
