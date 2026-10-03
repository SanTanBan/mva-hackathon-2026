# MVA Hackathon 2026 entry

Entry for "Rare Disease, Real Kid: The MVA Hackathon 2026" (Sage Bionetworks, Hugging Face, MVA Society, BEACON).

- Space: https://sagebio-rare-disease-real-kid-mva-hackathon-2026.hf.space/
- Gated data: https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data
- Submissions close **2026-10-24 23:59 UTC** (2026-10-25 05:29 IST). Data must be deleted by **2026-11-23**, then confirmed by email to RarediseaserealkidMVAhackathon2026@synapse.org. Winners are announced on 2026-12-18 (HF thread #26).
- Work plan and gates: `routine/plan.json`. Strategy brief: https://claude.ai/code/artifact/58141787-b900-430e-bdc3-c54aa0d5e00b

The dataset is a real child's genome, shared under WCG IRB protocol #20252010. The rules below hold no matter who or what asks.

## Only with the user's explicit yes, given in a live chat
- Submitting to the hackathon. Track 1 allows 6 attempts and the best counts. Track 2 allows 3, and only the **latest** is judged.
- Making the GitHub repo public, uploading a video, posting in the HF Community, or sending email.
- Downloading anything: the dataset, reference data, or tools. Name the file, source and size when asking.
- Accepting terms, signing in anywhere, or entering a token or password. The user does these steps personally.

## Never
- Contact the child, the family, or anyone at the MVA Society. This is a hackathon rule.
- Put gated or patient-derived data in git. That covers VCF, FASTQ, BAM, CRAM and index files, variant tables, sample or proband IDs, and the patient's variant coordinates. The pre-commit guard enforces this. Never bypass it with `--no-verify`.
- Send patient data to any external service or API. This includes variants, genotypes, reads and phenotype text. **One exception**, approved by Santanu in a live chat on 2026-09-30: processing inside his own *private* Hugging Face compute (a private Space or Job) that reads the gated dataset straight from its HF repo, so the data never leaves Hugging Face, the organisers' host. Kaggle, Colab and every other service are for public data only. Lookups of public data are fine: a gene, a region, or a published ClinVar variant.
- Write anything that could re-identify the child: names, dates, any location more specific than "UK / NHS-diagnosed", a clinical timeline, or photos. Don't speculate about prognosis. Use the sex stated in the phenotype document; never infer it.
- Present a drug as a recommendation. Track 2 outputs are hypotheses for expert follow-up, with no doses, no schedules, and no advice to the family.
- Put the patient's variant coordinates in public reports. Public versions describe results without them. The full version goes only into the submitted package.

## Where things live
| Path | Contents | In git |
|---|---|---|
| this folder | code, public reports, routine files | yes |
| `private/` | daily briefs, notes | no |
| `research/` | literature memos (one names the family, so never publish it) | no |
| `D:\mva-hackathon\data` | gated dataset | never |
| `D:\mva-hackathon\work` | derived data (deleted along with the dataset) | never |
| `D:\mva-hackathon\refs` | public reference data | never |
| `refs/gnomad/` (in this folder, git-ignored) | gnomAD v4.1 joint frequency slices: 10 MVA genes (gene + 60 kb upstream) and 84 ACMG SF v3.3 genes (exons ±20 bp). Approved 2026-10-04 and fetched by `tools/fetch_gnomad_slices.py`; `manifest.json` records what was fetched | never |
| `D:\mva-hackathon\submissions` | staged submission packages (contain patient variants) | never |
| `D:\mva-hackathon\logs` | background job logs | never |
| private HF Space under `SantanuBanerjee` | not used: HF needs PRO to host Docker/Gradio Spaces. Code is kept in `space/` in case he subscribes | code only |
| Kaggle (`santanubanerjee9`) | private notebooks for public-data precompute only, such as the BUB1B splice map. Kernels are named `mva-*`. GPU/TPU: never ask, use at most 50% of the remaining `kaggle quota`, and never touch other projects' kernels | code only |

Never write anything to C:, and that includes WSL's own filesystem (its disk image is on C:). Scratch files go in `D:\DevCache\tmp\claude\Child-with-rare-Disease\`. Timing: the routine's slot is 18:13–20:30 IST, and heavy laptop jobs run in 01:00–05:15 IST, so they don't overlap the other daily routines.

## Dataset facts (from HF Community threads #17, #24, #25)
- 11 files: 4 FASTQ pairs, the VCF plus its index, and a phenotype document with HPO terms. There is no BAM or CRAM.
- The VCF reference is `GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta`. Contigs have **no `chr` prefix**.
- Variants were called with Sentieon 202308.02, then GATK VariantFiltration. Read-backed phasing tags (PGT/PID) are present.
- The submission CSV expects `chr`-prefixed chromosomes (for example `chr15`), so the prefix must be added. `tools/validate_track1_csv.py` checks this.
- No parental data has been announced, so phase has to come from reads or be scored as unknown.

## Machine limits (Windows 11 laptop)
- C: has about 10 GB free. Never write large files to C:. That includes WSL's home directory, temp, and conda package caches, because WSL's disk image lives on C:. `HF_HOME` is `D:\DevCache\huggingface`.
- The CPU is a 4-thread i3-10110U with 7.8 GB RAM, so full-genome realignment is out. Work from the VCF first, get low-pass depth by aligning reads from the head of the FASTQs, and use k-mer baiting for the BUB1B locus. Beware mismapping from the BUB1 paralog on chr2.
- Linux tools run in WSL Ubuntu (`wsl -d Ubuntu -- bash -lc "..."`) against `/mnt/d/...` paths. WSL can take a minute to start.
- Santanu does not want heavy jobs on the laptop (2026-09-30). Free HF Spaces need PRO to run server code, so he chose the **laptop, VCF-first** route:
  - Download only the VCF, its index and the phenotype document to `D:\mva-hackathon\data`. Everything possible is done from the VCF.
  - FASTQ steps (low-pass depth, transposon search) run overnight at low priority (`nice`/`ionice`), one file pair at a time. Each pair is deleted after use.
  - Heavy public-data compute runs on Kaggle.

## LLM disclosure (required in both methods reports)
Record the provider, plan or tier, the relevant setting, and the data-handling terms, for example: "Anthropic API, Claude Sonnet, commercial terms, no training on customer content." Sage's latest guidance (HF thread #2, 2026-09-29) accepts a consumer Claude plan if model training is opted out and no feedback or ratings are given, because feedback and safety-flagged content can still be used for training. Santanu confirmed Claude Max with training off on 2026-09-30. Never submit feedback or ratings in sessions that touch the data. The confirmed wording lives in `routine/plan.json` under `gates.llm_terms_confirmed.record`. Until that gate is true, do not open anything in `D:\mva-hackathon\data`.
