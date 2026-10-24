# Rare Disease, Real Kid: MVA Hackathon 2026 entry (SantanuBanerjee)

Code, pre-registered rules and public reports for both tracks of *Rare Disease, Real Kid: The MVA Hackathon 2026* (Sage Bionetworks, Hugging Face, the MVA Society and BEACON). The case is a child with mosaic variegated aneuploidy type 1 (MVA1), caused by biallelic *BUB1B* variants.

- **Track 1, variant prediction.** A VCF-first, inheritance-aware ranking with a rule fixed before data access. The top-ranked *BUB1B* compound-heterozygous pair matched the clinically confirmed answer at rank 1 (100/100 rank points, F-max 1.000).
- **Track 2, drug repurposing.** Five ranked repurposing hypotheses for expert follow-up, with safety filters for a cancer-predisposed child, re-ranked with the child's alleles under rules fixed in advance.

**No patient data is in this repository.** The gated dataset (a real child's genome, shared under WCG IRB protocol #20252010) was processed only in private compute: private Kaggle notebooks named `mva-*`. The public reports describe results without variant names, coordinates or identifiers; the full versions went only to the organisers. A pre-commit guard (`tools/guard_commit.py`, run by `.githooks/pre-commit`) refuses VCF lines, variant rows and protected identifiers. All copies of the data are deleted within 30 days of the close, as the data-use terms require.

## Reports

| Track | Public report | Supporting documents |
|---|---|---|
| 1 | [reports/track1/SantanuBanerjee_track1_report.md](reports/track1/SantanuBanerjee_track1_report.md) (also PDF) | [locus dossier](reports/track1/locus_dossier.md), [methods template](reports/track1/methods_report.md) |
| 2 | [reports/track2/SantanuBanerjee_track2_report.md](reports/track2/SantanuBanerjee_track2_report.md) (also PDF) | [allele-classification rules](reports/track2/d08_worksheet.md), [stability-map rules](reports/track2/stability_map_rules.md), [pitch script](reports/track2/pitch_script.md), [methods form](reports/track2/methods_description.md) |

Pitch video: https://youtu.be/OpTEsDaAQnM

## Rules fixed before the results

Each interpretive rule was committed before the data or result it applies to existed, so the git history is the record:

- the Track 1 ranking rule and evidence weights (`pipeline/track1_rank.py`, methods §B9), before data access;
- the PanelApp gene list (`tools/build_panel_genes.py`), before data access;
- the mechanism classes, antisense tiers and Track 2 re-ranking rules (`reports/track2/d08_worksheet.md`, `pipeline/allele_mechanism.py`), before data access;
- the reading of the BUBR1 stability map (`reports/track2/stability_map_rules.md`), before the map was computed.

## Layout

| Path | What it is |
|---|---|
| `pipeline/` | Analysis modules: annotation table, compound-het pairing, Track 1 ranking, allele mechanism classes, splice-map and stability-map lookups, B-allele fraction and low-pass depth, ACMG points, the ACMG SF v3.3 screen, the drug-triage pipeline. Tests in `pipeline/tests/`. |
| `tools/` | Command-line steps: data intake, annotation (`annotate_vcf.sh`), Kaggle helpers, the Track 1 CSV validator, report drafting and rendering, the commit guard. |
| `kaggle/` | Kaggle notebooks. Public data only: `mva_refs_build` (references and tests), `bub1b_splice_map` (Pangolin map of every *BUB1B* SNV), `mva_bub1b_stability_map` (ThermoMPNN and ESM-2 for every BUBR1 substitution), `mva_pitch_video`. Patient data, run privately: `mva_d02` (Track 1 fast path, SF screen), `mva_d03` (low-pass depth). `mva_d02_synth` is the same pipeline on a synthetic VCF. |
| `reports/` | Public reports and the documents above; `reports/track2/stability_map/` holds the public stability map. |
| `routine/plan.json` | The work plan, gates and dated notes kept during the hackathon. |
| `space/` | An unused Hugging Face Space variant of the pipeline. |

## Reproducing

1. **References and tests (public data):** run `kaggle/mva_refs_build`. It fetches the GRCh38 analysis set, ClinVar and gnomAD v4.1 slices, checks their md5 sums, and runs the unit tests and the end-to-end test on a synthetic VCF (`tools/test_annotate_e2e.py`).
2. **Synthetic dry run:** `kaggle/mva_d02_synth` runs the Track 1 pipeline end to end on public alleles only.
3. **Gated data** (requires approved access to `SageBio/mva-hackathon-2026-data`): `tools/fetch_vcf_first.py` → `tools/d01_intake.py` → `tools/kaggle_patient_upload.py` → `tools/kaggle_d02.py` (runs `kaggle/mva_d02`) → `tools/validate_track1_csv.py` → `tools/draft_track1_report.py` → `tools/render_reports.py`. Low-pass depth: `kaggle/mva_d03`, started from the Kaggle website with the user's own read-only Hugging Face token as a Kaggle secret.
4. **Public maps:** `kaggle/bub1b_splice_map` and `kaggle/mva_bub1b_stability_map`; look a variant up locally with `python -m pipeline.splice_map` or `python -m pipeline.stability_map`.

Unit tests: `python -m pytest pipeline/tests` and the `tools/test_*.py` scripts.

## Generative AI

Anthropic Claude Max (consumer plan), used through Claude Code; model-training opt-out enabled; no feedback or ratings given in sessions involving the data.

## Licence

Code and reports: CC BY 4.0 (see [LICENSE](LICENSE)), as the hackathon releases participant outputs. Third-party tools (Pangolin, ThermoMPNN, ESM-2, bcftools, minimap2) are fetched at run time under their own licences and are not included.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.
