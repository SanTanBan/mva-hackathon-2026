# BUBR1 stability map: interpretation rules, fixed before the run

**Rare Disease, Real Kid: The MVA Hackathon 2026 · package R13 · rules written 2026-10-23, before the map was computed · public version**

> **Why this exists.** The D08 worksheet classifies an untested *BUB1B* substitution by position alone: anything at residue 727 or later is "abundance class (provisional)", because every substitution tested there was unstable (Suijkerbuijk 2010, PMID 20516114). A structure-based stability prediction can test that provisional call. To keep the test honest, we fix here, before computing anything, how the prediction will be read and what it may change.

## What is computed (public data only)

The private Kaggle notebook `kaggle/mva_bub1b_stability_map/` scores every possible amino-acid substitution in BUBR1 (UniProt O60566, 1,050 residues; 19,950 substitutions). No patient data enters it; the child's substitution is looked up afterwards, locally.

- **ΔΔG prediction:** ThermoMPNN (Dieckhaus 2024, PNAS 121:e2314853121; MIT licence) on the AlphaFold DB model of O60566. The sign is calibrated inside the run so that positive means destabilising: buried hydrophobic residues replaced by charged residues or proline must score as destabilising on average.
- **Sequence score:** ESM-2 650M masked-marginal log-likelihood ratio (Lin 2023, Science 379:1123; MIT licence). Windows of 1,022 residues are centred on each position, because the protein is longer than the model's context.
- **Structure features:** AlphaFold pLDDT, relative solvent accessibility and contact count for each residue.

## Validity gate (checked before the child's substitution is looked up)

The published substitutions serve as controls (Suijkerbuijk 2010):
- **Abundance class (unstable):** R727C, R814H, L844F, I909T, L1012P.
- **No defect detected (stable):** R550Q, Q921H.
- **Functional class (stable, defective):** Y155C, which is reported but not used in the gate.

**The map counts as informative for *BUB1B*** only if at least 3 of the 5 abundance-class controls have a higher predicted ΔΔG than both no-defect controls. Otherwise the report says the predictor does not separate the tested classes for this protein, and the child's substitution is not interpreted with it.

## Reading the child's substitution (only if the gate passes)

The reference set is every possible substitution at residues 727–1,050 with pLDDT ≥ 70.

| Predicted ΔΔG of the child's substitution | Call |
|---|---|
| At or above the 75th percentile of the reference set | **supports** the abundance class |
| At or below the 25th percentile, at a solvent-exposed residue (relative accessibility > 0.25) | **argues against** the abundance class: a stable protein, possibly functional class |
| Anything else | **inconclusive** |

The ESM-2 score is reported alongside as context, but it does not change the call.

## What the call may change

Computation does not move the ranking. The pre-registered rules in `reports/track2/d08_worksheet.md` (step 4) key the ranking on the mechanism class, and only assay 1 (BUBR1 protein in patient-derived cells) settles the class.
- **"Supports":** the report says the provisional abundance class has independent computational support. The ranking stays H1–H5.
- **"Argues against":** the report flags the functional-class possibility. It states the worksheet's contingency (H3 first, then H1, H2, H5) as the expected order if assay 1 shows a stable protein. The ranking stays H1–H5 until then.
- **"Inconclusive" or a failed gate:** the report says so, and the class stays positional.

The public version of the report gives the call and the controls. It does not give the child's substitution or its scores; those appear only in the submitted version.
