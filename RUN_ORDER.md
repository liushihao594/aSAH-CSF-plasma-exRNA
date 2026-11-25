# Run order

Run from the repository root in a disposable working copy. Never point the scripts at the only copy of source data.

## 0. Prerequisites

1. Replace the provisional license and citation authors only after author approval.
2. Acquire and validate the source inputs described in `metadata/INPUT_FILE_CONTRACTS.csv`.
3. Recreate the original stage-relative data layout locally, or make path-only adaptations that do not change model parameters.
4. Use the R/package versions in `environment/R_sessionInfo_stage1b.txt` where feasible.

## 1. Pairing and source QC

Run `scripts/01_sample_reconstruction/build_paired_qc_sets.ps1` with explicit input/output paths. This is the public, parameterized copy of the pairing/QC construction logic. Its outputs are local sample-level files and remain ignored by Git.

## 2. Original models and Stage 1 inputs

1. `scripts/02_original_models/legacy_minimal_check.R`
2. `scripts/03_compositional_analysis/stage1_analysis.R pilot`
3. `scripts/04_random_slope_and_QC/finalize_after_pilot_gate.R`

The legacy script retains original stage-relative input names to document the exact reproduction. It requires locally obtained source RPM tables and the provenance crosswalk.

## 3. Boundary diagnostics and full M3/M4 refit

1. `run_stage1a_boundary_diagnostics.R`
2. `write_stage1a_reports.R`
3. `validate_stage1a_outputs.R`
4. Run `run_stage1b_worker.R` for the fixed task grid/checkpoints.
5. `finalize_stage1b_results.R`
6. `build_stage1b_reports.R`
7. `validate_stage1b_outputs.R`

The full historical Stage 1B run comprised 6,344 model tasks and 64 checkpoints. Do not choose optimizers, features or random-effect structures by significance.

## 4. Composition and biofluid analyses

1. `scripts/05_biotype_analysis/run_stage2_composition_analysis.R`
2. `scripts/05_biotype_analysis/validate_stage2_outputs.R`
3. Optional feasibility-only check: `xdec_technical_checks_no_deconvolution.R`; do not interpret this as XDec deconvolution.
4. `scripts/06_biofluid_trajectories/run_stage3_biofluid_analysis.R`
5. `scripts/06_biofluid_trajectories/validate_stage3_outputs.R`

## 5. Cross-layer analysis

1. `scripts/07_cross_layer_analysis/run_stage4_crosslayer_analysis.R`
2. `regenerate_stage4_figure_b.R` if the diagnostic figure must be regenerated.
3. `validate_stage4_outputs.R`

The primary bootstrap used 2,000 patient-cluster resamples. The two annotation layers share libraries and the result is a technical-consistency analysis, not independent biological validation.

## 6. Figures and tables

1. `scripts/08_figures_and_tables/stage5_build.R`
2. `scripts/08_figures_and_tables/validate_stage5.R`

## 7. Revision supplements

Run `scripts/09_revision_supplements/generate_s8_s9.py` only after Stage 1B outputs and the local paired-QC file exist. The script creates model-implied within-subject correlation summaries and mapped-read/QC summaries. Do not commit its sample-level intermediate output.

## 8. Compare expected results

Compare aggregate results with `expected_outputs/EXPECTED_KEY_RESULTS.csv` and the validation CSV files. A successful static check is not the same as a complete rerun. Record the R version, package versions, random seeds, input checksums and output checksums for any independent rerun.
