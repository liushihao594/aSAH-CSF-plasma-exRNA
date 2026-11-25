# aSAH paired CSF–plasma exRNA analysis

Publication candidate for the analysis code supporting:

> **Longitudinal extracellular RNA gradients between cerebrospinal fluid and plasma after aneurysmal subarachnoid hemorrhage**

Journal: BMC Neurology  
Manuscript ID: `8fbcb807-d9ec-4860-a22b-e3f021f14a4d`

**Status:** author-approved metadata have been incorporated into this local release candidate. It has not been committed, pushed, tagged, released, or archived by Zenodo. A separate explicit approval is required before the first public commit and push.

## Research overview

The revision reconstructs longitudinal paired CSF–plasma extracellular-RNA gradients from public study GSE121868. The primary D01–D07 analysis contains 30 patients, 153 paired CSF–plasma observations and 306 biological samples. The fixed testing sets contain 507 miRNAs and 1,079 GENCODE-annotated exRNA features derived from the same underlying small-RNA sequencing libraries.

The original paired gradient is:

```text
log2(CSF RPM + 1) - log2(plasma RPM + 1)
```

## Analytical workflow

The revised analysis includes:

- cohort reconstruction and complete-pair QC rules;
- original limma/duplicateCorrelation models (M1);
- fixed-reference CLR compositional sensitivity (M2/M4);
- patient random-intercept/random-slope models (M3/M4);
- paired QC PASS-only and CR2 patient-cluster sensitivity analyses;
- high-abundance RNA and RNA-biotype composition analyses;
- separate CSF/plasma trajectories and biofluid-by-label interactions;
- endpoint-complete and patient-influence analyses;
- cross-layer LOSO scores and patient-cluster bootstrap;
- final figure/table generation and validation;
- revision supplements for model-implied within-subject correlations and source-reported mapped-read summaries.

## Important scientific boundaries

- The available analysis input is RPM, not validated feature-level raw counts.
- TMM–voom was **not** performed, and RPM was never rescaled or rounded into pseudo-counts.
- No reliable sample-level XDec cargo proportions were recovered or estimated.
- RNA-biotype composition is not equivalent to physical vesicle, lipoprotein or RNP carrier composition.
- miRNA and GENCODE annotations arise from shared sequencing libraries and are not independent biological replication.
- Collection labels D01–D07 are ordered identifiers, not verified exact post-bleed days.
- The 16 miRNA and 33 GENCODE strict-core features are an internally robust subset, not an externally validated biomarker panel.

## Public data source and redistribution

Source study records are indexed at GEO accession [GSE121868](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121868) and BioProject PRJNA498685. This repository deliberately does not redistribute full RPM matrices, FASTQ/BAM files, analysis-ready RDS objects, detailed patient crosswalks, CFDE/exRNA Atlas exports, or third-party XDec materials.

See `DATA_AVAILABILITY.md` and `metadata/INPUT_FILE_CONTRACTS.csv`. Users must obtain source assets legally, verify checksums and create the local inputs expected by the scripts. Private/local inputs belong under `data_private/`, which is ignored by Git.

## Software environment

The finalized Stage 1B run used R 4.6.1 on Windows 11, with key versions recorded in `environment/R_sessionInfo_stage1b.txt`, including limma 3.68.5, lme4 2.0-6, lmerTest 3.2-1, clubSandwich 0.7.0, statmod 1.5.2 and ggplot2 4.0.3.

The release candidate was statically audited on 25 November 2026. An R executable was not available in the packaging shell, so the reorganized candidate was not claimed to have passed a fresh end-to-end rerun. The copied scripts are the paper-specific scripts used for the completed stage analyses; path-only release adaptations and new packaging helpers are identified in the manifest.

## Repository layout

```text
scripts/01_sample_reconstruction/   path-parameterized pairing/QC reconstruction
scripts/02_original_models/         original limma/duplicateCorrelation reproduction
scripts/03_compositional_analysis/  fixed-reference CLR and Stage 1 preparation
scripts/04_random_slope_and_QC/     boundary diagnostics, M3/M4, QC and CR2
scripts/05_biotype_analysis/        high-abundance/biotype analyses and XDec feasibility only
scripts/06_biofluid_trajectories/   CSF/plasma and interaction analyses
scripts/07_cross_layer_analysis/    LOSO scores, bootstrap and influence analyses
scripts/08_figures_and_tables/      final figure/table build and validation
scripts/09_revision_supplements/    Supplementary Tables S8/S9 generation
metadata/                           input contracts and non-identifying constants only
validation/                         non-sensitive stage validation summaries
expected_outputs/                   aggregate expected results only
```

## Models and strict-core rule

- **M1:** original RPM gradient with limma + duplicateCorrelation.
- **M2:** CLR gradient with the original repeated-measures framework.
- **M3:** original RPM gradient with `centered_label + (1 + centered_label || patient)`.
- **M4:** CLR gradient with the same random-effect structure.

Strict-core features had to be significant in M1, M3 and M4; maintain a concordant slope direction; have evaluable random-slope inference in M3 and M4; and show no material direction dependence on the paired QC PASS-only sensitivity branch. This within-cohort intersection does not provide a separate joint-FDR guarantee.

## Running the analysis

Read `RUN_ORDER.md` before running anything. The stage scripts preserve the model parameters, random seeds and decision rules used for the paper. Several scripts expect the original stage-relative directory structure and locally reconstructed source data. Do not silently substitute another processing version or infer counts from RPM.

## Validation level

- Historical stage validation: available and included as non-sensitive CSV summaries.
- Candidate file/path/secret/size checks: performed.
- Candidate R syntax/end-to-end execution: not rerun in the packaging shell because R was unavailable.
- Fresh public-download-to-final-figures reconstruction: not automated; the provenance/QC normalization and some historical comparison inputs must be rebuilt locally as described in `DATA_AVAILABILITY.md`.
- Publication, GitHub push, release and Zenodo archive: not performed.

**End-to-end execution not independently rerun in the release preparation environment.**

## Citation and license

The software authors are Ruina Cui, Zhao Gao and Shihao Liu, in that order. Citation metadata are provided in `CITATION.cff`; no release date or DOI is included before a real release and archive exist.

Original code in this repository is released under the MIT License, copyright (c) 2026 Ruina Cui. This license does not relicense public study data, R packages, XDec materials or other third-party resources; see `THIRD_PARTY_NOTICES.md`.
