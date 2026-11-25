# Data availability and local input preparation

## Public source

- GEO: GSE121868
- BioProject: PRJNA498685
- Related exRNA/CFDE project identifier used during provenance auditing: EXR-KJENS1aSAHPR-ST

The two processed expression files used by the analysis were identified as:

- `GSE121868_KJENS1-aSAH_Project-2016-10-13_exceRpt_miRNA_ReadsPerMillion.txt.gz`
- `GSE121868_KJENS1-aSAH_Project-2016-10-13_exceRpt_gencode_ReadsPerMillion.txt.gz`

Users should obtain files from the public source records, verify that the study/project identifiers and feature annotations match, and retain source checksums. The provenance audit also used file-to-biosample-to-subject relations and source-reported mapped-read/QC fields from the CFDE/exRNA Atlas project. Portal availability and exported filenames can change, so the project identifier and relationship fields—not a filename prefix alone—must be used to restrict the download to this study.

## Not redistributed in this repository

The repository excludes raw sequencing data, complete RPM matrices, counts, RDS/RData objects, detailed sample crosswalks, patient-level metadata, CFDE/exRNA Atlas exports and XDec reference/input files. Public indexing does not by itself establish unrestricted redistribution rights for every downloaded or derived file.

## Local-only preparation

After legally obtaining the source assets, place them under `data_private/` or in the original stage-relative input structure described in `metadata/INPUT_FILE_CONTRACTS.csv`. `data_private/` is excluded from Git. Verify study/project identity, feature annotations, sample/library identity and checksums before analysis.

Authentic feature-level raw counts and matched library-size denominators for the 306 analysis samples were not validated. Consequently, this code release does not include or claim a TMM–voom workflow. Source-reported transcriptome-mapped reads are QC descriptors only and must not be substituted for feature-count library sizes or RPM denominators.

## Reproduction boundary

The repository documents the exact analysis scripts and expected input contracts, but it is not currently a one-command reconstruction from a fresh public download. In particular:

- the detailed file-to-library-to-biosample crosswalk and source-QC export are not redistributed and must be rebuilt from the public relationship tables;
- `prepared_stage1_inputs.rds` is a derived local object created after the verified pairing/QC preparation and is not redistributed;
- the historical cross-layer comparison tables are not redistributed, although the Stage 4 analysis code is provided;
- final figure/table scripts consume validated outputs from preceding stages and therefore require those stages to have been completed locally.

The public-source assets make reconstruction feasible in principle, but the acquisition and provenance-normalization step is not fully automated in this release candidate. The included historical validation summaries support inspection of expected results; they are not a substitute for a fresh end-to-end rerun.
