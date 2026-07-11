# PC Task: GSE235863 Sample-Timepoint Rebuild and PRJNA932556 CRC Raw Processing

Created: 2026-07-11

## Goal

Complete two independent stages:

1. Rebuild GSE235863 sample-level metadata and patient-level composition without collapsing tissue or treatment timepoint.
2. Validate the PRJNA932556 raw-read structure and, only if technically supported, generate six CRC single-cell count matrices and an R/Seurat QC object.

Do not run CellChat, LIANA, NicheNet, differential communication, manuscript figures, or final cross-cohort statistical testing in this task.

## Required Language

Use R for metadata analysis, composition analysis, Seurat processing, QC, doublet assessment, clustering, marker export, and annotation.

SRA Toolkit and Cell Ranger or another explicitly validated read-processing tool may be used for raw-read download and quantification. Record exact tool versions and commands. Do not use Scanpy for downstream single-cell analysis.

## Stage A: GSE235863 Sample and Timepoint Rebuild

### Inputs

- PC-local `GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad.gz`.
- PC-local `GSE235863_five_patients_scRNAseq_cd8t_raw_counts.h5ad.gz`.
- `uploads/20260709_metadata_download/GSE235863_family.soft.gz`, if present locally from the metadata download inventory.
- Existing response and composition tables under `uploads/20260709_preprocess_priority_response`, `uploads/20260709_response_object_construction`, and `uploads/20260710_gse189926_r_only_rerun`.

### Mandatory source-conflict audit

GEO row `GSM7510911` has conflicting source fields:

- sample title: `Pre-treatment blood of patient P27 (R/PR)`
- `characteristics_ch1::patient`: `P18`
- response: `PR, responder`

Do not overwrite the raw GEO fields. Resolve the analytical mapping only when the object `sample` value, GEO accession/title, neighboring source rows, and any article supplement or object metadata provide an exact link. Record every observed field and the evidence used. If the object-to-GEO link remains unresolved, exclude only the conflicting sample from analytical summaries and report it.

### Required mapping table

Create `gse235863_sample_mapping_audit.tsv` with these exact columns:

- `dataset_id`
- `object_sample`
- `object_patient_raw`
- `source_gsm`
- `source_sample_title`
- `source_patient_raw`
- `source_response_raw`
- `source_tissue_raw`
- `derived_patient`
- `derived_response_group`
- `derived_tissue`
- `derived_timepoint`
- `mapping_status`
- `mapping_evidence`
- `exclusion_reason`

Do not infer a mapping from capitalization, substring similarity, row order, or an undocumented code. Unresolved rows must remain unresolved.

### Required composition tables

Create tables that retain patient, sample, tissue, and timepoint:

- `gse235863_cd45_sample_level_composition.tsv`
- `gse235863_cd45_patient_timepoint_composition.tsv`
- `gse235863_cd45_patient_timepoint_subcluster.tsv`
- `gse235863_cd8_sample_level_doublet_sensitivity.tsv`
- `gse235863_cd8_patient_timepoint_doublet_sensitivity.tsv`
- `gse235863_analysis_stratum_counts.tsv`

The CD45 tables must include `major_cluster` and `sub_cluster`. The CD8 tables must include `sub_cluster` and `doublet_filter_state`.

Required analysis strata:

- pre-treatment liver tumor
- post-treatment liver tumor
- pre-treatment blood
- post-treatment blood
- paired pre/post within liver tumor
- paired pre/post within blood

Never combine pre-treatment and post-treatment rows or blood and tumor rows to increase sample size.

### Exploratory statistics

Create `gse235863_patient_level_composition_statistics.tsv` with one row per stratum and cell label. Include patient counts by response group, median proportions, median difference, exact test name, raw p value, adjusted p value, and analysis status.

If either response group has fewer than 3 patients, set the inferential fields to missing, retain descriptive values, and set `analysis_status` to `descriptive_only_group_n_below_3`.

Create compact QC figures with all patient points visible:

- `gse235863_cd45_timepoint_composition_qc.png`
- `gse235863_cd8_timepoint_doublet_qc.png`

These are QC figures, not manuscript figures.

## Stage B: PRJNA932556 Raw-Read Feasibility and Processing

### Exact source rows

Use the six runs and response mapping already recorded in:

- `uploads/20260709_metadata_download/PRJNA932556_sra_runinfo.csv`
- `uploads/20260709_small_file_download/response_metadata__PRJNA932556_sample_response_mapping.tsv`

The recorded run sizes total 52,709 MB. Process sequentially to avoid exhausting disk, memory, or I/O.

### Mandatory pilot before full download

Use one run as a pilot. Do not assume the 10x chemistry or barcode/UMI layout from average read length alone.

Required pilot outputs:

- `prjna932556_read_structure_audit.tsv`
- `prjna932556_pilot_command_log.txt`
- `prjna932556_pilot_fastq_metrics.tsv`
- `prjna932556_quantification_feasibility.tsv`

Record the exact observed read files, read lengths, pairing, headers, barcode/UMI evidence, cDNA read evidence, chemistry evidence, quantification tool, reference build, and decision.

If barcode/UMI structure or chemistry cannot be validated, stop Stage B after the pilot and set the decision to `blocked_read_structure_not_validated`. Do not force bulk RNA-seq alignment and do not fabricate a cell-by-gene matrix.

### Full processing only after the pilot passes

For each of S1, S2, S3, R1, R2, and R3:

1. Download and verify the run.
2. Convert to FASTQ while preserving pairing and read identifiers.
3. Quantify with the validated single-cell workflow and exact reference build.
4. Record raw and filtered matrix dimensions, reads assigned to cells, median genes per cell, median UMI per cell, and any tool warnings.
5. Keep SRA, FASTQ, Cell Ranger output, and large objects PC-local.

Create:

- `prjna932556_run_processing_status.tsv`
- `prjna932556_cellranger_or_equivalent_metrics.tsv`
- `prjna932556_matrix_inventory.tsv`
- `prjna932556_r_qc_summary.tsv`
- `prjna932556_r_doublet_summary.tsv`
- `prjna932556_r_cluster_resolution_summary.tsv`
- `prjna932556_r_marker_top50_no_technical.tsv`
- `prjna932556_r_label_audit.tsv`
- `prjna932556_r_patient_celltype_counts.tsv`
- `prjna932556_r_umap_qc.png`
- `prjna932556_r_marker_dotplot.png`

Build the downstream object in R/Seurat. Preserve sample name and raw response labels. Do not pool cells before QC, doublet assessment, or sample-level summaries.

## General Required Files

- `STATUS.md`
- `run_status.tsv`
- `warning_log.tsv`
- `object_inventory.tsv`
- `session_info.txt`
- `SHA256SUMS.txt`

Large files must remain PC-local and be recorded in `object_inventory.tsv` with absolute PC paths, file sizes, and SHA256 values.

## Output Directory

Upload lightweight outputs to:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260711_gse235863_timepoint_prjna932556`

## Completion Signal

Commit and push to branch:

`pan-gi-ici-metadata-request-20260709`

Commit message:

`Add Pan-GI ICI timepoint rebuild and CRC raw processing audit`
