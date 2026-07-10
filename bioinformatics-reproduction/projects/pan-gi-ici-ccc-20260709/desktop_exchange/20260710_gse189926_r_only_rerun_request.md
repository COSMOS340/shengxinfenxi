# PC Task: R-only GSE189926 Rebuild and Annotation

Created: 2026-07-10

## Decision

Redo the current `GSE189926` task entirely in R.

This does not mean restarting the whole Pan-GI ICI project. It means rebuilding the failed `GSE189926` UMAP/QC/annotation step from the raw 22 matrices using R/Seurat, then comparing unintegrated, QC-pass, and sample-aware corrected embeddings.

Do not run CellChat, LIANA, NicheNet, differential communication, final statistical testing, or manuscript figures in this task.

## Required Language

Use R for the full `GSE189926` rebuild.

Allowed R packages include:

- Seurat
- SeuratObject
- Matrix
- data.table or readr
- dplyr
- ggplot2
- patchwork
- harmony, if available
- SingleR or celldex, only if available and clearly recorded

Do not use Scanpy, Python UMAP, or Python-based clustering for `GSE189926` in this task.

If a small helper is needed for file copying or checksums, that is acceptable, but analysis, QC, UMAP, clustering, marker export, and labeling must be R-based.

## Inputs

Use the raw PC-local `GSE189926` matrices recorded in:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_matrix_download/matrix_file_inventory.tsv`

Use sample metadata from:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_preprocess_priority_response/gse189926_sample_qc.tsv`

Use sample parse audit from:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_preprocess_priority_response/gse189926_sample_parse_audit.tsv`

Use this manifest:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/20260710_gse189926_r_only_rerun_manifest.tsv`

Do not use the previous Scanpy-generated h5ad as the analysis source for `GSE189926`. It may be used only for comparison in notes if helpful.

## Output Directory

Upload lightweight outputs to:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260710_gse189926_r_only_rerun`

Large Seurat objects must remain PC-local and be listed in `object_inventory.tsv`.

## Required Upload Files

General:

- `STATUS.md`
- `run_status.tsv`
- `warning_log.tsv`
- `object_inventory.tsv`
- `session_info.txt`

`GSE189926` R rebuild:

- `gse189926_r_import_audit.tsv`
- `gse189926_r_feature_alignment_audit.tsv`
- `gse189926_r_qc_summary.tsv`
- `gse189926_r_qc_flag_by_sample.tsv`
- `gse189926_r_integration_strategy_audit.tsv`
- `gse189926_r_mixing_metrics.tsv`
- `gse189926_r_cluster_resolution_summary.tsv`
- `gse189926_r_marker_top50_raw.tsv`
- `gse189926_r_marker_top50_no_technical.tsv`
- `gse189926_r_refined_label_audit.tsv`
- `gse189926_r_label_by_cluster_sample.tsv`
- `gse189926_r_celltype_counts_by_response_timepoint.tsv`
- `gse189926_r_umap_unintegrated_qc.png`
- `gse189926_r_umap_qc_pass.png`
- `gse189926_r_umap_harmony.png`
- `gse189926_r_marker_dotplot.png`

`GSE235863` R summaries:

- `gse235863_cd45_r_patient_level_composition.tsv`
- `gse235863_cd45_r_subcluster_patient_proportions.tsv`
- `gse235863_cd45_r_response_composition_qc.png`
- `gse235863_cd8_r_patient_level_doublet_sensitivity.tsv`
- `gse235863_cd8_r_doublet_excluded_composition_qc.png`

## Required Status Columns

`run_status.tsv`:

- dataset_id
- step
- status
- language
- package_route
- input_files
- output_files
- cell_count_input
- cell_count_used
- notes
- error_message

`object_inventory.tsv`:

- dataset_id
- object_name
- pc_local_path
- object_format
- file_size_bytes
- sha256
- upload_status
- notes

Use `pc_local_only_large_file` for large RDS objects.

`warning_log.tsv`:

- dataset_id
- warning_type
- affected_field
- affected_count
- explanation
- action_taken

## GSE189926 R Workflow Requirements

### 1. Raw Matrix Import

Read all 22 gzip text matrices in R.

Requirements:

- import as sparse matrices
- do not create a dense full expression matrix
- prefix cell IDs with `sample_accession`
- preserve all raw metadata fields from `gse189926_sample_qc.tsv`
- align features across samples using the union of observed feature names
- record per-sample feature counts and missing feature counts in `gse189926_r_feature_alignment_audit.tsv`

If a file cannot be imported without dense expansion, stop that sample and report the exact error. Do not silently downsample.

### 2. Seurat Object Construction

Build a Seurat object from the raw counts.

Required metadata:

- `sample_accession`
- `sample_title`
- `patient_id`
- `timepoint`
- `timepoint_simple`
- `characteristics_ch1::outcome`
- `characteristics_ch1::treatment`
- `characteristics_ch1::mmr`
- `characteristics_ch1::cell type`
- `characteristics_ch1::tissue`

Preserve `NE` as a raw outcome label. Do not include `NE` in binary response comparisons.

### 3. QC Audit

Compute at least:

- `nCount_RNA`
- `nFeature_RNA`
- mitochondrial percentage
- ribosomal percentage
- hemoglobin percentage

Create QC flags:

- low detected genes
- high detected genes
- high total counts
- high mitochondrial percentage
- high ribosomal percentage
- high hemoglobin percentage

Do not remove cells from the unfiltered object. Create a QC-pass subset for visualization and marker review.

Upload:

- `gse189926_r_qc_summary.tsv`
- `gse189926_r_qc_flag_by_sample.tsv`

The QC summary must include exact thresholds.

### 4. Normalization and Dimensionality Reduction

Run at least two R/Seurat routes:

1. unintegrated QC-pass route
2. sample-aware corrected route

Preferred route:

- `NormalizeData`
- `FindVariableFeatures`
- remove technical genes from PCA-driving variable features if they dominate
- `ScaleData`
- `RunPCA`
- `FindNeighbors`
- `FindClusters`
- `RunUMAP`
- Harmony correction by `sample_accession` if `harmony` is available

If using SCTransform, also upload a clear note in `gse189926_r_integration_strategy_audit.tsv`. Do not only upload SCT output without a log-normalized route unless memory makes the log-normalized route fail.

Technical genes to evaluate for exclusion from PCA-driving feature set:

- mitochondrial genes
- ribosomal genes
- hemoglobin genes
- MALAT1-like high-abundance genes
- immunoglobulin constant-region genes

Keep these genes in the object and raw marker table. Exclude them only from the PCA-driving feature set when needed.

### 5. UMAP Comparison

Upload three figures:

1. `gse189926_r_umap_unintegrated_qc.png`
2. `gse189926_r_umap_qc_pass.png`
3. `gse189926_r_umap_harmony.png`

Each figure should include panels for:

- sample accession
- patient
- raw outcome
- timepoint
- broad label
- refined label

If Harmony is unavailable, replace `gse189926_r_umap_harmony.png` with the best available R-based sample-aware correction and record the method in `gse189926_r_integration_strategy_audit.tsv`.

Figure requirements:

- no large blank panel
- no clipped labels
- legends outside plotting areas
- consistent colors for the same labels across related UMAPs

### 6. Mixing Metrics

Upload `gse189926_r_mixing_metrics.tsv` with one row per embedding route:

- route_name
- k_neighbors
- cells_used
- mean_same_sample_neighbor_fraction
- median_same_sample_neighbor_fraction
- mean_same_patient_neighbor_fraction
- median_same_patient_neighbor_fraction
- interpretation

The goal is to quantify whether UMAP is dominated by sample or patient structure. Do not erase all patient biology blindly.

### 7. Clustering and Marker Tables

Run multiple clustering resolutions, at minimum:

- 0.3
- 0.5
- 0.8

Upload `gse189926_r_cluster_resolution_summary.tsv`.

Export:

- `gse189926_r_marker_top50_raw.tsv`
- `gse189926_r_marker_top50_no_technical.tsv`

The technical-gene-excluded marker table should remove mitochondrial, ribosomal, hemoglobin, MALAT1-like, and immunoglobulin constant-region genes before selecting top markers.

### 8. Refined Immune Labels

Assign refined labels only when marker evidence supports the label.

Use this label space if supported:

- CD4_T_naive_memory
- CD4_Treg
- CD8_T_GZMK
- CD8_T_cytotoxic
- NK_or_ILC
- B_naive_memory
- Plasma
- Monocyte_classical
- Monocyte_FCGR3A
- Macrophage_C1QC_APOE
- Macrophage_SPP1
- cDC1
- cDC2
- pDC
- DC_LAMP3
- Mast
- Cycling
- Unresolved_low_signal

Do not force labels. Use `Unresolved_low_signal` when evidence is weak.

Upload `gse189926_r_refined_label_audit.tsv` with:

- cluster
- resolution
- broad_label
- refined_label
- supporting_marker_genes
- conflicting_marker_genes
- marker_summary
- label_confidence_level
- label_decision_note

### 9. Marker Dotplot

Upload `gse189926_r_marker_dotplot.png`.

Requirements:

- readable marker groups
- not one very wide compressed strip
- split into panels if needed
- include T, NK/ILC, B/plasma, monocyte/macrophage/DC, mast, and cycling marker groups

### 10. Response and Timepoint Counts

Upload `gse189926_r_celltype_counts_by_response_timepoint.tsv`.

Required columns:

- raw outcome
- binary response group
- timepoint
- timepoint_simple
- broad_label
- refined_label
- cell_count
- proportion_within_outcome_timepoint

Binary response helper rule:

- `PR` -> responder
- `SD` and `PD` -> non_responder
- `NE` -> not_in_binary_response

Preserve raw outcome labels.

## GSE235863 R Summary Requirements

Do not rerun heavy expression processing for `GSE235863`.

Use existing uploaded tables or read metadata in R from the PC-local h5ad if convenient.

Create patient-level composition outputs:

- `gse235863_cd45_r_patient_level_composition.tsv`
- `gse235863_cd45_r_subcluster_patient_proportions.tsv`
- `gse235863_cd45_r_response_composition_qc.png`
- `gse235863_cd8_r_patient_level_doublet_sensitivity.tsv`
- `gse235863_cd8_r_doublet_excluded_composition_qc.png`

Preserve raw response labels and the `P18` ambiguity flag.

## Stop Conditions

Stop and upload a partial `STATUS.md` plus `warning_log.tsv` if:

- the raw matrices cannot be imported sparsely in R
- Seurat object creation fails
- Harmony or other R correction fails and no R-based fallback is available
- memory requires downsampling for visualization

If downsampling is used only for a figure, report retained counts and keep the full object unchanged.

## Completion Signal

After uploading outputs, commit and push to the same branch:

`pan-gi-ici-metadata-request-20260709`

Commit message:

`Add Pan-GI ICI GSE189926 R-only rerun outputs`
