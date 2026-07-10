# PC Task: GSE189926 Annotation and QC Refinement

Created: 2026-07-10

## Goal

Refine `GSE189926` immune cell labels and audit QC/sample structure before any communication analysis.

Do not run CellChat, LIANA, NicheNet, differential communication, final statistical testing, or manuscript figures in this task.

## Inputs

Use the PC-local object recorded in:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_response_object_construction/object_inventory.tsv`

Use the prior response object construction outputs from:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_response_object_construction`

Use this manifest:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/20260710_gse189926_annotation_qc_refinement_manifest.tsv`

## Output Directory

Upload lightweight outputs to:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260710_gse189926_annotation_qc_refinement`

Large h5ad/RDS objects should remain PC-local and be listed in `object_inventory.tsv` if modified or newly created.

## Required Upload Files

General files:

- `STATUS.md`
- `run_status.tsv`
- `warning_log.tsv`
- `object_inventory.tsv`
- `session_info.txt`

`GSE189926` required files:

- `gse189926_refined_label_audit.tsv`
- `gse189926_cluster_marker_top50_raw.tsv`
- `gse189926_cluster_marker_top50_no_technical.tsv`
- `gse189926_qc_flag_summary.tsv`
- `gse189926_qc_flag_by_sample.tsv`
- `gse189926_sample_patient_mixing_audit.tsv`
- `gse189926_refined_celltype_counts_by_response_timepoint.tsv`
- `gse189926_refined_label_by_cluster_sample.tsv`
- `gse189926_refined_umap_qc.png`
- `gse189926_refined_marker_dotplot.png`

`GSE235863` required files:

- `gse235863_cd45_patient_level_composition.tsv`
- `gse235863_cd45_subcluster_patient_proportions.tsv`
- `gse235863_cd45_response_composition_qc.png`
- `gse235863_cd8_patient_level_doublet_sensitivity.tsv`
- `gse235863_cd8_doublet_excluded_composition_qc.png`

## Required Status Columns

`run_status.tsv`:

- dataset_id
- step
- status
- input_files
- output_files
- cell_count_input
- cell_count_used_for_tables
- notes
- error_message

`warning_log.tsv`:

- dataset_id
- warning_type
- affected_field
- affected_count
- explanation
- action_taken

`object_inventory.tsv`:

- dataset_id
- object_name
- pc_local_path
- object_format
- file_size_bytes
- sha256
- upload_status
- notes

## GSE189926 Requirements

### 1. Preserve Raw Fields

Preserve all existing raw fields:

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

Do not overwrite broad labels from the previous task. Add refined labels in new fields.

### 2. QC Flag Audit

Do not remove cells by default.

Create QC flags and summarize how many cells are affected:

- low detected genes
- high detected genes
- high total counts
- high mitochondrial percentage
- high ribosomal percentage
- high hemoglobin percentage

Use explicit numeric thresholds and write them to `gse189926_qc_flag_summary.tsv`.

Also upload:

- cell counts by QC flag and sample
- cell counts by QC flag and raw outcome
- cell counts by QC flag and broad cell label

### 3. Marker Export

Export two marker tables for the working clustering field:

1. `gse189926_cluster_marker_top50_raw.tsv`
2. `gse189926_cluster_marker_top50_no_technical.tsv`

The technical-gene-excluded table should remove obvious mitochondrial genes, ribosomal genes, hemoglobin genes, and MALAT1-like technical high-abundance genes before selecting top markers. The removed gene pattern must be written in `warning_log.tsv` or `run_status.tsv`.

### 4. Refined Immune Label Audit

Refine labels only when marker evidence supports it. Suggested label space:

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

Use `Unresolved_low_signal` when the marker evidence is not strong enough. Do not force labels.

Marker evidence fields in `gse189926_refined_label_audit.tsv`:

- cluster
- broad_label_previous
- refined_label
- supporting_marker_genes
- conflicting_marker_genes
- mean_expression_summary
- fraction_expressing_summary
- label_confidence_level
- label_decision_note

### 5. Sample and Patient Structure Audit

The previous UMAP showed strong sample and patient structure. Quantify it instead of relying only on visual inspection.

Upload `gse189926_sample_patient_mixing_audit.tsv` with:

- cluster
- refined_label
- total_cells
- patient_count
- largest_patient
- largest_patient_fraction
- sample_count
- largest_sample
- largest_sample_fraction
- outcome_count
- largest_outcome_fraction
- timepoint_count
- largest_timepoint_fraction

Flag clusters with high single-patient or single-sample dominance.

### 6. Readable QC Figures

Regenerate compact QC figures:

- `gse189926_refined_umap_qc.png`
- `gse189926_refined_marker_dotplot.png`

Figure requirements:

- no large blank panel
- no clipped axis labels or titles
- marker dotplot split into readable marker groups if needed
- show broad and refined labels separately
- include raw outcome and timepoint panels

These are still QC figures, not manuscript figures.

## GSE235863 Requirements

### 1. CD45 Patient-Level Composition

Create patient-level composition tables rather than only pooled cell-level summaries:

- patient
- raw_response
- derived_response_group
- response_ambiguity_flag
- tissue
- major_cluster
- sub_cluster
- cell_count
- proportion_within_patient_tissue

Preserve the `P18` ambiguity flag.

### 2. CD8 Doublet Sensitivity

Summarize CD8 subcluster proportions before and after predicted doublet exclusion at patient level:

- patient
- raw_response
- derived_response_group
- tissue
- doublet_filter_state
- sub_cluster
- cell_count
- proportion_within_patient_tissue

Do not overwrite raw h5ad objects.

## Scientific Guardrails

- Do not infer response labels for `GSE236581`.
- Do not remove `NE` samples from raw tables.
- Do not pool cell-level counts as independent biological replicates in any statistical statement.
- Do not force fine labels when marker evidence is weak.
- Preserve raw labels alongside refined labels.
- Record package versions in `session_info.txt`.

## Completion Signal

After uploading results, commit and push to the same branch:

`pan-gi-ici-metadata-request-20260709`

Commit message:

`Add Pan-GI ICI annotation QC refinement outputs`
