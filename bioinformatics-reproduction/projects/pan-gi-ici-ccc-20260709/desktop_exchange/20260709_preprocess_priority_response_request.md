# PC Task: Priority Response Dataset Preprocessing

Created: 2026-07-09

## Goal

Run preprocessing and QC for the first formal Pan-GI ICI cell communication analysis set.

This task is preprocessing only. Do not run CellChat, LIANA, NicheNet, differential communication, manuscript figures, or final statistical testing in this task.

## Inputs

Use the PC-local matrix files recorded in:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_matrix_download/matrix_file_inventory.tsv`

Also use:

- `bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_matrix_download/download_status.tsv`
- `bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_matrix_download/lightweight_inspection.tsv`
- `/Volumes/research/pan_gi_ici_ccc_20260709/99_logs/pc_tasks/20260709_preprocess_priority_response_manifest.tsv`

Use existing matrix files from the PC local paths in `matrix_file_inventory.tsv`. Do not redownload files unless a checksum verification fails or a file is missing.

## Output Directory

Upload lightweight outputs to:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_preprocess_priority_response`

Large processed objects should remain PC-local and be listed in `processed_object_inventory.tsv`.

## Required Upload Files

Upload these files:

- `STATUS.md`
- `preprocess_status.tsv`
- `checksum_verification.tsv`
- `processed_object_inventory.tsv`
- `dataset_qc_summary.tsv`
- `sample_cell_counts.tsv`
- `response_timepoint_counts.tsv`
- `celltype_counts.tsv`
- `response_join_audit.tsv`
- `response_grouping_rules.tsv`
- `warning_log.tsv`
- `session_info.txt`

Upload small QC figures if generated:

- `gse189926_import_qc_overview.png`
- `gse235863_cd45_umap_overview.png`
- `gse235863_cd8_umap_overview.png`
- `liver_ma_metadata_overview.png`
- `omix001073_reference_counts.png`

Do not upload large h5ad, RDS, MTX, or full expression tables to GitHub.

## Required Status Columns

`preprocess_status.tsv`:

- dataset_id
- step
- status
- input_files
- output_files
- cell_count_before_filtering
- cell_count_after_filtering
- gene_count
- notes
- error_message

`processed_object_inventory.tsv`:

- dataset_id
- object_name
- pc_local_path
- object_format
- file_size_bytes
- sha256
- upload_status
- notes

Use `pc_local_only_large_file` for large processed objects.

`warning_log.tsv`:

- dataset_id
- warning_type
- affected_field
- affected_count
- explanation
- action_taken

## Dataset-Specific Instructions

### 1. GSE189926

Input:

- 22 `GSE189926` gzip text matrices.
- Sample metadata: `/Volumes/research/pan_gi_ici_ccc_20260709/00_metadata/extracted_samples/GSE189926_samples.tsv`.

Required handling:

- Map each matrix to `sample_accession` using the GSM accession in the file name.
- Preserve the raw columns exactly:
  - `sample_title`
  - `characteristics_ch1::outcome`
  - `characteristics_ch1::treatment`
  - `characteristics_ch1::mmr`
  - `characteristics_ch1::cell type`
  - `characteristics_ch1::tissue`
- Parse `patient_id` and `timepoint` from `sample_title`, then upload a parse audit table.
- Preserve `NE` as its own outcome label. Do not include `NE` in binary response comparisons unless we explicitly approve that choice later.
- If a binary response helper column is generated, write the rule to `response_grouping_rules.tsv`.

Suggested binary rule for audit only:

- `PR` -> `responder`
- `SD` and `PD` -> `non_responder`
- `NE` -> `not_evaluable`

Required outputs:

- `gse189926_sample_qc.tsv`
- `gse189926_cell_count_by_sample.tsv`
- `gse189926_cell_count_by_outcome_timepoint.tsv`
- `gse189926_sample_parse_audit.tsv`
- `gse189926_response_grouping_rules.tsv`
- `gse189926_import_status.tsv`
- optional PC-local processed h5ad or RDS object listed in `processed_object_inventory.tsv`

Memory note:

If full dense import of the 22 text matrices is too large, use a streaming sparse import route and record the exact failed command, memory error, and partial output status. Do not silently downsample.

### 2. GSE235863 CD45

Input:

- `GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad.gz`.
- Sample metadata: `/Volumes/research/pan_gi_ici_ccc_20260709/00_metadata/extracted_samples/GSE235863_samples.tsv`.

Required handling:

- Read the h5ad object and preserve existing metadata fields, including `patient`, `sample`, `tissue`, `major_cluster`, `sub_cluster`, doublet-related fields, QC fields, and UMAP coordinates.
- Join sample-level response metadata from `characteristics_ch1::response`.
- Output a join audit showing matched, unmatched, duplicated, and ambiguous records.
- Preserve raw response strings. A derived grouping column may be generated only if the mapping is written to `response_grouping_rules.tsv`.
- Keep patient `P18` visible in the audit because the GEO metadata contains both `CR, responder` and `PR, responder` entries for this patient.

Suggested grouping rule for audit only:

- `CR, responder`, `PR, responder`, and `R` -> `responder`
- `PD, non-responder`, `SD, non-responder`, `NR`, and `NR/SD` -> `non_responder`

Required outputs:

- `gse235863_cd45_obs_summary.tsv`
- `gse235863_cd45_response_join_audit.tsv`
- `gse235863_cd45_celltype_counts.tsv`
- `gse235863_cd45_qc_summary.tsv`
- `gse235863_response_grouping_rules.tsv`
- optional PC-local processed h5ad or RDS object listed in `processed_object_inventory.tsv`

### 3. GSE235863 CD8

Input:

- `GSE235863_five_patients_scRNAseq_cd8t_raw_counts.h5ad.gz`.

Required handling:

- Run the same metadata preservation and response join audit as for CD45.
- Summarize `sub_cluster`, tissue, patient, response, and doublet/QC fields.
- Do not merge the CD8 object into the CD45 object in this task.

Required outputs:

- `gse235863_cd8_obs_summary.tsv`
- `gse235863_cd8_response_join_audit.tsv`
- `gse235863_cd8_subcluster_counts.tsv`
- `gse235863_cd8_qc_summary.tsv`

### 4. ICB_Zenodo_Liver_Ma

Input:

- `All_Celltypes_seurat_Liver_seurat_set_1_2_epi_subset_revised.RDS`.

Required handling:

- Read metadata only unless full object loading is already trivial.
- Summarize `cell_types`, `Treatment`, `outcome`, `pre_post`, and `Cancer_type`.
- Keep HCC and iCCA labels separate.

Required outputs:

- `liver_ma_metadata_summary.tsv`
- `liver_ma_celltype_outcome_counts.tsv`
- `liver_ma_pre_post_counts.tsv`

### 5. OMIX001073

Input:

- `OMIX001073-20-03.zip`
- `OMIX001073-20-12.zip`
- `OMIX001073-20-14.zip`
- Annotation summaries already uploaded in the small-file task.

Required handling:

- Verify `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` dimensions for all three archives.
- Compare barcode counts against available annotation tables.
- Summarize myeloid and stromal cluster counts.
- Do not run full integration with ICI response datasets in this task.

Required outputs:

- `omix001073_matrix_dimension_audit.tsv`
- `omix001073_annotation_overlap.tsv`
- `omix001073_myeloid_stroma_counts.tsv`

### 6. GSE236581

Input:

- `GSE236581_CRC-ICB_metadata.txt.gz`
- `GSE236581_counts.mtx.gz`
- `GSE236581_features.tsv.gz`
- `GSE236581_barcodes.tsv.gz`

Required handling:

- Summarize metadata only unless this is fast.
- Confirm again that no response-related column exists.
- Summarize `Treatment`, `Tissue`, `MajorCellType`, and `SubCellType`.
- Keep this dataset out of response comparisons in this task.

Required outputs:

- `gse236581_metadata_summary.tsv`
- `gse236581_treatment_tissue_celltype_counts.tsv`
- `gse236581_response_column_audit.tsv`

## Scientific Guardrails

- Do not remove clusters or samples in this task.
- Do not remove doublets unless the raw and filtered counts are both reported and the filtering rule is written to `dataset_qc_summary.tsv`.
- Do not downsample without a written reason and exact retained cell counts.
- Do not infer response labels for `GSE236581`.
- Preserve raw labels even when derived labels are generated for convenience.
- Record all package versions in `session_info.txt`.

## Completion Signal

After uploading results, commit and push to the same branch:

`pan-gi-ici-metadata-request-20260709`

Commit message:

`Add Pan-GI ICI priority response preprocessing outputs`
