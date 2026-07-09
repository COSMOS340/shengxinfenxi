# PC Task: Response-Ready Object Construction and Cell Type Audit

Created: 2026-07-09

## Goal

Build response-ready analysis objects and audit tables for the first formal Pan-GI ICI cell-state analysis.

This task must not run CellChat, LIANA, NicheNet, differential communication, final statistical testing, or manuscript figures. The purpose is to create defensible objects, QC summaries, UMAP overviews, marker tables, and cell type label audit files.

## Input Context

Use the PC-local matrices recorded in:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_matrix_download/matrix_file_inventory.tsv`

Use the previous preprocessing outputs from:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_preprocess_priority_response`

Use this manifest:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/20260709_response_object_construction_manifest.tsv`

Do not redownload input files unless a checksum verification fails.

## Output Directory

Upload lightweight outputs to:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_response_object_construction`

Large objects must remain PC-local and be listed in `object_inventory.tsv`.

## Required Upload Files

General files:

- `STATUS.md`
- `run_status.tsv`
- `object_inventory.tsv`
- `warning_log.tsv`
- `session_info.txt`

`GSE189926` files:

- `gse189926_sparse_import_audit.tsv`
- `gse189926_feature_alignment_audit.tsv`
- `gse189926_qc_summary.tsv`
- `gse189926_sample_qc_summary.tsv`
- `gse189926_cluster_resolution_summary.tsv`
- `gse189926_cluster_marker_top20.tsv`
- `gse189926_celltype_label_audit.tsv`
- `gse189926_celltype_counts_by_response_timepoint.tsv`
- `gse189926_umap_overview.png`
- `gse189926_marker_dotplot.png`

`GSE235863` files:

- `gse235863_cd45_response_celltype_counts.tsv`
- `gse235863_cd45_response_subcluster_counts.tsv`
- `gse235863_cd45_tissue_response_counts.tsv`
- `gse235863_cd45_umap_overview.png`
- `gse235863_cd45_abundance_overview.png`
- `gse235863_cd8_response_subcluster_counts.tsv`
- `gse235863_cd8_doublet_sensitivity_counts.tsv`
- `gse235863_cd8_umap_overview.png`
- `gse235863_cd8_abundance_overview.png`

Reference and exclusion files:

- `omix001073_reference_marker_support.tsv`
- `omix001073_reference_annotation_summary.tsv`
- `gse236581_response_exclusion_note.tsv`

## Required Status Columns

`run_status.tsv`:

- dataset_id
- step
- status
- input_files
- output_files
- cell_count_input
- cell_count_in_object
- gene_count
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

Use `pc_local_only_large_file` for objects not uploaded to GitHub.

`warning_log.tsv`:

- dataset_id
- warning_type
- affected_field
- affected_count
- explanation
- action_taken

## Dataset-Specific Requirements

### 1. GSE189926 Full Sparse Import

Inputs:

- 22 gzip text matrices from the PC-local `GSE189926` matrix folder.
- `uploads/20260709_preprocess_priority_response/gse189926_sample_qc.tsv`.
- `uploads/20260709_preprocess_priority_response/gse189926_sample_parse_audit.tsv`.

Required handling:

- Import all 22 text matrices as a sparse object. Do not create a dense full matrix.
- Prefix cell barcodes with `sample_accession` or another exact sample identifier so all cell IDs are globally unique.
- Preserve raw sample metadata exactly, including:
  - `sample_accession`
  - `sample_title`
  - `patient_id`
  - `timepoint`
  - `characteristics_ch1::outcome`
  - `characteristics_ch1::treatment`
  - `characteristics_ch1::mmr`
  - `characteristics_ch1::cell type`
  - `characteristics_ch1::tissue`
- Create a derived `timepoint_simple` field for overview tables, but preserve the raw `timepoint` field. `post-treatment1` and `post-treatment2` should map to `post-treatment` only in the derived field.
- Preserve `NE` as its own outcome label. Do not include `NE` in binary response comparisons.
- If a derived binary response field is created, write its exact rule to `gse189926_celltype_counts_by_response_timepoint.tsv` metadata notes or `warning_log.tsv`.

Feature alignment:

- Audit gene identifiers across all 22 files.
- If features differ by sample, align using the union of gene names and report per-sample missing feature counts in `gse189926_feature_alignment_audit.tsv`.
- Record whether row names are gene symbols, gene IDs, or another exact observed identifier format. Do not guess identifier type if not directly observable.

QC:

- Compute per-cell `n_counts`, `n_genes`, `pct_mt`, `pct_rp`, and `pct_hb` if gene symbols allow those calculations.
- Do not remove cells by default.
- If a QC-pass flag is generated, keep it as a metadata flag and upload the exact threshold rule. Do not overwrite the unfiltered object.

Dimensionality reduction and clustering:

- Normalize counts and calculate highly variable genes.
- Run PCA, neighborhood graph, UMAP, and clustering.
- Run at least three clustering resolutions, for example `0.3`, `0.5`, and `0.8`, or the closest equivalent in the chosen framework.
- Upload `gse189926_cluster_resolution_summary.tsv` showing cluster counts for each resolution.
- Choose one working resolution for marker export and explain the choice in `run_status.tsv`.

Marker and cell type label audit:

- Export top 20 marker genes per cluster to `gse189926_cluster_marker_top20.tsv`.
- Assign immune cell type labels using marker evidence and write the evidence to `gse189926_celltype_label_audit.tsv`.
- Do not use labels unsupported by marker evidence.
- Use `Other_or_low_signal` for clusters that cannot be defensibly labeled.

Marker groups to inspect:

- T cells: `CD3D`, `CD3E`, `TRAC`, `IL7R`, `CCR7`, `LTB`
- CD4 T cells: `CD4`, `IL7R`, `CCR7`, `ICOS`, `CTLA4`
- Treg: `FOXP3`, `IL2RA`, `CTLA4`, `IKZF2`
- CD8 T cells: `CD8A`, `CD8B`, `GZMK`, `GZMB`, `NKG7`, `PRF1`
- NK or ILC: `NKG7`, `GNLY`, `KLRD1`, `KLRF1`, `FCGR3A`
- B cells: `MS4A1`, `CD79A`, `CD79B`, `BANK1`
- Plasma cells: `MZB1`, `JCHAIN`, `XBP1`, `IGHG1`
- Monocytes: `LYZ`, `S100A8`, `S100A9`, `FCN1`, `LST1`, `FCGR3A`
- Macrophages: `C1QA`, `C1QB`, `C1QC`, `APOE`, `SPP1`, `MSR1`
- Dendritic cells: `FCER1A`, `CLEC10A`, `CD1C`, `CLEC9A`, `XCR1`, `LAMP3`, `LILRA4`
- Mast cells: `TPSAB1`, `TPSB2`, `KIT`, `CPA3`
- Cycling cells: `MKI67`, `TOP2A`, `STMN1`

Required figures:

- `gse189926_umap_overview.png`: UMAP colored by sample, patient, raw outcome, derived timepoint, and working cell type label.
- `gse189926_marker_dotplot.png`: dot plot for the marker groups above by working cell type label or cluster.

The figures are QC figures, not manuscript figures.

### 2. GSE235863 CD45 Response-Aware Summaries

Inputs:

- `GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad.gz`.
- Previous response join audit and cell type count tables.

Required handling:

- Preserve raw response strings and derived response grouping.
- Preserve the `P18` ambiguity flag.
- Summarize cell counts and proportions by:
  - raw response
  - derived response group
  - patient
  - sample
  - tissue
  - `major_cluster`
  - `sub_cluster`
- Export response-aware abundance tables and small overview figures.
- Do not modify or re-label existing `major_cluster` and `sub_cluster` fields.

Required figures:

- `gse235863_cd45_umap_overview.png`: UMAP colored by major cluster, raw response, derived response group, tissue, and patient.
- `gse235863_cd45_abundance_overview.png`: stacked composition overview by response and tissue.

### 3. GSE235863 CD8 Doublet Sensitivity

Inputs:

- `GSE235863_five_patients_scRNAseq_cd8t_raw_counts.h5ad.gz`.

Required handling:

- Preserve raw response strings and derived response grouping.
- Report subcluster counts before and after excluding cells where `predicted_doublets` is true.
- Do not overwrite the original object.
- Upload the exact retained and excluded cell counts.

Required figures:

- `gse235863_cd8_umap_overview.png`: UMAP colored by subcluster, raw response, derived response group, tissue, and predicted doublet status.
- `gse235863_cd8_abundance_overview.png`: subcluster composition overview by response, before and after excluding predicted doublets.

### 4. OMIX001073 Reference Support

Inputs:

- Previous OMIX001073 annotation summaries.
- OMIX001073 all-cell, myeloid, and stromal annotation tables.

Required handling:

- Create a compact reference support table for gastric myeloid and stromal labels.
- Include cluster labels, cell counts, source annotation file, and marker genes if already present in the source files.
- Do not run integration with `GSE189926` in this task.

### 5. GSE236581 Exclusion Note

Inputs:

- `gse236581_response_column_audit.tsv`.

Required handling:

- Upload `gse236581_response_exclusion_note.tsv` stating that no response-related metadata column was detected and that the dataset remains excluded from response analysis.
- Do not infer response labels.

## Scientific Guardrails

- Preserve all raw labels alongside derived labels.
- Do not remove samples or cell clusters without an explicit audit table.
- Do not downsample unless the run fails without downsampling. If downsampling is unavoidable for a QC figure, upload the retained counts and explain that the object itself remains full-size.
- Keep `GSE189926` `NE` samples separate from binary response comparisons.
- Keep `GSE236581` out of response analyses.
- Record all package versions in `session_info.txt`.

## Completion Signal

After uploading results, commit and push to the same branch:

`pan-gi-ici-metadata-request-20260709`

Commit message:

`Add Pan-GI ICI response object construction outputs`
