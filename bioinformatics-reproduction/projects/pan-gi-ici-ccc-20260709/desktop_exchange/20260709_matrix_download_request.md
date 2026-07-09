# PC Task: Matrix Download and Lightweight Inspection

Created: 2026-07-09

## Goal

Download priority 1-2 matrix files for the first formal pan-GI ICI communication analysis set.

This task is download and lightweight inspection only. Do not run normalization, clustering, integration, CellChat, LIANA, or plotting in this task.

## Manifest

Use this manifest:

`/Volumes/research/pan_gi_ici_ccc_20260709/99_logs/pc_tasks/20260709_matrix_download_priority12.tsv`

It contains:

- priority 1: `GSE189926` 22 gastric cancer response matrices and `GSE235863` CD45 h5ad.
- priority 2: `GSE235863` CD8 h5ad, `GSE236581` metadata/count matrix set, `OMIX001073` all-cell/myeloid/stromal matrices, and the integrated Liver Ma all-cell RDS.

## Output Root

Download files to the exact `suggested_output_path` in the manifest:

`/Volumes/research/pan_gi_ici_ccc_20260709/03_matrices`

## Upload Back Through GitHub

Do not upload large matrix files to GitHub.

Upload only lightweight audit files to:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_matrix_download`

Required upload files:

- `download_status.tsv`
- `SHA256SUMS.txt`
- `matrix_file_inventory.tsv`
- `lightweight_inspection.tsv`
- `STATUS.md`

## Required Columns

`download_status.tsv`:

- priority
- dataset_id
- file_name
- status
- output_path
- file_size_bytes
- sha256
- error_message

`matrix_file_inventory.tsv`:

- dataset_id
- file_name
- local_path
- file_size_bytes
- sha256
- upload_status

Use `local_only_large_file` for files not uploaded to GitHub because they are large.

`lightweight_inspection.tsv`:

- dataset_id
- file_name
- inspection_status
- detected_format
- n_obs_or_cells
- n_vars_or_genes
- obs_columns_or_notes
- error_message

## Lightweight Inspection Guidance

For text matrices:

- report first few lines and inferred delimiter in `obs_columns_or_notes`.
- for `GSE189926` matrix text files, check whether rows appear to be genes and columns are cells.

For `GSE236581_CRC-ICB_metadata.txt.gz`:

- inspect column names.
- report whether any response-related field exists.

For h5ad:

- if Python can install packages, use `anndata` to read `.h5ad.gz` after decompressing or with a temporary file.
- report `n_obs`, `n_vars`, and `obs` column names.
- do not export full expression matrices.

For RDS:

- if R can install packages, use `readRDS()` only to report object class, dimensions, metadata column names, and unique disease/cancer labels if easy.
- do not save a converted full object.

## Scientific Priority

The most important outputs are:

1. `GSE189926` matrices, because response metadata is already directly available.
2. `GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad.gz`, because response metadata is directly available.
3. `GSE236581_CRC-ICB_metadata.txt.gz`, because it may resolve CRC response metadata.
4. `OMIX001073-20-14.zip` and `OMIX001073-20-12.zip`, because they support gastric myeloid and stromal reference states.

If time or storage is limited, complete priority 1 first and record skipped priority 2 rows as `skipped_storage_or_time`.
