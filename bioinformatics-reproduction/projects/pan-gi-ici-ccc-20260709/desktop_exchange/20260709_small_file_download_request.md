# PC Task: Small Metadata and Annotation Download

Created: 2026-07-09

## Why PC Is Needed

This Mac again failed at DNS resolution when trying to download the Zenodo small Excel index. GitHub upload/download works, but direct external metadata downloads remain unreliable.

## Goal

Download only small metadata, annotation, and response-mapping files needed to lock the first formal analysis set.

Do not download:

- large RDS files from Zenodo
- OMIX matrix zip files
- FASTQ files
- BAM or CRAM files
- h5ad objects
- spatial image files

## Manifest

`/Volumes/research/pan_gi_ici_ccc_20260709/99_logs/pc_tasks/20260709_small_file_download_manifest.tsv`

## Required Output Folders

Downloaded small files:

`/Volumes/research/pan_gi_ici_ccc_20260709/02_input`

PC upload folder in GitHub:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_small_file_download`

## Required Completion Files

Create:

- `download_status.tsv`
- `SHA256SUMS.txt`
- `upload_manifest.tsv`
- `STATUS.md`

Required columns for `download_status.tsv`:

- resource_id
- status
- output_path
- file_size_bytes
- sha256
- error_message

For `manual` rows, use one of these statuses:

- `manual_done`
- `manual_blocked`
- `manual_not_found`

## Manual Items

For manual rows, the key need is not a full paper download. We need exact sample-level mapping tables:

- `STAD-PRJEB25780`: sample-level ICI response metadata from CIDE.
- `GSE205506`: subject or sample to pCR/non-pCR response mapping.
- `GSE236581`: patient or sample to response mapping.
- `PRJNA932556`: mapping between `S1-S3`/`R1-R3` and sensitive/resistant groups, plus any processed matrix links.

If a mapping is only visible in PDF supplementary material, export it as TSV and include the source note in `STATUS.md`.

## Important

This task should keep the next step lightweight. The goal is to decide exactly which datasets enter matrix-level preprocessing, not to run preprocessing yet.
