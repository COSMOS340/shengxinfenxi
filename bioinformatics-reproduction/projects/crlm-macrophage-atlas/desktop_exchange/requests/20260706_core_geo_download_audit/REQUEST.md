# PC request: CRLM core GEO download and audit

Created: 2026-07-06

## Purpose

Download and audit the core public GEO files for the CRLM macrophage atlas project. The Mac side can download some files, but large NCBI GEO transfers are currently slow and unstable. Do not upload raw matrices to GitHub. Upload only small audit tables and logs.

## Core accessions

- `GSE164522`
- `GSE178318`
- `GSE225857`

## Required actions

1. Run `scripts/pc_download_and_audit_core_geo.py` from this request directory.
2. Keep downloaded raw files on the PC.
3. Upload the generated `outputs/` directory to a new upload folder.
4. Do not commit raw `.gz`, `.mtx`, `.tar`, image, or matrix files.

## Required return files

- `outputs/core_geo_download_status.tsv`
- `outputs/core_geo_file_integrity.tsv`
- `outputs/gse164522_expression_dimensions.tsv`
- `outputs/gse178318_matrix_dimensions.tsv`
- `outputs/gse225857_tar_filelist.tsv`
- `outputs/STATUS.md`

## Notes for interpretation

- `GSE164522` metadata cell IDs use a terminal hyphen-number suffix, while expression matrix columns use a dot-number suffix. The Mac audit confirmed this mapping for `MN` and `MT`.
- `GSE178318` has data usage terms on the GEO page. Record that note in `STATUS.md`.
- `GSE225857_RAW.tar` contains processed scRNA and spatial files. The tar file list is needed before deciding which files to parse.
