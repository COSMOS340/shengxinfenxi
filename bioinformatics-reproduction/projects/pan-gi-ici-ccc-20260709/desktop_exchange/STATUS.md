# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-09

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

Run response-ready object construction and cell type audit.

Do not run CellChat, LIANA, NicheNet, differential communication, manuscript figures, or final statistical testing in this task. This step should build the `GSE189926` sparse object, export QC overviews, clustering/marker tables, and cell type label audit files, then summarize `GSE235863` response-aware abundance patterns.

## Current Request Files

- `20260709_response_object_construction_manifest.tsv`
- `20260709_response_object_construction_request.md`

## Input Audits From Previous PC Steps

The matrix download audit is uploaded at:

`uploads/20260709_matrix_download`

The priority response preprocessing audit is uploaded at:

`uploads/20260709_preprocess_priority_response`

Use `matrix_file_inventory.tsv` to locate PC-local matrix files. Large matrices were intentionally not uploaded to GitHub.

## Required Upload Directory For This Step

After object construction and label audit, upload lightweight outputs to:

`uploads/20260709_response_object_construction`

Large processed objects should remain PC-local and be listed in `object_inventory.tsv`.

## Earlier Request Files

- `20260709_metadata_download_manifest.tsv`
- `20260709_metadata_download_request.md`
- `20260709_small_file_download_manifest.tsv`
- `20260709_small_file_download_request.md`
- `20260709_matrix_download_priority12.tsv`
- `20260709_matrix_download_request.md`
- `20260709_preprocess_priority_response_manifest.tsv`
- `20260709_preprocess_priority_response_request.md`
