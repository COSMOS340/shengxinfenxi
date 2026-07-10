# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-10

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

Run `GSE189926` annotation and QC refinement.

Do not run CellChat, LIANA, NicheNet, differential communication, manuscript figures, or final statistical testing in this task. This step should refine `GSE189926` immune labels, audit QC flags, quantify sample/patient structure, redraw readable QC figures, and produce patient-level `GSE235863` composition summaries.

## Current Request Files

- `20260710_gse189926_annotation_qc_refinement_manifest.tsv`
- `20260710_gse189926_annotation_qc_refinement_request.md`

## Input Audits From Previous PC Steps

The matrix download audit is uploaded at:

`uploads/20260709_matrix_download`

The priority response preprocessing audit is uploaded at:

`uploads/20260709_preprocess_priority_response`

The response object construction output is uploaded at:

`uploads/20260709_response_object_construction`

Use `object_inventory.tsv` in the response object construction output to locate the PC-local `GSE189926` h5ad object. Large objects were intentionally not uploaded to GitHub.

## Required Upload Directory For This Step

After annotation and QC refinement, upload lightweight outputs to:

`uploads/20260710_gse189926_annotation_qc_refinement`

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
- `20260709_response_object_construction_manifest.tsv`
- `20260709_response_object_construction_request.md`
