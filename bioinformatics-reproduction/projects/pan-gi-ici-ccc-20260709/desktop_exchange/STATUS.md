# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-10

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

Run the `GSE189926` task again using R only.

Do not continue the Scanpy-based repair as the current task. Use R/Seurat to rebuild `GSE189926` from the raw 22 matrices, rerun QC, UMAP, sample-aware correction, clustering, marker export, and immune label audit. Do not run CellChat, LIANA, NicheNet, differential communication, manuscript figures, or final statistical testing in this task.

## Current Request Files

- `20260710_gse189926_r_only_rerun_manifest.tsv`
- `20260710_gse189926_r_only_rerun_request.md`

## Input Audits From Previous PC Steps

The matrix download audit is uploaded at:

`uploads/20260709_matrix_download`

The priority response preprocessing audit is uploaded at:

`uploads/20260709_preprocess_priority_response`

The response object construction output is uploaded at:

`uploads/20260709_response_object_construction`

Use `matrix_file_inventory.tsv` in the matrix download output to locate the PC-local raw `GSE189926` matrices. The previous h5ad object may be used only for comparison notes, not as the analysis source for the R-only rerun.

## Required Upload Directory For This Step

After the R-only rerun, upload lightweight outputs to:

`uploads/20260710_gse189926_r_only_rerun`

Large Seurat RDS objects should remain PC-local and be listed in `object_inventory.tsv`.

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
- `20260710_gse189926_annotation_qc_refinement_manifest.tsv`
- `20260710_gse189926_annotation_qc_refinement_request.md`
- `20260710_gse189926_umap_repair_policy.md`
