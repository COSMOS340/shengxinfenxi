# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-11

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

No open PC-side request is pending in this exchange directory.

The most recent request, `20260711_gse235863_timepoint_prjna932556_request.md`, is complete. `GSE235863` was rebuilt without collapsing tissue or treatment timepoint. `PRJNA932556` stopped after the mandatory pilot with `blocked_read_structure_not_validated`, so no full count matrices or Seurat QC object were fabricated.

## Previous Status

PC-side requests are complete through the `GSE235863` sample-timepoint rebuild and `PRJNA932556` pilot audit.

Most recent completion commits:

- `c31eff2be250591944f4a8c6079d2b897ec1171c`: added the `GSE235863` sample-timepoint rebuild and `PRJNA932556` pilot audit outputs.
- `b89f9b811ac17ec1877c94f4db370cfe4cf00336`: refreshed manual-source audit statuses for metadata and small-file downloads.
- `2327677ec877f904c6aea0938a10936f621ff007`: added the R-only `GSE189926` rerun outputs.

## Completed Request Files

- `20260711_gse235863_timepoint_prjna932556_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_request.md`
- `20260710_gse189926_r_only_rerun_manifest.tsv`
- `20260710_gse189926_r_only_rerun_request.md`
- `20260710_gse189926_annotation_qc_refinement_manifest.tsv`
- `20260710_gse189926_annotation_qc_refinement_request.md`
- `20260710_gse189926_umap_repair_policy.md`

## Completed Upload Directories

- `uploads/20260711_gse235863_timepoint_prjna932556`
- `uploads/20260710_gse189926_r_only_rerun`
- `uploads/20260710_gse189926_annotation_qc_refinement`
- `uploads/20260709_metadata_download`
- `uploads/20260709_small_file_download`
- `uploads/20260709_matrix_download`
- `uploads/20260709_preprocess_priority_response`
- `uploads/20260709_response_object_construction`

The `GSE189926` R-only rerun used R/Seurat from the raw 22 matrices, reran QC, UMAP, sample-aware correction, clustering, marker export, and immune label audit. CellChat, LIANA, NicheNet, differential communication, manuscript figures, and final statistical testing were not run.

Manual metadata rows that could not be resolved from public reproducible sources are marked as `verified_unavailable_public_sources` with evidence tables; unknown response labels were not inferred.

Large RDS/h5ad objects remain PC-local and are recorded in the relevant `object_inventory.tsv` files.

## Input Audits From Previous PC Steps

The matrix download audit is uploaded at:

`uploads/20260709_matrix_download`

The priority response preprocessing audit is uploaded at:

`uploads/20260709_preprocess_priority_response`

The response object construction output is uploaded at:

`uploads/20260709_response_object_construction`

Use `matrix_file_inventory.tsv` in the matrix download output to locate the PC-local raw `GSE189926` matrices. Use the relevant `object_inventory.tsv` files to locate PC-local large processed objects. Large objects were intentionally not uploaded to GitHub.

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
