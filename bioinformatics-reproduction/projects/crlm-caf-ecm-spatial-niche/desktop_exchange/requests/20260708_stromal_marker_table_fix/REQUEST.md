# CRLM CAF/ECM Stromal Marker Table Fix Request

Date: 2026-07-08

## Purpose

Repair the uploaded stromal subset all-gene marker review tables from the previous clustering step.

The previous upload successfully produced per-dataset Seurat objects, UMAP plots, marker dot plots, composition tables, marker score tables, and PC-local full marker RDS files. However, the uploaded `stromal_subset_cluster_markers_top10.tsv` and `stromal_subset_cluster_markers_top50.tsv` files did not include `dataset_id`. Because cluster identifiers are reused across `GSE178318` and `GSE245552`, those two TSV files are ambiguous for final label review.

This request only re-exports marker TSV files with explicit `dataset_id`.

## Inputs

Use the PC-local marker RDS files already written by the previous request:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_core_geo_download/derived_stromal_subset_clustering_20260708/`

Expected files:

- `gse178318_stromal_subset_all_markers_res0.4.rds`
- `gse245552_stromal_subset_all_markers_res0.4.rds`

If either file is missing, stop with an explicit error. Do not rerun Seurat clustering or marker testing for this repair step.

## Command

Run from the repository root:

```bash
Rscript bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_stromal_marker_table_fix/scripts/pc_fix_stromal_marker_tables.R
```

The script uses:

- `data.table`

The script can install the missing CRAN package on the PC if needed.

## Expected GitHub Upload Directory

Upload only small tables and status files to:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/uploads/20260708_stromal_marker_table_fix/`

Do not commit PC-local Seurat objects or full marker RDS files.

## Expected Uploaded Files

- `STATUS.md`
- `stromal_subset_cluster_marker_count_summary.tsv`
- `stromal_subset_cluster_markers_top10_by_dataset.tsv`
- `stromal_subset_cluster_markers_top50_by_dataset.tsv`
- `stromal_subset_marker_table_fix_manifest.tsv`

## Interpretation Boundary

This request repairs evidence tables only. It does not finalize cell labels, filter cells, merge datasets, draw publication figures, or change the previous clustering.
