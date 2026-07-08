# CRLM CAF/ECM Stromal Subset Clustering Request

Date: 2026-07-08

## Purpose

Run review-grade clustering and all-gene marker testing for stromal-supported cells.

This request uses the prior stromal first-pass review to subset cells with these review classes:

- `caf_ecm_review`
- `pericyte_review`
- `endothelial_review`

The output should support cell-type interpretation and contamination review. It is still not final publication annotation.

## Inputs

Use existing PC-local inputs:

- Working QC sparse objects:
  - `bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_core_geo_download/derived_working_qc_sparse_objects_20260708/`
- Stromal first-pass review RDS objects:
  - `bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_core_geo_download/derived_stromal_first_pass_review_20260708/`

The review RDS objects are expected to contain:

- `dataset_id`
- `object_label`
- `review_table`

The `review_table` is expected to contain:

- `review_class`
- `top_marker_program`
- `score_margin`

If any required field is absent, stop with an explicit error.

## Command

Run from the repository root:

```bash
Rscript bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_stromal_subset_clustering/scripts/pc_stromal_subset_clustering.R
```

The script uses these R packages:

- `Matrix`
- `Seurat`
- `data.table`
- `ggplot2`
- `patchwork`

The script can install missing CRAN packages on the PC if needed.

## Analysis Scope

Analyze `GSE178318` and `GSE245552` separately.

For each dataset:

1. Extract cells with `review_class` in `caf_ecm_review`, `pericyte_review`, or `endothelial_review`.
2. Build a Seurat object from the PC-local working sparse counts.
3. Normalize data.
4. Select variable features.
5. Run PCA.
6. Run UMAP.
7. Run graph clustering at resolutions 0.2, 0.4, and 0.6.
8. Use resolution 0.4 as the main review clustering.
9. Run all-gene marker testing for resolution 0.4 clusters.
10. Summarize cluster composition by review class, sample, tissue, and top marker program.
11. Produce UMAP QC plots and marker dot plots.

## Expected GitHub Upload Directory

Upload only small tables, images, and status files to:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/uploads/20260708_stromal_subset_clustering/`

Do not commit PC-local Seurat RDS objects or all-cell metadata tables.

## Expected Uploaded Files

- `STATUS.md`
- `stromal_subset_run_manifest.tsv`
- `stromal_subset_dataset_summary.tsv`
- `stromal_subset_cluster_counts.tsv`
- `stromal_subset_cluster_composition_by_review_class.tsv`
- `stromal_subset_cluster_composition_by_sample.tsv`
- `stromal_subset_cluster_composition_by_tissue.tsv`
- `stromal_subset_cluster_composition_by_top_marker_program.tsv`
- `stromal_subset_cluster_markers_top50.tsv`
- `stromal_subset_cluster_markers_top10.tsv`
- `stromal_subset_marker_score_by_cluster.tsv`
- `gse178318_stromal_subset_umap_by_cluster.png`
- `gse178318_stromal_subset_umap_by_review_class.png`
- `gse178318_stromal_subset_marker_dotplot.png`
- `gse245552_stromal_subset_umap_by_cluster.png`
- `gse245552_stromal_subset_umap_by_review_class.png`
- `gse245552_stromal_subset_umap_by_tissue.png`
- `gse245552_stromal_subset_marker_dotplot.png`

## PC-Local Outputs

Write PC-local RDS objects under:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_core_geo_download/derived_stromal_subset_clustering_20260708/`

These RDS objects can include Seurat objects, full marker tables, and per-cell metadata. They should not be committed to GitHub.

## Interpretation Boundary

This request creates evidence for label review. It does not finalize cell labels.

Final labels should be assigned only after reviewing:

- all-gene cluster markers
- marker score patterns
- review-class composition
- sample and tissue distribution
- epithelial marker carryover
- immune marker carryover
- doublet-risk patterns

