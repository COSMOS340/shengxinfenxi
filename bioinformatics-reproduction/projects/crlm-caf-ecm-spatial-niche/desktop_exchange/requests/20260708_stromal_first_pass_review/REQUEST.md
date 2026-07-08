# CRLM CAF/ECM Stromal First-Pass Review Request

Date: 2026-07-08

## Purpose

Run a first-pass cell-level review of CAF/ECM-supported cells using the PC-local working QC sparse objects and the PC-local marker-program prepass score objects.

This request is a review and evidence-generation step. It must not be treated as final cell-type annotation.

## Inputs

Use existing PC-local inputs from the prior requests:

- Working QC sparse objects:
  - `bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_core_geo_download/derived_working_qc_sparse_objects_20260708/`
- Marker-program prepass score objects:
  - `bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_core_geo_download/derived_marker_program_prepass_20260708/`

The score RDS objects are expected to contain:

- `dataset_id`
- `object_label`
- `sample_label`
- `tissue_site`
- `score_matrix`
- `top_marker_program`

## Command

Run from the repository root:

```bash
Rscript bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_stromal_first_pass_review/scripts/pc_stromal_first_pass_review.R
```

The script only requires the R package `Matrix`.

## Review Classes

The script maps broad marker-program winners into review classes:

- `caf_ecm_review`: `fibroblast_ecm`, `caf_activation`, `inflammatory_caf`
- `pericyte_review`: `pericyte_like`
- `endothelial_review`: `endothelial`
- `epithelial_review`: `epithelial_tumor`
- `immune_review`: `immune_context`

These are review labels only. They are not final annotations.

## Expected GitHub Upload Directory

Upload only small tables, images, and status files to:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/uploads/20260708_stromal_first_pass_review/`

Do not commit PC-local RDS objects.

## Expected Uploaded Files

- `STATUS.md`
- `review_class_counts_by_group.tsv`
- `review_class_score_margin_summary.tsv`
- `sample_caf_ecm_review_rank.tsv`
- `marker_expression_by_review_class.tsv`
- `dataset_review_class_marker_expression.tsv`
- `stromal_first_pass_object_manifest.tsv`
- `gse178318_review_class_fraction_by_suffix.png`
- `gse245552_review_class_fraction_by_tissue.png`
- `gse245552_review_class_fraction_by_sample.png`
- `dataset_review_class_marker_expression_heatmap.png`

## PC-Local Outputs

Write PC-local RDS objects under:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_core_geo_download/derived_stromal_first_pass_review_20260708/`

These RDS files may contain per-cell review labels, scores, and margins for downstream use. They should not be committed to GitHub.

## Constraints

- Do not run final annotation.
- Do not run integration, clustering, UMAP, marker testing, or publication figure generation in this request.
- Do not remove cells based on this review step.
- Do not upload per-cell score tables or barcode-level tables to GitHub.
- Use exact object fields from the input RDS files. If a required field is absent, stop with an explicit error.

