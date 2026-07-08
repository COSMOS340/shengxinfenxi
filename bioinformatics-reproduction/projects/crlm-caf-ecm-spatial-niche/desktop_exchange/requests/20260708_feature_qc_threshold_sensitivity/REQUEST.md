# PC Request: CRLM CAF/ECM Feature QC and Threshold Sensitivity

Created: 2026-07-08

## Purpose

Use the PC-side raw GEO files to resolve the GSE245552 feature-column issue and generate objective QC threshold sensitivity tables before any final filtering or integration.

This request is still pre-analysis. It should not make final biological claims.

## Existing PC Raw Input

Use the raw directory from the previous PC download task:

`I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-caf-ecm-spatial-niche\desktop_exchange\requests\20260708_core_geo_download\raw_core_geo`

The initial QC preflight found:

- `GSE225857`: 41,892 author non-immune metadata cells.
- `GSE178318`: 140,281 matrix cells across 15 barcode suffix groups.
- `GSE245552`: 39 matrices and 190,190 total cells.

Important issue to resolve: the first GSE245552 preflight produced mitochondrial percentages that were mostly zero in tumor/metastasis samples but nonzero in several adjacent samples. Do not interpret this as biology. Inspect the exact feature files first.

## Run

```bash
Rscript bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_feature_qc_threshold_sensitivity/scripts/pc_feature_qc_threshold_sensitivity.R
```

Required R package:

- `Matrix`

Do not install packages inside the script. If `Matrix` is unavailable, install it on PC and rerun.

## Required Outputs

Commit only small tables/logs under:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/uploads/20260708_feature_qc_threshold_sensitivity/`

Expected outputs:

- `STATUS.md`
- `gse245552_feature_column_audit.tsv`
- `gse245552_mito_qc_by_feature_column.tsv`
- `gse178318_threshold_retention_by_suffix.tsv`
- `gse245552_threshold_retention_by_sample.tsv`

Do not commit raw files, extracted matrices, Seurat objects, RDS objects, or per-cell QC tables.

## Boundaries

- Do not filter cells.
- Do not run doublet detection.
- Do not run normalization, PCA, Harmony, UMAP, or clustering.
- Do not assign cell types.
- Do not infer patient IDs from partial strings.
- Preserve exact raw sample labels, barcode suffixes, and GEO sample fields.
- For GSE245552, report every feature column exactly. Do not silently choose a gene-symbol column.
