# PC Request: CRLM CAF/ECM Initial QC Preflight

Created: 2026-07-08

## Purpose

Use the PC-side raw GEO files to generate objective QC and metadata summaries before any filtering, integration, clustering, cell-type annotation, or figure generation.

This is a pre-analysis audit step. It should not make final biological claims.

## Existing PC Raw Input

Use the raw directory from the previous PC download task:

`I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-caf-ecm-spatial-niche\desktop_exchange\requests\20260708_core_geo_download\raw_core_geo`

The raw structure audit has already confirmed:

- `GSE225857`: non-immune counts/meta plus 6 spatial samples.
- `GSE178318`: one matrix with 33,694 genes and 140,281 cells.
- `GSE245552`: 39 sample matrices with 190,190 total cells.

## Run

```bash
Rscript bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_initial_qc_preflight/scripts/pc_initial_qc_preflight.R
```

Required R package:

- `Matrix`

Do not install packages inside the script. If `Matrix` is unavailable, install it on PC and rerun.

## Required Outputs

Commit only small tables/logs under:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/uploads/20260708_initial_qc_preflight/`

Expected outputs:

- `STATUS.md`
- `gse225857_nonimmune_meta_columns.tsv`
- `gse225857_nonimmune_meta_value_counts.tsv`
- `gse225857_nonimmune_qc_by_field.tsv`
- `gse178318_barcode_suffix_counts.tsv`
- `gse178318_qc_by_suffix.tsv`
- `gse245552_sample_qc_summary.tsv`
- `gse245552_sample_cell_counts.tsv`

Do not commit raw files, extracted matrices, Seurat objects, RDS objects, or per-cell QC tables.

## Boundaries

- Do not filter cells.
- Do not run doublet detection.
- Do not run normalization, PCA, Harmony, UMAP, or clustering.
- Do not assign cell types.
- Do not infer patient IDs from partial strings.
- Preserve exact raw sample labels, barcode suffixes, and GEO sample fields.

