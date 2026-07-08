# PC Request: CRLM CAF/ECM Working QC Sparse Objects

Created: 2026-07-08

## Purpose

Use the PC-side raw GEO files to create traceable working QC sparse objects and small review tables for the CRLM CAF/ECM project.

This request should prepare the next analysis step. It should not perform final filtering, integration, clustering, cell-type annotation, or publication figure generation.

## Existing PC Raw Input

Use the raw directory from the previous PC download task:

`I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-caf-ecm-spatial-niche\desktop_exchange\requests\20260708_core_geo_download\raw_core_geo`

Prior QC decisions for this working pass:

- `GSE225857`: use author non-immune metadata and exact source cluster labels; do not re-filter cells in this request.
- `GSE178318`: working sensitivity checkpoint is `nFeature_RNA >= 500` and `percent.mt <= 15`.
- `GSE245552`: working sensitivity checkpoint is `nFeature_RNA >= 500` and `nCount_RNA >= 1000`; do not use mitochondrial percentage because feature files differ between tumor/metastasis and adjacent samples.

## Run

```bash
Rscript bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_working_qc_sparse_objects/scripts/pc_working_qc_sparse_objects.R
```

Required R package:

- `Matrix`

Do not install packages inside the script. If `Matrix` is unavailable, install it on PC and rerun.

## Expected Local-Only Object Outputs

The script writes RDS objects outside the upload folder under:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_core_geo_download/derived_working_qc_sparse_objects_20260708/`

These RDS files must stay PC-local. Do not commit them.

## Required Upload Outputs

Commit only small tables/logs/PNG files under:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/uploads/20260708_working_qc_sparse_objects/`

Expected outputs:

- `STATUS.md`
- `working_qc_thresholds.tsv`
- `working_qc_object_manifest.tsv`
- `gse225857_nonimmune_source_cluster_summary.tsv`
- `gse225857_nonimmune_fibroblast_source_clusters.tsv`
- `gse178318_working_qc_retention_by_suffix.tsv`
- `gse245552_working_qc_retention_by_sample.tsv`
- `caf_ecm_marker_gene_presence.tsv`
- `gse178318_working_qc_retention_by_suffix.png`
- `gse245552_working_qc_retention_by_sample.png`

## Boundaries

- Do not commit raw files or RDS objects.
- Do not run integration, PCA, Harmony, UMAP, clustering, or marker testing.
- Do not assign final cell types.
- Do not infer patient IDs beyond exact raw labels and GEO metadata.
- Preserve exact source cluster strings from `GSE225857`.
