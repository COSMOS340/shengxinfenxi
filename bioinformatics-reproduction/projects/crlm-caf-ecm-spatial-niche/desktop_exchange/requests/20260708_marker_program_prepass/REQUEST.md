# PC Request: CRLM CAF/ECM Marker Program Prepass

Created: 2026-07-08

## Purpose

Use the PC-local working QC sparse RDS objects to generate broad marker-program review summaries for the CRLM CAF/ECM project.

This is a review prepass. It should not perform final cell-type annotation, integration, clustering, marker testing, or publication figure generation.

## Existing PC Inputs

Working QC RDS objects from the previous request are PC-local under:

`I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-caf-ecm-spatial-niche\desktop_exchange\requests\20260708_core_geo_download\derived_working_qc_sparse_objects_20260708`

Raw input directory:

`I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-caf-ecm-spatial-niche\desktop_exchange\requests\20260708_core_geo_download\raw_core_geo`

Important issue to fix:

`gse225857_nonimmune_source_cluster_summary.tsv` includes `F06_cycling_MKI67` with 227 cells, but the detailed fibroblast table from the previous request omitted it because that script matched `fibroblast/fibrblast` in the cluster name. This request must regenerate the `GSE225857` fibroblast-prefix table using the exact `F` source prefix.

## Run

```bash
Rscript bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_marker_program_prepass/scripts/pc_marker_program_prepass.R
```

Required R package:

- `Matrix`

Do not install packages inside the script. If `Matrix` is unavailable, install it on PC and rerun.

## Required Upload Outputs

Commit only small tables/logs/PNG files under:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/uploads/20260708_marker_program_prepass/`

Expected outputs:

- `STATUS.md`
- `gse225857_f_prefix_source_clusters_by_organ.tsv`
- `marker_program_gene_coverage.tsv`
- `marker_program_score_summary_by_group.tsv`
- `marker_program_top_counts_by_group.tsv`
- `stromal_review_counts_by_group.tsv`
- `marker_program_prepass_object_manifest.tsv`
- `gse178318_marker_program_top_counts_by_suffix.png`
- `gse245552_marker_program_top_counts_by_tissue.png`

Do not commit working RDS objects or per-cell score tables.

## Boundaries

- Do not assign final cell types.
- Do not run integration, PCA, Harmony, UMAP, clustering, or marker testing.
- Do not upload per-cell score tables.
- Do not remove samples or cells in this request.
- Treat marker-program summaries as review evidence only.
