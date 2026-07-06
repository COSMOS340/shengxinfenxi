# PC Request: GSE164522 Full Macrophage/Monocyte Analysis

## Context

The previous PC job completed the core GEO download audit for the CRLM macrophage atlas project. The full raw files are already present on PC under the prior request directory:

`bioinformatics-reproduction/projects/crlm-macrophage-atlas/desktop_exchange/requests/20260706_core_geo_download_audit/raw_core_geo`

This request runs the first full all-tissue `GSE164522` analysis needed for Fig.1. The local Mac only has complete MN/MT expression files, so the all-tissue analysis should run on PC.

## Scientific Scope

- Dataset: `GSE164522`.
- Expression files: exact six files in the previous raw directory:
  - `GSE164522_CRLM_LN_expression.csv.gz`
  - `GSE164522_CRLM_MN_expression.csv.gz`
  - `GSE164522_CRLM_MT_expression.csv.gz`
  - `GSE164522_CRLM_PBMC_expression.csv.gz`
  - `GSE164522_CRLM_PN_expression.csv.gz`
  - `GSE164522_CRLM_PT_expression.csv.gz`
- Metadata file: `GSE164522_CRLM_metadata.csv.gz`.
- Primary cell selection: exact `celltype_major` values `Macrophage` and `Monocyte`.
- Do not mix `DC` or `MAST` into this primary run. Those labels can be checked later as boundary controls.

## Command

From this request directory:

```powershell
Rscript scripts/pc_gse164522_full_mac_mono_analysis.R
python scripts/pc_prepare_upload_package.py
```

If the script cannot find the previous raw directory automatically, set the raw path explicitly first:

```powershell
$env:RAW_CORE_GEO_DIR = "I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-macrophage-atlas\desktop_exchange\requests\20260706_core_geo_download_audit\raw_core_geo"
Rscript scripts/pc_gse164522_full_mac_mono_analysis.R
python scripts/pc_prepare_upload_package.py
```

## Expected Outputs

The R script writes analysis outputs to:

`requests/20260706_gse164522_full_mac_mono_analysis/outputs/`

The packaging script copies only those outputs to:

`uploads/20260706_gse164522_full_mac_mono_analysis/`

Expected upload files include:

- `outputs/STATUS.md`
- `outputs/sessionInfo.txt`
- `outputs/gse164522_full_mac_mono_run_summary.tsv`
- `outputs/gse164522_metadata_celltype_counts.tsv`
- `outputs/gse164522_full_mac_mono_selection_audit.tsv`
- `outputs/gse164522_full_mac_mono_gene_audit.tsv`
- `outputs/gse164522_full_mac_mono_selected_cells.tsv`
- `outputs/gse164522_full_mac_mono_counts_by_tissue_patient.tsv`
- `outputs/gse164522_full_mac_mono_umap_coordinates.tsv.gz`
- `outputs/gse164522_full_mac_mono_cluster_markers.tsv.gz`
- `outputs/gse164522_full_mac_mono_cluster_top50_markers.tsv`
- `outputs/gse164522_full_mac_mono_marker_panel_coverage.tsv`
- `outputs/gse164522_full_mac_mono_cluster_tissue_composition.tsv`
- `outputs/fig1_gse164522_full_mac_mono_umap_layout.png`
- `outputs/fig1_gse164522_full_mac_mono_umap_clusters.png`
- `outputs/fig1_gse164522_full_mac_mono_umap_tissue.png`
- `outputs/fig1_gse164522_full_mac_mono_umap_author_major.png`
- `outputs/fig1_gse164522_full_mac_mono_marker_dotplot.png`
- `outputs/fig1_gse164522_full_mac_mono_cluster_tissue_composition.png`
- `outputs/fig1_gse164522_full_mac_mono_top_marker_heatmap.png` if heatmap generation succeeds.
- `upload_manifest.tsv`

## Do Not Commit

Do not commit any raw GEO files.

Do not commit `local_objects_not_for_github/`. The R script saves the Seurat RDS there for PC-side reuse only.

Do not commit partial files, package caches, or temporary R session files.

## Minimum QC Before Upload

1. Confirm `outputs/STATUS.md` exists and has normal line breaks.
2. Confirm `outputs/gse164522_full_mac_mono_selection_audit.tsv` has six expression groups.
3. Confirm the selected cell count in `outputs/gse164522_full_mac_mono_run_summary.tsv` is not zero.
4. Open every PNG locally and check that the plot is not blank, labels are not clipped, and legends do not overflow.
5. Run `python scripts/pc_prepare_upload_package.py` after the R script finishes.
6. Commit only `bioinformatics-reproduction/projects/crlm-macrophage-atlas/desktop_exchange/uploads/20260706_gse164522_full_mac_mono_analysis/`.

## If The Full Seurat Step Fails

Still upload any completed audit tables from `outputs/`, especially:

- `gse164522_metadata_celltype_counts.tsv`
- `gse164522_full_mac_mono_selection_audit.tsv`
- `gse164522_full_mac_mono_gene_audit.tsv`
- `gse164522_full_mac_mono_selected_cells.tsv` if created
- `STATUS.md` or the terminal error log

The selection audit is useful even if clustering needs to be rerun with a smaller memory strategy.
