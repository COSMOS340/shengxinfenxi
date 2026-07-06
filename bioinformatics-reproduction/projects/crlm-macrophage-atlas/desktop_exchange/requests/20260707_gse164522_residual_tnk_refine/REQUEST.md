# PC Request: GSE164522 Residual T/NK Refined Macrophage/Monocyte Recluster

## Context

The cleaned recluster run completed, but post-clean marker review still found a small residual T/NK-like population inside the macrophage/monocyte analysis:

- `broad_clean_singlets` clean cluster `13`: 105 cells; top marker evidence includes `TRAC`, `CD3D`, `CD247`, `TRBC2`, `CD3E`, `PRF1`, and `TRBC1` in the top50 marker scan.
- `core_clean_singlets` clean cluster `11`: 35 cells; top10 markers include `CD3E`, `CD3D`, `TRBC2`, and `CD247`.

This request performs a narrow refinement only. Do not remove additional clusters unless a separate request is issued.

## Input Object

Use the same PC-side Seurat RDS from the previous full GSE164522 macrophage/monocyte run:

`requests/20260706_gse164522_full_mac_mono_analysis/local_objects_not_for_github/gse164522_full_mac_mono_seurat.rds`

If needed, set:

```powershell
$env:GSE164522_MAC_MONO_RDS = "I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-macrophage-atlas\desktop_exchange\requests\20260706_gse164522_full_mac_mono_analysis\local_objects_not_for_github\gse164522_full_mac_mono_seurat.rds"
```

## Refinement Rule

The input cell lists are explicit:

- `inputs/gse164522_broad_refined_singlet_cells.tsv`
  - Source: `broad_clean_singlets`
  - Remove only post-clean `clean_cluster == 13`
  - Retained cells: 9,834
  - Removed cells: 105
- `inputs/gse164522_core_refined_singlet_cells.tsv`
  - Source: `core_clean_singlets`
  - Remove only post-clean `clean_cluster == 11`
  - Retained cells: 3,443
  - Removed cells: 35
- `inputs/gse164522_residual_tnk_excluded_cells.tsv`
  - Audit list of removed cells and their source metadata.
- `inputs/gse164522_residual_tnk_refinement_input_summary.tsv`
  - Summary of source cells, removed clusters, removed cells, and retained cells.

The 35 cells removed from `core_clean_singlets` are included within the 105 cells removed from `broad_clean_singlets`.

## Command

From this request directory:

```powershell
Rscript scripts/pc_gse164522_residual_tnk_refine.R
python scripts/pc_prepare_upload_package.py
```

## Expected Outputs

The R script writes to:

`requests/20260707_gse164522_residual_tnk_refine/outputs/`

The packaging script copies those files to:

`uploads/20260707_gse164522_residual_tnk_refine/`

Expected upload files include:

- `outputs/STATUS.md`
- `outputs/sessionInfo.txt`
- `outputs/gse164522_residual_tnk_refined_combined_summary.tsv`
- For each refined set:
  - `outputs/gse164522_<set>_run_summary.tsv`
  - `outputs/gse164522_<set>_input_cell_audit.tsv`
  - `outputs/gse164522_<set>_umap_coordinates.tsv.gz`
  - `outputs/gse164522_<set>_cluster_markers.tsv.gz`
  - `outputs/gse164522_<set>_cluster_top50_markers.tsv`
  - `outputs/gse164522_<set>_cluster_tissue_composition.tsv`
  - `outputs/gse164522_<set>_cluster_sample_composition.tsv`
  - `outputs/gse164522_<set>_cluster_patient_composition.tsv`
  - `outputs/gse164522_<set>_cluster_original_cluster_composition.tsv`
  - `outputs/gse164522_<set>_cluster_cleanup_label_composition.tsv`
  - `outputs/gse164522_<set>_cluster_celltype_sub_composition.tsv`
  - UMAP and marker PNGs.
- `upload_manifest.tsv`

Set names:

- `broad_refined_singlets`
- `core_refined_singlets`

## Minimum QC Before Upload

1. Confirm both refined sets ran to completion.
2. Confirm `outputs/gse164522_residual_tnk_refined_combined_summary.tsv` has one row per refined set.
3. Confirm the input audits report:
   - `broad_refined_singlets`: 9,834 requested cells, 9,834 present cells, 0 missing cells.
   - `core_refined_singlets`: 3,443 requested cells, 3,443 present cells, 0 missing cells.
4. Confirm every uploaded PNG is nonblank, not clipped, and has readable legends.
5. Run `python scripts/pc_prepare_upload_package.py`.
6. Commit only `bioinformatics-reproduction/projects/crlm-macrophage-atlas/desktop_exchange/uploads/20260707_gse164522_residual_tnk_refine/`.

## Do Not Commit

Do not commit raw GEO files.

Do not commit `.rds`, `.rda`, `.RData`, package caches, or temporary files.

Do not commit the request `outputs/` directory directly. Use the packaging script and commit only the upload directory.
