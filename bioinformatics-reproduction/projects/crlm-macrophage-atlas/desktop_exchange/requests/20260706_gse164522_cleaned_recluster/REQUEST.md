# PC Request: GSE164522 Cleaned Macrophage/Monocyte Recluster

## Context

The doublet and batch audit for `GSE164522` completed. We should not directly delete many clusters without a reproducible cleanup step. This request rebuilds cleaned macrophage/monocyte objects from the PC-side Seurat RDS and two explicit cell lists.

## Input Object

The script first looks for the previous local-only Seurat object:

`requests/20260706_gse164522_full_mac_mono_analysis/local_objects_not_for_github/gse164522_full_mac_mono_seurat.rds`

If needed, set:

```powershell
$env:GSE164522_MAC_MONO_RDS = "I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-macrophage-atlas\desktop_exchange\requests\20260706_gse164522_full_mac_mono_analysis\local_objects_not_for_github\gse164522_full_mac_mono_seurat.rds"
```

## Cleanup Sets

- `broad_clean_singlets`: remove strong lineage-boundary clusters, hold out mixed-boundary clusters, and remove per-cell `scDblFinder.class == doublet` from retained clusters. This is the next main-analysis preview.
- `core_clean_singlets`: retain only the most stable singlet clusters from the audit. This is a robustness control, not the primary figure by itself.

The input cell lists were generated from:

- `gse164522_audit_local_cleanup_interpretation.tsv`
- `gse164522_audit_scDblFinder_per_cell.tsv.gz`

## Command

From this request directory:

```powershell
Rscript scripts/pc_gse164522_cleaned_recluster.R
python scripts/pc_prepare_upload_package.py
```

## Expected Outputs

The R script writes to:

`requests/20260706_gse164522_cleaned_recluster/outputs/`

The packaging script copies those files to:

`uploads/20260706_gse164522_cleaned_recluster/`

Expected upload files include:

- `outputs/STATUS.md`
- `outputs/sessionInfo.txt`
- `outputs/gse164522_cleaned_recluster_combined_summary.tsv`
- For each cleanup set:
  - `outputs/gse164522_<set>_run_summary.tsv`
  - `outputs/gse164522_<set>_input_cell_audit.tsv`
  - `outputs/gse164522_<set>_umap_coordinates.tsv.gz`
  - `outputs/gse164522_<set>_cluster_markers.tsv.gz`
  - `outputs/gse164522_<set>_cluster_top50_markers.tsv`
  - `outputs/gse164522_<set>_cluster_tissue_composition.tsv`
  - `outputs/gse164522_<set>_cluster_sample_composition.tsv`
  - `outputs/gse164522_<set>_cluster_patient_composition.tsv`
  - `outputs/gse164522_<set>_cluster_original_cluster_composition.tsv`
  - UMAP and marker PNGs.
- `upload_manifest.tsv`

## Minimum QC Before Upload

1. Confirm both `broad_clean_singlets` and `core_clean_singlets` ran to completion.
2. Confirm `outputs/gse164522_cleaned_recluster_combined_summary.tsv` has one row per cleanup set.
3. Confirm every uploaded PNG is nonblank, not clipped, and has readable legends.
4. Run `python scripts/pc_prepare_upload_package.py`.
5. Commit only `bioinformatics-reproduction/projects/crlm-macrophage-atlas/desktop_exchange/uploads/20260706_gse164522_cleaned_recluster/`.

## Do Not Commit

Do not commit raw GEO files.

Do not commit `.rds`, `.rda`, `.RData`, package caches, or temporary files.

Do not commit the request `outputs/` directory directly. Use the packaging script and commit only the upload directory.
