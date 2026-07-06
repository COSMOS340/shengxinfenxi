# PC Request: GSE164522 Macrophage/Monocyte Doublet And Batch Audit

## Context

The previous request `20260706_gse164522_full_mac_mono_analysis` produced a full `GSE164522` Macrophage/Monocyte Seurat object with 15,853 selected cells and 23 resolution 0.6 clusters.

Mac-side marker review found that several clusters have strong non-myeloid boundary signals, but direct removal of many clusters would be weak for peer review unless we also check doublets and batch or sample effects. This request is an audit only. Do not delete clusters in this job.

## Input Object

The script will first look for the local-only Seurat object produced by the previous PC job:

`requests/20260706_gse164522_full_mac_mono_analysis/local_objects_not_for_github/gse164522_full_mac_mono_seurat.rds`

If the object was moved, set:

```powershell
$env:GSE164522_MAC_MONO_RDS = "I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-macrophage-atlas\desktop_exchange\requests\20260706_gse164522_full_mac_mono_analysis\local_objects_not_for_github\gse164522_full_mac_mono_seurat.rds"
```

## Command

From this request directory:

```powershell
Rscript scripts/pc_gse164522_doublet_batch_audit.R
python scripts/pc_prepare_upload_package.py
```

The R script can install missing CRAN/Bioconductor packages by default. To disable package installation:

```powershell
$env:INSTALL_MISSING_R_PACKAGES = "0"
```

## Scientific Checks

The audit must separate four explanations:

- Non-myeloid lineage boundary signal from top markers and module scores.
- Doublet enrichment from `scDblFinder`.
- Sample, patient, tissue, or author subtype dominance.
- QC complexity differences, including `n_genes`, `n_counts`, `percent_mito`, `nFeature_RNA`, `nCount_RNA`, and calculated mitochondrial percentage when available.

The output should support a conservative decision table. It must not output a final deletion set.

## Expected Outputs

The R script writes to:

`requests/20260706_gse164522_doublet_batch_audit/outputs/`

The packaging script copies those files to:

`uploads/20260706_gse164522_doublet_batch_audit/`

Expected upload files include:

- `outputs/STATUS.md`
- `outputs/sessionInfo.txt`
- `outputs/gse164522_audit_marker_set_coverage.tsv`
- `outputs/gse164522_audit_per_cell_lineage_scores.tsv.gz`
- `outputs/gse164522_audit_cluster_lineage_score_summary.tsv`
- `outputs/gse164522_audit_cluster_qc_summary.tsv`
- `outputs/gse164522_audit_cluster_composition_summary.tsv`
- `outputs/gse164522_audit_cluster_knn_batch_mixing.tsv` if nearest-neighbor mixing succeeds.
- `outputs/gse164522_audit_scDblFinder_per_cell.tsv.gz` if `scDblFinder` succeeds.
- `outputs/gse164522_audit_scDblFinder_cluster_summary.tsv` if `scDblFinder` succeeds.
- `outputs/gse164522_audit_harmony_sample_umap_coordinates.tsv.gz` and Harmony UMAP PNGs if `harmony` succeeds.
- `outputs/gse164522_audit_cluster_integrated_review.tsv`
- `outputs/gse164522_audit_cluster_review_status.tsv`
- PNG figures with UMAP, QC, dominance, doublet, and lineage score summaries.
- `upload_manifest.tsv`

## Minimum QC Before Upload

1. Confirm `outputs/STATUS.md` exists and has normal line breaks.
2. Confirm `gse164522_audit_cluster_integrated_review.tsv` has all 23 clusters.
3. If `scDblFinder` fails, keep the error message in `STATUS.md` and still upload all other audit outputs.
4. Open every PNG locally and check that the plot is not blank, labels are not clipped, and legends do not overflow.
5. Run `python scripts/pc_prepare_upload_package.py` after the R script finishes.
6. Commit only `bioinformatics-reproduction/projects/crlm-macrophage-atlas/desktop_exchange/uploads/20260706_gse164522_doublet_batch_audit/`.

## Do Not Commit

Do not commit raw GEO files.

Do not commit any `.rds`, `.rda`, `.RData`, package cache, or temporary file.

Do not commit `outputs/` from the request directory. Use the packaging script and commit only the upload directory.
