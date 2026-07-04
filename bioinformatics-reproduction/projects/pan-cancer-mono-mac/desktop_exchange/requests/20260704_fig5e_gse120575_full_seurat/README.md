# 20260704 Fig5E GSE120575 Full Seurat Request

## Purpose

Run a PC-side full Seurat reproduction for Fig. 5E using GSE120575 melanoma single-cell TPM data.

The Mac current-scope marker/reference-correlation version does not reproduce the original posttreatment statistics. This request asks the PC to run the heavier paper-style workflow:

1. Download GSE120575 TPM and cell metadata from GEO.
2. Convert the dense TPM text file to sparse MatrixMarket by streaming the gzip file.
3. Build a Seurat object from TPM values without an extra normalization step.
4. Cluster all cells at resolution 1.0.
5. Subset clusters `8`, `12`, and `15`, matching the paper's stated CD68/LYZ myeloid clusters.
6. Recluster that subset at resolution 0.4.
7. Annotate myeloid subclusters by correlation to the supplied Fig. 2 monocyte/macrophage reference average-expression table.
8. Compute sample-level proportions for `C1QC+ TAMs` and `THBS1+ MDSCs + SPP1+ TAMs`, split by pretreatment/posttreatment and responder/nonresponder.
9. Upload tables, plots, logs, and the final Seurat object if feasible.

## Run

From this request directory:

```bash
Rscript scripts/run_fig5e_gse120575_full_seurat.R
```

The R script calls the Python sparse converter automatically if the MatrixMarket files are not present.

## Inputs Supplied

- `inputs/fig2_monocyte_macrophage_cluster_average_expression_common_features.tsv.gz`
- `inputs/fig2_monocyte_macrophage_cluster_subtype_annotation_current_scope.tsv`
- `inputs/fig5_signature_gene_sets_unique_current_scope.tsv`

## Expected Outputs

Upload the generated `outputs/` directory under:

`desktop_exchange/uploads/20260704_fig5e_gse120575_full_seurat/`

Important files:

- `outputs/tables/fig5e_full_cluster_marker_audit.tsv`
- `outputs/tables/fig5e_full_myeloid_subcluster_annotation.tsv`
- `outputs/tables/fig5e_full_sample_proportions.tsv`
- `outputs/tables/fig5e_full_wilcoxon.tsv`
- `outputs/figures/fig5e_full_seurat_proportions.png`
- `outputs/logs/run_fig5e_gse120575_full_seurat.log`
- `outputs/files/gse120575_myeloid_reclustered_seurat.rds` if size is manageable

## Notes for PC

- Do not install or upload R packages into the repository.
- If cluster IDs `8`, `12`, and `15` do not exist after resolution 1.0, stop and upload the cluster marker audit instead of substituting other cluster IDs silently.
- If those clusters exist but do not show CD68/LYZ enrichment, still upload the audit and do not force a figure.
