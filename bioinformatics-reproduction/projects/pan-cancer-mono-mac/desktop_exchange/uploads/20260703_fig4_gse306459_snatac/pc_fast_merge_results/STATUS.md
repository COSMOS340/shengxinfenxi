# PC Fast-Merge Result Status

Task: `20260703_fig4_gse306459_snatac`

Run scope: `FIG4_RUN_SCOPE=fig4c`

Output directory: `fig4_snatac/outputs_fast_merge`

## Completed Outputs

- Fig4A UMAP label transfer: `figures/fig4a_snatac_umap_label_transfer.png` and `.pdf`
- Fig4B tissue proportion boxplots: `figures/fig4b_tissue_proportion_boxplots.png` and `.pdf`
- Fig4C ChIPSeeker broad annotation pie: `figures/fig4c_chipseeker_annotation_pie.png` and `.pdf`
- Core tables, checks, log, and scripts are included in this upload package.

## Run Metrics

- Cells: 122218
- Peaks after ChromatinAssay creation: 1307563
- Libraries: 15
- Samples: 45
- Transferred labels: 8
- Transfer features: 2675
- Differential peak rows: 32372
- Unique differential peaks: 32059
- ChIPSeeker broad annotation rows: 7

## Method Notes

- Raw GEO H5 files were extracted by HTTP Range from `GSE306459_RAW.tar`; the full 87.55 GB archive was not downloaded.
- `fig4c` scope used H5 peak matrices only. Fragments were not downloaded.
- Matrix merging used `SeuratObject:::RowMergeSparseMatrices` before creating one combined ChromatinAssay, avoiding slow object-by-object Seurat merge.
- Annotation used local `Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)` plus local seqlevel renaming, avoiding online UCSC seqlevel lookup.
- Gene activity used fragments-free peak-gene overlap against local gene annotations.
- Label transfer used `FindTransferAnchors(..., reduction="pcaproject", reference.reduction="pca")`.
- Differential peaks for Fig4C were selected by sparse fold-change filtering and top 5000 rows per transferred label. The `p_val` and `p_val_adj` fields are `NA` by design; this is not an LR p-value result.
- Fig4C plot uses broad ChIPSeeker categories. Full transcript-level ChIPSeeker annotations remain in `tables/fig4c_chipseeker_peak_annotation_full.tsv.gz`; the original detailed count table is preserved as `tables/fig4c_chipseeker_annotation_counts_detailed.tsv`.

## Checks

- `tables/fig4_postrun_checks.tsv`: all checks PASS.
- UMAP rows equal 122218.
- ChIPSeeker broad annotation fractions sum to 1.
- Figure visual QA:
  - Fig4A: title, axes, legend, and canvas passed.
  - Fig4B: title, axes, facet labels, legend, and canvas passed. It displays the labels with valid Wilcoxon comparisons.
  - Fig4C: original detailed annotation legend failed visual QA; it was regenerated using broad annotation categories and passed.

## Large Object Policy

The full processed Seurat RDS was not included. Whole-object `saveRDS()` produced a multi-GB file and was not suitable for this GitHub handoff package. The patched script defaults to `FIG4_SAVE_RDS=false`; set `FIG4_SAVE_RDS=true` to explicitly save the large object locally.
