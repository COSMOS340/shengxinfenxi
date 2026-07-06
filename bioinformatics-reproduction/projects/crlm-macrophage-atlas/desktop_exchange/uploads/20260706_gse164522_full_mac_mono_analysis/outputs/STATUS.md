# GSE164522 full Macrophage/Monocyte analysis

Completed: 2026-07-06T20:36:22+0800
Selected cells: 15853
Genes: 24662
Resolution 0.6 clusters: 23

Raw matrices and the Seurat RDS object should not be committed to GitHub.
Commit only the packaged upload directory.

PC local fixes applied before the successful run:
- Preserved expression group names after `file.path()` so LN/MN/MT/PBMC/PN/PT were iterated correctly.
- Avoided duplicated `expression_cell_id` when exporting UMAP coordinates.
- Regenerated the marker dotplot and heatmap on wider canvases after visual QC.

PC visual QC:
- PNG outputs were opened locally.
- Plots were not blank.
- Axes, legends, and main labels were not clipped.
- The heatmap rightmost cluster labels are tight because several clusters are small, but the matrix and labels remain visible.
