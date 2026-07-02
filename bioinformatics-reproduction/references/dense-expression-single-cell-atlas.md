# Dense Expression Single-Cell Atlas First Pass

Use this route when public single-cell data are distributed as dense `csv.gz` expression matrices rather than 10X `matrix.mtx` directories.

## Required Input Contract

- One metadata file per dataset or sample group.
- One expression file per matching dataset or sample group.
- Expression file shape: cells as rows, genes as columns.
- First expression column: `index`.
- Metadata cell identifier column: `index`.
- Metadata columns for atlas auditing: `cancer`, `tissue`, `MajorCluster`, `UMAP1`, `UMAP2`.
- File naming convention handled by the template script:
  - `<prefix><dataset>_metadata.csv.gz`
  - `<prefix><dataset>_normalized_expression.csv.gz`

Do not use this route for 10X sparse matrices, raw Cell Ranger folders, AnnData objects, or spatial data without adapting the reader.

## Recommended First Pass

1. Stream-inspect expression files before matrix loading.
2. Verify row orientation, first-column name, metadata identifier alignment, row counts, gene counts, and common genes.
3. Use common genes across expression files for the first integrated audit.
4. If all-cell dense loading is too large for local memory, use deterministic stratified sampling.
5. Stratify on biologically relevant metadata such as `cancer + MajorCluster`, not only on file name.
6. Select highly variable genes from the sampled common-gene matrix.
7. Run PCA on scaled top-variable genes.
8. Build a kNN graph from PCs and run Louvain or Leiden clustering.
9. Run UMAP from PCs for visualization.
10. Compare expression-derived clusters to author or reference labels with adjusted Rand index and a contingency heatmap.

## Scientific Boundary

This route is an audit layer, not a final replacement for the paper's full integration. A stratified dense-matrix first pass can show whether expression structure is broadly compatible with published labels, but final lineage labels require marker review, batch correction or integration, repeated sampling or full-cell stability checks, and figure-by-figure comparison.

## Script

Use:

```bash
Rscript bioinformatics-reproduction/scripts/stratified_dense_expression_atlas.R \
  /path/to/raw_dir \
  /path/to/output_dir \
  GSE154763_
```

The third argument is optional and defaults to an empty prefix. Output folders are created under the supplied output directory.
