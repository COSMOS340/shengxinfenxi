# Single-Cell, Spatial, Trajectory, Communication, and Perturbation Workflows

Use this reference when a reproduction task uses single-cell RNA-seq, spatial transcriptomics, cell trajectory analysis, cell-cell communication, target-gene localization, or virtual gene perturbation.

## Route Selection

| Study claim | Minimum route |
|---|---|
| A bulk-derived gene is expressed in a disease-relevant cell type | Seurat QC, clustering, annotation, marker check, target-gene expression plots |
| A target gene varies by disease group within a cell type | Verified sample metadata, cell-type annotation, per-cell or pseudo-bulk differential testing |
| A gene is linked to cell-state ordering | Monocle3 or another trajectory method with a justified root state |
| A tissue-localized signal is claimed | Spatial transcriptomics with image-coordinate alignment and spatial feature plots |
| Intercellular communication is claimed | CellChat, NicheNet, CellPhoneDB, or another complete ligand-receptor workflow |
| A gene perturbation effect is claimed from single-cell data | scTenifoldKnk or another declared simulation workflow, reported as virtual perturbation |

## Data Inputs

### 10X Directory

Use this route when each sample has `barcodes.tsv.gz`, `features.tsv.gz` or `genes.tsv.gz`, and `matrix.mtx.gz`.

Minimum steps:

1. Read each sample with `Read10X`.
2. Prefix cell barcodes with sample IDs before merging.
3. Store sample ID, disease group, tissue, batch, platform, and accession in metadata.
4. Create a Seurat object with explicit `min.cells` and `min.features`.
5. Save a raw Seurat object or sparse matrix before filtering.

### Dense CSV or TSV

Use this route when authors provide an expression table.

Minimum checks:

1. Verify whether rows are cells or genes.
2. Verify gene identifiers and species.
3. Convert to sparse matrix when the table is large.
4. Resolve duplicated gene identifiers with a declared rule.
5. Preserve original cell IDs.

### Spatial Transcriptomics

For Visium-style data, verify:

1. Expression matrix and barcodes.
2. High-resolution and low-resolution tissue images.
3. Scale-factor JSON.
4. Tissue-position file.
5. Barcode alignment between expression and image coordinates.

## Seurat QC

Minimum QC outputs:

1. Cell count and gene count per sample before filtering.
2. `nFeature_RNA`, `nCount_RNA`, and mitochondrial fraction distributions.
3. Scatter plots for UMI count vs mitochondrial fraction and UMI count vs detected genes.
4. Cells retained per sample after filtering.
5. Doublet summary when doublet removal is used.
6. Ambient RNA or contamination summary when correction is used.

Common steps:

```r
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > min_features & percent.mt < max_mt)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_hvg)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = dims_use)
obj <- FindClusters(obj, resolution = resolution)
obj <- RunUMAP(obj, dims = dims_use)
```

Guardrails:

1. Choose filtering thresholds from the current dataset distribution.
2. Do not copy thresholds from another study without justification.
3. Re-run normalization and PCA after removing doublets or contaminated cells.
4. Compare batches and samples on UMAP before biological interpretation.

## Clustering Resolution

Do not report a single arbitrary clustering resolution without diagnostics.

Recommended checks:

1. Run multiple resolutions.
2. Record cluster count and smallest cluster size.
3. Plot a clustree or equivalent cluster transition graph.
4. Check known marker genes per cluster.
5. Merge or relabel clusters only with marker support.

## Annotation

Reference annotation methods such as SingleR are useful first-pass tools.

Minimum route:

1. Run cluster-level annotation when clusters are stable.
2. Save cluster-to-cell-type mapping.
3. Save per-cell mapping.
4. Validate labels with canonical marker genes.
5. Mark weak or mixed clusters explicitly.

Example:

```r
ref <- celldex::HumanPrimaryCellAtlasData()
expr <- GetAssayData(obj, layer = "data")
ann <- SingleR(test = expr, ref = ref, labels = ref$label.main, clusters = obj$seurat_clusters)
```

Do not rely on reference labels alone for manuscript claims.

## Marker Genes and Cell-Type Differential Testing

Cluster markers:

```r
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = logfc_cutoff)
markers_sig <- markers[markers$p_val_adj < padj_cutoff, ]
```

Cell-type differential expression:

1. Use verified sample metadata to define groups.
2. Avoid deriving disease labels from barcode strings.
3. Use pseudo-bulk by sample and cell type for manuscript-grade differential claims when sample replication exists.
4. Reserve per-cell Wilcoxon tests for exploratory screens or visualization unless the study design justifies cell-level testing.
5. Correct for multiple testing across genes and cell types.

## Target-Gene Localization

For genes selected by bulk analysis, MR, ML, WGCNA, network toxicology, or docking:

1. Verify gene symbols exist in the single-cell object.
2. Plot expression on UMAP.
3. Plot expression by cluster and annotated cell type.
4. Export average expression by cluster and cell type.
5. Export per-cell expression with sample, group, cluster, and cell type.
6. Check whether expression is driven by one sample or batch.

High/low splits by median expression are visualization groupings. State this boundary in manuscripts.

## Monocle3 Trajectory

Use trajectory analysis only when there is a biological ordering question.

Minimum route:

```r
cds <- new_cell_data_set(expr_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
cds <- preprocess_cds(cds, num_dim = n_dim)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)
```

Guardrails:

1. Choose root cells from biology, marker genes, time point, or known lineage start.
2. Record the root-cell rule.
3. Do not interpret pseudotime as chronological time.
4. Plot trajectory by pseudotime, cell type, sample, and disease group.
5. Check whether inferred ordering is driven by batch.

## Spatial Transcriptomics

Minimum Seurat route:

```r
img <- Read10X_Image(image.dir = "spatial", filter.matrix = TRUE)
sp <- CreateSeuratObject(counts = expr_matrix, assay = "Spatial", project = project_id)
img <- img[Cells(sp)]
sp[["slice"]] <- img
sp <- subset(sp, subset = nCount_Spatial > 0 & nFeature_Spatial > 0)
sp <- SCTransform(sp, assay = "Spatial")
sp <- RunPCA(sp, assay = "SCT")
sp <- FindNeighbors(sp, reduction = "pca", dims = 1:30)
sp <- FindClusters(sp, resolution = resolution)
sp <- RunUMAP(sp, reduction = "pca", dims = 1:30)
```

Minimum outputs:

1. Spatial QC feature plots.
2. UMAP clusters.
3. Spatial cluster maps.
4. Spatial marker genes.
5. Target-gene spatial feature plots.
6. Exported metadata and coordinates.

Spatial annotation guardrails:

1. Spots are not single cells in standard Visium.
2. Cell-type labels assigned to spots are deconvolution or annotation results.
3. Check tissue image alignment before interpreting spatial enrichment.
4. Report whether a cluster label is dominant or mixed.

## Spatial Pseudotime

Spatial pseudotime can be used for exploratory tissue gradients.

Minimum checks:

1. Add pseudotime values back to the spatial object.
2. Plot pseudotime on tissue coordinates.
3. Compare pseudotime with clusters and anatomy.
4. Avoid directional claims without histology, time course, or lineage evidence.

## Cell-Cell Communication

Loading a package is not communication analysis. A complete workflow must include ligand-receptor inference and network outputs.

Minimum CellChat route:

```r
cellchat <- createCellChat(object = expr, meta = meta, group.by = "cell_type")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
```

Minimum outputs:

1. Communication probability table.
2. Pathway-level communication table.
3. Network circle or heatmap.
4. Sender/receiver role analysis.
5. Ligand-receptor pairs supporting each claimed interaction.

Guardrails:

1. Cell types with very few cells should be filtered.
2. Compare groups only after aligning cell-type composition and sample metadata.
3. Do not claim protein-level signaling from transcript-only ligand-receptor inference.

## scTenifoldKnk Virtual Perturbation

Use scTenifoldKnk to simulate gene knockout effects on inferred regulatory networks.

Minimum route:

1. Confirm the target gene is present.
2. Select highly variable genes.
3. Build a count matrix containing the target gene and selected genes.
4. Run scTenifoldKnk with recorded network and cell-sampling parameters.
5. Export all differential regulation results and significant results.
6. Plot adjusted p values and top regulated genes.

Interpretation boundary:

1. This is virtual perturbation.
2. It supports mechanism prioritization.
3. It does not replace CRISPR, RNAi, or drug perturbation experiments.

## Publication Checklist

Before reporting single-cell or spatial results:

1. Dataset accessions and sample metadata are verified.
2. Raw and processed object paths are recorded.
3. QC thresholds and retained-cell counts are reported.
4. Doublet and contamination handling is documented.
5. Batch structure is inspected.
6. Annotation has marker support.
7. Group comparisons use sample-level replication when claims are inferential.
8. Trajectory root state is justified.
9. Spatial image and coordinates align.
10. Communication analysis includes real ligand-receptor modeling.
11. Virtual perturbation is labeled as simulation.
12. All figures have been visually inspected for overlap, clipping, and unreadable labels.
