# PC request: rebuild GSE205506 with the source-paper QC and UMAP palette

Status: `REQUESTED_GSE205506_PAPER_QC_REBUILD`

Date: 2026-07-17

## Purpose

Rebuild the formal GSE205506 R/Seurat object from the 40 author matrices because the previous 238,934-cell object used a per-sample mitochondrial filter that is not the method reported in the source article. The rebuilt object will be the CRC input for Figure 1.

Do not perform response differential testing, cell communication analysis, pseudobulk analysis, or final manuscript composition figures in this request.

## Inputs that have already passed review

- Root: `F:/pan-gi-ici-ccc-20260709/GSE205506`
- GEO archive: `F:/pan-gi-ici-ccc-20260709/GSE205506/raw/GSE205506_RAW.tar`
- Table S1: `F:/pan-gi-ici-ccc-20260709/GSE205506/supplement/mmc2.xlsx`
- Author marker workbook: `F:/pan-gi-ici-ccc-20260709/GSE205506/supplement/mmc3.xlsx`
- Existing exact response mapping: use the already audited 40-sample mapping from the completed request.
- Previous upload for reusable metadata and marker references: `uploads/20260716_stop_prjna932556_and_build_gse205506`

Keep the previous RDS files and label them as rejected QC version 1 in the new inventory. Do not overwrite or delete them.

## Mandatory language and provenance rules

1. All expression analysis, QC, integration, clustering, marker calculation, and plotting must run in R.
2. Do not use Scanpy, Python UMAP, Python clustering, or inherited coordinates.
3. Record exact R, Seurat, Matrix, SingleCellExperiment, scDblFinder, ggplot2, and integration package versions.
4. Preserve exact `geo_accession`, subject, tissue, treatment, timepoint, and Table S1 response fields for every cell.
5. Do not infer or rename a source field without showing the exact input field and transformation in an audit table.

## 1. Re-import and source-paper basic QC

Start from all 40 author matrix, barcode, and feature triplets. Do not start from `gse205506_seurat_qc_singlets.rds`, `gse205506_seurat_integrated_clusters.rds`, or either annotated RDS.

Apply the numerical rules reported in Li et al., Cancer Cell 2023, DOI `10.1016/j.ccell.2023.04.011`:

- remove genes detected in fewer than three cells;
- retain cells with `nFeature_RNA >= 500` and `nFeature_RNA <= 5000`;
- retain cells with `nCount_RNA >= 400` and `nCount_RNA <= 25000`.

The article does not state whether the fewer-than-three-cells gene rule was applied per sample or after merging. Use the merged-object interpretation for the formal route and report the per-sample interpretation as a gene-count sensitivity check. Do not describe either interpretation as author code. This ambiguity does not directly change cell inclusion because it is a gene filter.

Export one row per cell before filtering with the exact QC metrics and pass/fail reason. Also export per-sample and project-total stage counts.

## 2. Doublet handling must match the article

The article does not report scDblFinder, DoubletFinder, Scrublet, or another explicit computational caller. It states that the upper gene and UMI limits remove most barcodes associated with doublet cells.

Therefore:

- the formal primary object must use the article's gene and UMI bounds without scDblFinder-based exclusion;
- run scDblFinder per sample only as a secondary sensitivity analysis;
- export the overlap between article-QC pass cells and scDblFinder calls;
- compare broad-cell and cluster proportions before and after the sensitivity exclusion;
- never call the scDblFinder-filtered object the source-paper primary object.

## 3. First-round integration and six broad compartments

Reproduce the reported route as closely as the available Seurat version allows:

1. `NormalizeData` with default normalization and scale factor 10,000.
2. Select 2,000 variable genes using `FindVariableFeatures(selection.method = "vst")`.
3. Use Seurat RPCA integration through `FindIntegrationAnchors(reduction = "rpca")` and `IntegrateData`; record all arguments and dimensions.
4. Regress `nCount_RNA` in scaling, matching the paper's regression of UMI counts.
5. Run PCA and use PCs 1:20 for neighbors, first-round clustering, and UMAP.
6. Use `FindClusters(resolution = 1.2)` for the first-round broad-cell clustering.
7. Use the integrated assay only for clustering and broad-cell classification. Use RNA counts or normalized RNA expression for markers and downstream expression summaries.

Assign exactly the six broad compartments used in the article, based on cluster Top10 and Top50 markers plus the reported marker sets:

- T/I/NK: `CD3D`, `CD3E`, `TRAC`, `TRBC1`
- B: `CD79A`, `CD79B`, `MS4A1`, `TNFRSF17`, `MZB1`
- Myeloid: `CD14`, `CD68`
- Epithelial: `EPCAM`, `CD24`
- Fibroblast: `COL1A2`, `COL3A1`, `MYH11`, `ACTA2`
- Endothelial: `VWF`, `PECAM1`

For every first-round cluster, export positive evidence, negative evidence, Top10, Top50, sample/subject distribution, and broad-marker co-expression. Do not delete a cluster merely because it is small or because it has a high mitochondrial median.

The paper reports removal of extremely low-abundance and broad-marker-conflicting clusters but gives no numerical rule. Do not invent a threshold or silently remove these clusters. Flag them in a review table and retain them in a review object. Any exclusion from the formal object must have an explicit reproducible rule and a before/after cell list.

## 4. Compartment-aware mitochondrial QC

After the first-round broad assignment, apply the source-paper strategy at the cell level.

### Epithelial

- remove cells with `percent.mt > 75`;
- retain cells with `percent.mt <= 75`.

### Lymphoid, Myeloid, Fibroblast, and Endothelial

For the formal route, build four model groups: lymphoid is the combined T/I/NK plus B cells, followed by separate Myeloid, Fibroblast, and Endothelial groups. The article reports one combined lymphoid removal fraction for T and B cells; do not fit separate T and B primary models. A separate T-versus-B fit may be exported only as a clearly labeled sensitivity analysis.

For each of the four formal model groups, implement the paper's median-centered, MAD-variance normal model as the following reproducible operationalization:

```r
center <- median(percent.mt, na.rm = TRUE)
sigma <- mad(percent.mt, center = center, constant = 1.4826, na.rm = TRUE)
p_upper <- pnorm(percent.mt, mean = center, sd = sigma, lower.tail = FALSE)
p_bonferroni <- p.adjust(p_upper, method = "bonferroni")
mt_qc_pass <- p_bonferroni >= 0.05
```

This is a project implementation of the article's written method, not a claim that the unavailable author code used these exact R statements.

Required safeguards:

- stop with an explicit error if `sigma` is zero or non-finite; do not substitute a guessed threshold;
- export `center`, `sigma`, cell count, raw p value, adjusted p value, and pass/fail for every compartment and cell;
- report removal fractions and compare them with the article's reported 9.20% lymphoid, 12.84% myeloid, 8.11% fibroblast, and 8.50% endothelial values;
- report epithelial removal fraction and compare it with the article's 29.75%;
- do not tune parameters, downsample, or remove extra cells to force the article's final total of 155,397.

## 5. Rebuild the formal atlas after QC

From cells passing the source-paper basic and compartment-aware mitochondrial QC:

1. repeat the RPCA integration and 20-PC route with the same recorded parameters;
2. rebuild the all-cell UMAP and broad labels;
3. export project, sample, subject, response, timepoint, tissue, and broad-cell counts;
4. quantify sample structure before and after integration using the same KNN and PC eta-squared audits used previously;
5. export the exact cell IDs retained and excluded at each stage.

Compare the final total with the article's 155,397 cells. A difference must be explained by stage, sample, and compartment. Do not force equality.

## 6. Compartment-specific reclustering and annotation review

After the formal broad atlas is complete, recluster these five compartments separately, matching the structure of the author marker workbook:

- T/I/NK
- B
- Myeloid
- Endothelial
- Fibroblast

Keep epithelial cells as a separate sixth analysis object. The author workbook has no epithelial sheet, so do not assign epithelial subtypes from that workbook.

For each compartment:

- rerun variable-feature selection, PCA, neighbors, UMAP, and multiple clustering resolutions in R;
- calculate each cluster's own Top10 and Top50 markers from RNA expression;
- compare with the matching sheet in `mmc3.xlsx`;
- use canonical positive and negative markers as an independent check;
- export one evidence row per cluster with source markers, author overlap, conflicts, sample/subject support, and final reviewed label;
- do not transfer the 45 author labels directly from whole-object clusters;
- do not force a subtype label when the marker evidence conflicts.

At subtype level, the article reports removal of mitochondrial, ribosomal, or hybrid clusters but does not provide complete numerical rules for every compartment. Flag these states and export them for Mac review. Do not silently remove whole subtype clusters in this request.

## 7. Required UMAP palette and visual design

The previous dark palette is rejected. Use the exact bright color family already used for the first R/Seurat dataset.

```r
broad_palette <- c(
  "Epithelial" = "#B79F00",
  "T/I/NK" = "#0072B2",
  "B" = "#56B4E9",
  "Myeloid" = "#D55E00",
  "Endothelial" = "#17BECF",
  "Fibroblast" = "#8C564B",
  "Hybrid_or_unresolved_review" = "#7A5195"
)

response_palette <- c("pCR" = "#009E73", "non-pCR" = "#D55E00")
timepoint_palette <- c("pre-treatment" = "#0072B2", "post-treatment" = "#D55E00")
```

Use `#CC79A7` for plasma subtype, `#009E73` for macrophage subtype where needed, and `#BDBDBD` only for QC-failed review views. Preserve the same color for the same broad cell type across all datasets and panels.

Rendering requirements:

- white background;
- no dark brown, dark red, or dark purple all-cell palette;
- randomize plotting order with a recorded seed;
- rasterize points but export both PNG and vector-container PDF;
- make labels readable without opaque label boxes covering cell structure;
- inspect the final PNGs for clipping, text overlap, point washout, legend overflow, and insufficient color separation.

## 8. Required outputs

Upload a new directory:

`uploads/20260717_gse205506_paper_qc_rebuild`

At minimum include:

- `STATUS.md`
- `MAC_REVIEW_REQUEST.md`
- all R scripts used in this rebuild
- all stdout/stderr logs and `sessionInfo()` files
- `gse205506_paper_qc_parameters.tsv`
- `gse205506_paper_qc_cell_audit.tsv.gz`
- `gse205506_paper_qc_stage_counts.tsv`
- `gse205506_paper_qc_sample_counts.tsv`
- `gse205506_paper_qc_compartment_mt_model.tsv`
- `gse205506_paper_qc_removal_fraction_comparison.tsv`
- `gse205506_paper_qc_doublet_sensitivity.tsv`
- `gse205506_paper_qc_broad_cluster_evidence.tsv`
- `gse205506_paper_qc_hybrid_review.tsv`
- `gse205506_paper_qc_final_cell_ids.tsv.gz`
- `gse205506_paper_qc_final_object_summary.tsv`
- `gse205506_paper_qc_batch_knn_audit.tsv`
- `gse205506_paper_qc_batch_pc_eta_squared.tsv`
- all-cell Top10 and Top50 marker tables
- one Top10 and one Top50 table for every reclustered compartment
- one author-overlap and annotation-evidence table for every compartment
- broad UMAP, response UMAP, timepoint UMAP, sample UMAP, QC before/after plots, and five compartment UMAPs as PNG and PDF
- marker dotplots for the broad atlas and five compartments as PNG and PDF
- `output_manifest.tsv` with relative path, size, SHA-256, delivery status, and notes

Keep large RDS objects on F: and list their exact paths, sizes, and SHA-256 values:

- article-basic-QC object before broad mitochondrial filtering;
- first-round integrated broad-classification object;
- formal post-mitochondrial-QC all-cell object;
- each of the five compartment objects;
- separate epithelial object;
- scDblFinder sensitivity object if saved.

## Completion gate

The request is complete only when:

1. all 40 samples and all 19 subjects remain traceable;
2. the primary cell inclusion follows the article's numerical rules;
3. scDblFinder is sensitivity-only and is not mixed into the source-paper primary object;
4. the compartment-aware mitochondrial model has a cell-level audit;
5. the six broad compartments are supported by each cluster's own Top10/Top50 markers;
6. the five author compartments have been reclustered separately;
7. no whole cluster has been silently deleted;
8. all required figures use the approved bright palette and pass visual inspection;
9. the output manifest and local-only RDS inventory are complete.
