options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "Matrix", "Seurat", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    stop("Run this script with Rscript so --file is available.")
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

safe_write_png <- function(plot, path, width, height, dpi = 300) {
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

script_path <- get_script_path()
request_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
output_dir <- file.path(request_dir, "outputs")
object_dir <- file.path(request_dir, "local_objects_not_for_github")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

raw_root_default <- normalizePath(
  file.path(request_dir, "..", "20260706_core_geo_download_audit", "raw_core_geo"),
  mustWork = FALSE
)
raw_root <- Sys.getenv("RAW_CORE_GEO_DIR", unset = raw_root_default)
raw_root <- normalizePath(raw_root, mustWork = TRUE)
gse164522_dir <- file.path(raw_root, "GSE164522")
metadata_path <- file.path(gse164522_dir, "GSE164522_CRLM_metadata.csv.gz")
if (!file.exists(metadata_path)) {
  stop("Missing metadata file: ", metadata_path)
}

expression_files <- c(
  LN = "GSE164522_CRLM_LN_expression.csv.gz",
  MN = "GSE164522_CRLM_MN_expression.csv.gz",
  MT = "GSE164522_CRLM_MT_expression.csv.gz",
  PBMC = "GSE164522_CRLM_PBMC_expression.csv.gz",
  PN = "GSE164522_CRLM_PN_expression.csv.gz",
  PT = "GSE164522_CRLM_PT_expression.csv.gz"
)

expression_paths <- file.path(gse164522_dir, expression_files)
missing_expression <- expression_paths[!file.exists(expression_paths)]
if (length(missing_expression) > 0) {
  stop("Missing expression files: ", paste(missing_expression, collapse = "; "))
}

marker_path <- file.path(request_dir, "inputs", "mac_mono_marker_panels.tsv")
if (!file.exists(marker_path)) {
  stop("Missing marker panel input: ", marker_path)
}

meta <- fread(metadata_path, check.names = FALSE)
cell_id_col <- names(meta)[1]
required_meta_cols <- c(
  "n_genes", "percent_mito", "n_counts", "louvain", "sample", "patient",
  "tissue", "ID", "celltype_global", "celltype_major", "celltype_sub"
)
missing_meta_cols <- setdiff(required_meta_cols, names(meta))
if (length(missing_meta_cols) > 0) {
  stop("Metadata columns missing: ", paste(missing_meta_cols, collapse = ", "))
}

meta[, metadata_cell_id := get(cell_id_col)]
meta[, expression_cell_id := sub("-(\\d+)$", ".\\1", metadata_cell_id)]

fwrite(
  meta[, .N, by = .(celltype_major, celltype_sub, tissue)][order(celltype_major, celltype_sub, tissue)],
  file.path(output_dir, "gse164522_metadata_celltype_counts.tsv"),
  sep = "\t"
)

target_major <- c("Macrophage", "Monocyte")
selected_meta <- meta[celltype_major %in% target_major]
if (nrow(selected_meta) == 0) {
  stop("No cells found for target celltype_major values: ", paste(target_major, collapse = ", "))
}

matrices <- list()
metadata_parts <- list()
selection_audit <- list()
gene_audit <- list()

for (label in names(expression_paths)) {
  path <- expression_paths[[label]]
  headers <- names(fread(path, nrows = 0, check.names = FALSE))
  if (length(headers) < 2) {
    stop("Expression file has fewer than two columns: ", path)
  }
  gene_col <- headers[[1]]
  expression_columns <- headers[-1]
  overlap <- intersect(selected_meta$expression_cell_id, expression_columns)
  selection_audit[[label]] <- data.table(
    expression_group = label,
    file_name = basename(path),
    expression_columns = length(expression_columns),
    selected_metadata_cells = nrow(selected_meta),
    selected_overlap_cells = length(overlap)
  )
  if (length(overlap) == 0) {
    next
  }

  message("Reading ", basename(path), " selected columns: ", length(overlap))
  dt <- fread(path, select = c(gene_col, overlap), check.names = FALSE, showProgress = TRUE)
  genes <- dt[[gene_col]]
  duplicate_gene_count <- sum(duplicated(genes))
  if (duplicate_gene_count > 0) {
    genes <- make.unique(genes)
  }
  gene_audit[[label]] <- data.table(
    expression_group = label,
    file_name = basename(path),
    genes = length(genes),
    duplicate_gene_names = duplicate_gene_count
  )

  dense <- as.matrix(dt[, ..overlap])
  storage.mode(dense) <- "double"
  rownames(dense) <- genes
  sparse <- Matrix(dense, sparse = TRUE)
  rm(dense, dt)
  gc()

  matrices[[label]] <- sparse
  metadata_parts[[label]] <- selected_meta[match(overlap, expression_cell_id)]
  metadata_parts[[label]][, expression_group := label]
}

selection_audit_dt <- rbindlist(selection_audit, use.names = TRUE, fill = TRUE)
fwrite(selection_audit_dt, file.path(output_dir, "gse164522_full_mac_mono_selection_audit.tsv"), sep = "\t")
fwrite(rbindlist(gene_audit, use.names = TRUE, fill = TRUE), file.path(output_dir, "gse164522_full_mac_mono_gene_audit.tsv"), sep = "\t")

if (length(matrices) == 0) {
  stop("No selected cells overlapped expression columns in any expression file.")
}

common_genes <- Reduce(intersect, lapply(matrices, rownames))
if (length(common_genes) == 0) {
  stop("No common genes across selected expression matrices.")
}
matrices <- lapply(matrices, function(x) x[common_genes, , drop = FALSE])
combined <- do.call(cbind, matrices)
combined_meta <- rbindlist(metadata_parts, use.names = TRUE, fill = TRUE)
combined_meta <- as.data.frame(combined_meta)
rownames(combined_meta) <- combined_meta$expression_cell_id
combined_meta <- combined_meta[colnames(combined), , drop = FALSE]

fwrite(
  as.data.table(combined_meta, keep.rownames = "expression_cell_id"),
  file.path(output_dir, "gse164522_full_mac_mono_selected_cells.tsv"),
  sep = "\t"
)

count_by_tissue <- as.data.table(combined_meta)[, .N, by = .(tissue, celltype_major, celltype_sub, patient, expression_group)]
fwrite(
  count_by_tissue[order(tissue, celltype_major, celltype_sub, patient)],
  file.path(output_dir, "gse164522_full_mac_mono_counts_by_tissue_patient.tsv"),
  sep = "\t"
)

obj <- CreateSeuratObject(
  counts = combined,
  project = "GSE164522_CRLM_full_mac_mono",
  meta.data = combined_meta,
  min.cells = 0,
  min.features = 0
)
rm(combined, matrices)
gc()

obj[["percent_mt_calculated"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- NormalizeData(obj, verbose = TRUE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = TRUE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = TRUE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = TRUE)
dims_use <- seq_len(min(30, ncol(Embeddings(obj, "pca"))))
obj <- FindNeighbors(obj, dims = dims_use, verbose = TRUE)
obj <- FindClusters(obj, resolution = c(0.4, 0.6, 0.8), verbose = TRUE)
cluster_col <- "RNA_snn_res.0.6"
if (!cluster_col %in% colnames(obj[[]])) {
  stop("Expected clustering column not found: ", cluster_col)
}
Idents(obj) <- cluster_col
obj <- RunUMAP(obj, dims = dims_use, seed.use = 20260706, verbose = TRUE)

umap <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "expression_cell_id")
if (all(c("UMAP_1", "UMAP_2") %in% names(umap))) {
  setnames(umap, old = c("UMAP_1", "UMAP_2"), new = c("umap_1", "umap_2"))
}
umap_meta <- as.data.table(obj[[]], keep.rownames = "expression_cell_id")
fwrite(
  merge(umap, umap_meta, by = "expression_cell_id", all.x = TRUE, sort = FALSE),
  file.path(output_dir, "gse164522_full_mac_mono_umap_coordinates.tsv.gz"),
  sep = "\t"
)

if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  obj <- JoinLayers(obj)
}
Idents(obj) <- cluster_col
markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
markers_dt <- as.data.table(markers)
fwrite(markers_dt, file.path(output_dir, "gse164522_full_mac_mono_cluster_markers.tsv.gz"), sep = "\t")
top50 <- markers_dt[order(cluster, p_val_adj, -avg_log2FC), head(.SD, 50), by = cluster]
fwrite(top50, file.path(output_dir, "gse164522_full_mac_mono_cluster_top50_markers.tsv"), sep = "\t")

marker_panels <- fread(marker_path)
marker_coverage <- copy(marker_panels)
marker_coverage[, present := gene %in% rownames(obj)]
fwrite(marker_coverage, file.path(output_dir, "gse164522_full_mac_mono_marker_panel_coverage.tsv"), sep = "\t")
dotplot_features <- unique(marker_coverage[present == TRUE, gene])
if (length(dotplot_features) > 0) {
  dot <- DotPlot(obj, features = dotplot_features, group.by = cluster_col, cols = c("#edf2f7", "#b91c1c")) +
    RotatedAxis() +
    theme_bw(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      axis.title = element_blank()
    )
  safe_write_png(dot, file.path(output_dir, "fig1_gse164522_full_mac_mono_marker_dotplot.png"), 14, 7)
}

umap_cluster <- DimPlot(obj, reduction = "umap", group.by = cluster_col, label = TRUE, repel = TRUE, pt.size = 0.18) +
  ggtitle("GSE164522 Macrophage/Monocyte clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 11))
umap_tissue <- DimPlot(obj, reduction = "umap", group.by = "tissue", pt.size = 0.18) +
  ggtitle("Tissue") +
  theme(plot.title = element_text(hjust = 0.5, size = 11))
umap_major <- DimPlot(obj, reduction = "umap", group.by = "celltype_major", pt.size = 0.18) +
  ggtitle("Author major label") +
  theme(plot.title = element_text(hjust = 0.5, size = 11))
layout_plot <- (umap_cluster | umap_tissue | umap_major) + plot_layout(widths = c(1, 1, 1))
safe_write_png(layout_plot, file.path(output_dir, "fig1_gse164522_full_mac_mono_umap_layout.png"), 15, 5)
safe_write_png(umap_cluster, file.path(output_dir, "fig1_gse164522_full_mac_mono_umap_clusters.png"), 6, 5)
safe_write_png(umap_tissue, file.path(output_dir, "fig1_gse164522_full_mac_mono_umap_tissue.png"), 6, 5)
safe_write_png(umap_major, file.path(output_dir, "fig1_gse164522_full_mac_mono_umap_author_major.png"), 6, 5)

cluster_tissue <- as.data.table(obj[[]])[, .N, by = c(cluster_col, "tissue")]
setnames(cluster_tissue, cluster_col, "cluster")
cluster_tissue[, fraction := N / sum(N), by = cluster]
fwrite(cluster_tissue, file.path(output_dir, "gse164522_full_mac_mono_cluster_tissue_composition.tsv"), sep = "\t")
comp_plot <- ggplot(cluster_tissue, aes(x = cluster, y = fraction, fill = tissue)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.15) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 9) +
  theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()) +
  labs(y = "Cell fraction", fill = "Tissue")
safe_write_png(comp_plot, file.path(output_dir, "fig1_gse164522_full_mac_mono_cluster_tissue_composition.png"), 8, 4.5)

top_heatmap_features <- unique(top50[order(cluster, p_val_adj, -avg_log2FC), head(.SD, 5), by = cluster]$gene)
top_heatmap_features <- top_heatmap_features[top_heatmap_features %in% rownames(obj)]
if (length(top_heatmap_features) > 1) {
  heat <- DoHeatmap(obj, features = top_heatmap_features, group.by = cluster_col, size = 2.3) +
    NoLegend()
  safe_write_png(heat, file.path(output_dir, "fig1_gse164522_full_mac_mono_top_marker_heatmap.png"), 12, 9)
}

object_path <- file.path(object_dir, "gse164522_full_mac_mono_seurat.rds")
saveRDS(obj, object_path, compress = "xz")
object_size <- file.info(object_path)$size

run_summary <- data.table(
  item = c(
    "raw_root",
    "metadata_path",
    "selected_major_labels",
    "selected_cells",
    "genes",
    "clusters_res_0_6",
    "object_path_local_only",
    "object_bytes_local_only"
  ),
  value = c(
    raw_root,
    metadata_path,
    paste(target_major, collapse = ","),
    as.character(ncol(obj)),
    as.character(nrow(obj)),
    as.character(length(levels(Idents(obj)))),
    object_path,
    as.character(object_size)
  )
)
fwrite(run_summary, file.path(output_dir, "gse164522_full_mac_mono_run_summary.tsv"), sep = "\t")

writeLines(
  c(
    "# GSE164522 full Macrophage/Monocyte analysis",
    "",
    paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    paste0("Selected cells: ", ncol(obj)),
    paste0("Genes: ", nrow(obj)),
    paste0("Resolution 0.6 clusters: ", length(levels(Idents(obj)))),
    "",
    "Raw matrices and the Seurat RDS object should not be committed to GitHub.",
    "Commit only the generated outputs directory and this request folder."
  ),
  file.path(output_dir, "STATUS.md")
)

sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
