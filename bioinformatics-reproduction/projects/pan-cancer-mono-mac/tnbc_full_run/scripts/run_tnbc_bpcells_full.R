#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BPCells)
  library(Seurat)
  library(harmony)
  library(data.table)
  library(ggplot2)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- normalizePath(sub(file_arg, "", script_args[grepl(file_arg, script_args)][1]), mustWork = FALSE)

find_project_root <- function(start_dir) {
  current <- normalizePath(start_dir, mustWork = FALSE)
  for (i in 1:8) {
    package_matrix <- file.path(current, "inputs", "GSE169246_TNBC_RNA.counts.mtx.gz")
    project_matrix <- file.path(current, "03_data_raw", "single_cell", "GSE169246", "GSE169246_TNBC_RNA.counts.mtx.gz")
    if (file.exists(package_matrix) || file.exists(project_matrix)) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  stop("Could not locate TNBC project/package root from: ", start_dir)
}

start_dir <- if (nzchar(script_path)) dirname(script_path) else getwd()
project_root <- Sys.getenv("TNBC_PROJECT_ROOT", unset = find_project_root(start_dir))

package_input_dir <- file.path(project_root, "inputs")
package_layout <- file.exists(file.path(package_input_dir, "GSE169246_TNBC_RNA.counts.mtx.gz"))
default_input_dir <- if (package_layout) package_input_dir else file.path(project_root, "03_data_raw", "single_cell", "GSE169246")
default_metadata_path <- if (package_layout) {
  file.path(package_input_dir, "tnbc_gse169246_sample_metadata.tsv")
} else {
  file.path(project_root, "00_metadata", "tnbc_gse169246_sample_metadata.tsv")
}
default_output_root <- if (package_layout) file.path(project_root, "outputs") else project_root

input_dir <- Sys.getenv("TNBC_INPUT_DIR", unset = default_input_dir)
metadata_path <- Sys.getenv("TNBC_SAMPLE_METADATA", unset = default_metadata_path)
out_dir <- Sys.getenv("TNBC_OUTPUT_DIR", unset = file.path(default_output_root, "04_data_processed", "single_cell", "full_scope_preprocess", "tnbc_bpcells_full"))
fig_dir <- Sys.getenv("TNBC_FIGURE_DIR", unset = file.path(default_output_root, "06_figures", "main", "full_scope_preprocess", "tnbc_bpcells_full"))
tab_dir <- Sys.getenv("TNBC_TABLE_DIR", unset = file.path(default_output_root, "07_tables", "main", "full_scope_preprocess", "tnbc_bpcells_full"))

matrix_path <- file.path(input_dir, "GSE169246_TNBC_RNA.counts.mtx.gz")
barcode_path <- file.path(input_dir, "GSE169246_TNBC_RNA.barcode.tsv.gz")
feature_path <- file.path(input_dir, "GSE169246_TNBC_RNA.feature.tsv.gz")
bpcells_dir <- file.path(out_dir, "GSE169246_TNBC_RNA_counts_bpcells")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tmp"), recursive = TRUE, showWarnings = FALSE)

set.seed(20260702)
threads <- as.integer(Sys.getenv("TNBC_THREADS", unset = "4"))
npcs <- as.integer(Sys.getenv("TNBC_NPCS", unset = "30"))
resolution <- as.numeric(Sys.getenv("TNBC_RESOLUTION", unset = "0.1"))
marker_max_cells <- as.integer(Sys.getenv("TNBC_MARKER_MAX_CELLS_PER_CLUSTER", unset = "0"))

required_files <- c(matrix_path, barcode_path, feature_path, metadata_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required files: ", paste(missing_files, collapse = "; "))
}

message("Reading TNBC barcodes, features, and sample metadata")
barcodes <- fread(barcode_path, header = FALSE, col.names = "cell_id")[["cell_id"]]
features <- fread(feature_path, header = FALSE, col.names = "gene")[["gene"]]
sample_metadata <- fread(metadata_path)
sample_metadata <- sample_metadata[is_atac_title == "FALSE" & Sample_library_strategy == "RNA-Seq"]
sample_metadata[, sample_id_from_title := as.character(sample_id_from_title)]

cell_metadata <- data.table(cell_id = barcodes)
cell_metadata[, sample_id_from_title := sub("^[^.]+\\.", "", cell_id)]
cell_metadata <- merge(
  cell_metadata,
  sample_metadata,
  by = "sample_id_from_title",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(cell_metadata[["Sample_geo_accession"]])) {
  missing_samples <- unique(cell_metadata[is.na(Sample_geo_accession), sample_id_from_title])
  stop("Barcode sample suffixes missing from GEO sample metadata: ", paste(missing_samples, collapse = "; "))
}
cell_metadata <- as.data.frame(cell_metadata)
rownames(cell_metadata) <- cell_metadata[["cell_id"]]
cell_metadata <- cell_metadata[barcodes, , drop = FALSE]

sample_summary <- as.data.frame(table(
  sample_id_from_title = cell_metadata[["sample_id_from_title"]],
  patient_from_title = cell_metadata[["patient_from_title"]],
  tissue = cell_metadata[["tissue"]],
  group = cell_metadata[["group"]]
), stringsAsFactors = FALSE)
sample_summary <- sample_summary[sample_summary$Freq > 0, ]
write.table(sample_summary, file.path(tab_dir, "tnbc_full_sample_cell_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Importing or opening BPCells matrix")
if (!dir.exists(bpcells_dir) || length(list.files(bpcells_dir, all.files = TRUE, no.. = TRUE)) == 0L) {
  imported <- import_matrix_market(
    mtx_path = matrix_path,
    outdir = bpcells_dir,
    row_names = features,
    col_names = barcodes,
    row_major = FALSE,
    tmpdir = file.path(out_dir, "tmp"),
    sort_bytes = 4 * 1024^3
  )
  rm(imported)
  gc()
}
counts <- open_matrix_dir(bpcells_dir)

message("Creating Seurat object")
tnbc <- CreateSeuratObject(counts = counts, project = "set2_TNBC", meta.data = cell_metadata, min.cells = 3)
tnbc$paper_dataset_id <- "set2_TNBC"
tnbc$paper_label <- "TNBC"
tnbc[["percent.mt"]] <- PercentageFeatureSet(tnbc, pattern = "^MT-")

qc_before <- data.frame(
  cell_id = colnames(tnbc),
  sample_id_from_title = tnbc$sample_id_from_title,
  patient_from_title = tnbc$patient_from_title,
  tissue = tnbc$tissue,
  group = tnbc$group,
  nCount_RNA = tnbc$nCount_RNA,
  nFeature_RNA = tnbc$nFeature_RNA,
  percent.mt = tnbc$percent.mt,
  stringsAsFactors = FALSE
)
write.table(qc_before, file.path(tab_dir, "tnbc_full_qc_before_filter.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pre_filter_n <- ncol(tnbc)
tnbc <- subset(tnbc, subset = percent.mt <= 20)
post_filter_n <- ncol(tnbc)

qc_summary <- aggregate(
  cbind(nCount_RNA, nFeature_RNA, percent.mt) ~ sample_id_from_title + patient_from_title + tissue + group,
  data = tnbc@meta.data,
  FUN = median
)
qc_summary$cells_after_mt_filter <- as.integer(table(tnbc$sample_id_from_title)[qc_summary$sample_id_from_title])
write.table(qc_summary, file.path(tab_dir, "tnbc_full_qc_summary_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Normalizing and clustering full TNBC RNA data")
tnbc <- NormalizeData(tnbc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
tnbc <- FindVariableFeatures(tnbc, selection.method = "vst", nfeatures = 3000, verbose = TRUE)
tnbc <- ScaleData(tnbc, features = VariableFeatures(tnbc), verbose = TRUE)
tnbc <- RunPCA(tnbc, features = VariableFeatures(tnbc), npcs = npcs, verbose = TRUE)
tnbc <- RunHarmony(tnbc, group.by.vars = "sample_id_from_title", reduction.use = "pca", dims.use = 1:npcs, verbose = TRUE)
tnbc <- RunUMAP(tnbc, reduction = "harmony", dims = 1:npcs, n.neighbors = 30, min.dist = 0.3, verbose = TRUE)
tnbc <- FindNeighbors(tnbc, reduction = "harmony", dims = 1:npcs, verbose = TRUE)
tnbc <- FindClusters(tnbc, resolution = resolution, algorithm = 1, verbose = TRUE)
tnbc <- JoinLayers(tnbc, assay = "RNA")

cluster_by_sample <- as.data.frame(table(cluster = Idents(tnbc), sample_id_from_title = tnbc$sample_id_from_title), stringsAsFactors = FALSE)
write.table(cluster_by_sample, file.path(tab_dir, "tnbc_full_cluster_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cluster_by_tissue <- as.data.frame(table(cluster = Idents(tnbc), tissue = tnbc$tissue), stringsAsFactors = FALSE)
write.table(cluster_by_tissue, file.path(tab_dir, "tnbc_full_cluster_by_tissue.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cluster_by_group <- as.data.frame(table(cluster = Idents(tnbc), group = tnbc$group), stringsAsFactors = FALSE)
write.table(cluster_by_group, file.path(tab_dir, "tnbc_full_cluster_by_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Finding TNBC cluster markers")
marker_args <- list(
  object = tnbc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)
if (marker_max_cells > 0L) {
  marker_args$max.cells.per.ident <- marker_max_cells
}
markers <- do.call(FindAllMarkers, marker_args)
write.table(markers, file.path(tab_dir, "tnbc_full_cluster_markers_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

top_markers <- do.call(
  rbind,
  lapply(split(markers, markers$cluster), function(x) {
    x <- x[order(x$p_val_adj, -x$avg_log2FC), , drop = FALSE]
    head(x, 10)
  })
)
write.table(top_markers, file.path(tab_dir, "tnbc_full_cluster_top10_markers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

umap_df <- data.frame(
  cell_id = colnames(tnbc),
  UMAP_1 = Embeddings(tnbc, "umap")[, 1],
  UMAP_2 = Embeddings(tnbc, "umap")[, 2],
  cluster = as.character(Idents(tnbc)),
  sample_id_from_title = tnbc$sample_id_from_title,
  patient_from_title = tnbc$patient_from_title,
  tissue = tnbc$tissue,
  group = tnbc$group,
  nCount_RNA = tnbc$nCount_RNA,
  nFeature_RNA = tnbc$nFeature_RNA,
  percent.mt = tnbc$percent.mt,
  stringsAsFactors = FALSE
)
write.table(umap_df, file.path(tab_dir, "tnbc_full_umap_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_dim <- function(group_by, label = FALSE) {
  DimPlot(tnbc, reduction = "umap", group.by = group_by, label = label, repel = TRUE, raster = TRUE) +
    theme_classic(base_size = 10) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2)) +
    theme(
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = 7),
      plot.title = element_text(size = 11, face = "bold")
    )
}

p_cluster <- plot_dim("seurat_clusters", label = TRUE) + ggtitle("TNBC full RNA clusters")
p_tissue <- plot_dim("tissue", label = FALSE) + ggtitle("TNBC full RNA by tissue")
p_group <- plot_dim("group", label = FALSE) + ggtitle("TNBC full RNA by treatment group")
ggsave(file.path(fig_dir, "tnbc_full_umap_by_cluster.png"), p_cluster, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "tnbc_full_umap_by_cluster.pdf"), p_cluster, width = 8, height = 6)
ggsave(file.path(fig_dir, "tnbc_full_umap_by_tissue.png"), p_tissue, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "tnbc_full_umap_by_tissue.pdf"), p_tissue, width = 8, height = 6)
ggsave(file.path(fig_dir, "tnbc_full_umap_by_group.png"), p_group, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "tnbc_full_umap_by_group.pdf"), p_group, width = 8, height = 6)

marker_genes <- c(
  "PTPRC", "LYZ", "LST1", "CST3", "S100A8", "S100A9", "FCN1", "VCAN",
  "FCGR3A", "MS4A7", "CD68", "CD163", "MRC1", "APOE", "C1QA", "C1QB",
  "TREM2", "SPP1", "THBS1", "FCER1A", "CLEC10A", "CD1C", "CLEC9A",
  "XCR1", "LILRA4", "IL3RA", "GZMB", "TPSAB1", "TPSB2", "CPA3", "KIT",
  "MKI67", "TOP2A", "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "MS4A1",
  "CD79A", "IGHG3", "EPCAM", "KRT19", "PECAM1", "VWF", "DCN", "COL1A1"
)
marker_genes <- intersect(marker_genes, rownames(tnbc))
if (length(marker_genes) > 0L) {
  p_dot <- DotPlot(tnbc, features = marker_genes, group.by = "seurat_clusters") +
    coord_flip() +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      plot.title = element_text(size = 11, face = "bold")
    ) +
    ggtitle("TNBC canonical marker audit")
  ggsave(file.path(fig_dir, "tnbc_full_canonical_marker_dotplot.png"), p_dot, width = 9, height = 10, dpi = 300)
  ggsave(file.path(fig_dir, "tnbc_full_canonical_marker_dotplot.pdf"), p_dot, width = 9, height = 10)
}

run_summary <- data.frame(
  pre_filter_cells = pre_filter_n,
  post_filter_cells = post_filter_n,
  removed_percent_mt_gt_20 = pre_filter_n - post_filter_n,
  clusters = length(levels(Idents(tnbc))),
  variable_features = length(VariableFeatures(tnbc)),
  samples = length(unique(tnbc$sample_id_from_title)),
  patients = length(unique(tnbc$patient_from_title)),
  npcs = npcs,
  resolution = resolution,
  marker_max_cells_per_cluster = marker_max_cells,
  stringsAsFactors = FALSE
)
write.table(run_summary, file.path(tab_dir, "tnbc_full_run_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(tnbc, file.path(out_dir, "tnbc_full_marker_clustering_seurat.rds"))
message("Done")
