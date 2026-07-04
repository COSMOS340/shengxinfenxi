#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
})

request_dir <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = TRUE)
input_dir <- file.path(request_dir, "inputs")
output_dir <- file.path(request_dir, "outputs")
raw_dir <- file.path(output_dir, "raw")
sparse_dir <- file.path(output_dir, "sparse")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
file_dir <- file.path(output_dir, "files")
log_dir <- file.path(output_dir, "logs")
log_path <- file.path(log_dir, "run_fig5e_gse120575_full_seurat.log")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sparse_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  message(line)
  write(line, file = log_path, append = TRUE)
}

assert_columns <- function(data, expected, label) {
  missing <- setdiff(expected, names(data))
  if (length(missing) > 0L) {
    stop(label, " missing columns: ", paste(missing, collapse = ", "),
         ". Observed columns: ", paste(names(data), collapse = ", "))
  }
}

format_p_label <- function(p_value) {
  vapply(p_value, function(one_p) {
    if (!is.finite(one_p)) return("NA")
    if (one_p < 1e-4) return(formatC(one_p, format = "e", digits = 1))
    as.character(signif(one_p, 2))
  }, character(1))
}

converter <- file.path(request_dir, "scripts", "convert_gse120575_tpm_to_mtx.py")
mtx_path <- file.path(sparse_dir, "gse120575_tpm.mtx.gz")
genes_path <- file.path(sparse_dir, "genes.tsv")
cells_path <- file.path(sparse_dir, "cells.tsv")
cell_meta_path <- file.path(sparse_dir, "cell_metadata.tsv")

if (!file.exists(mtx_path) || !file.exists(genes_path) || !file.exists(cells_path) || !file.exists(cell_meta_path)) {
  log_msg("Sparse files missing; running Python converter.")
  status <- system2("python3", converter)
  if (!identical(status, 0L)) stop("Python converter failed with status: ", status)
}

log_msg("Reading sparse TPM matrix.")
expr <- readMM(mtx_path)
genes <- fread(genes_path, header = FALSE)[[1]]
cells <- fread(cells_path, header = FALSE)[[1]]
if (nrow(expr) != length(genes)) stop("Gene count mismatch in sparse matrix.")
if (ncol(expr) != length(cells)) stop("Cell count mismatch in sparse matrix.")
rownames(expr) <- make.unique(genes)
colnames(expr) <- cells
expr <- as(expr, "dgCMatrix")

cell_meta <- fread(cell_meta_path)
assert_columns(cell_meta, c("cell_id", "sample_id", "timepoint", "response_group"), "GSE120575 cell metadata")
cell_meta <- as.data.frame(cell_meta)
rownames(cell_meta) <- cell_meta$cell_id
cell_meta <- cell_meta[colnames(expr), , drop = FALSE]

log_msg("Creating Seurat object from TPM matrix.")
obj <- CreateSeuratObject(counts = expr, meta.data = cell_meta, project = "GSE120575")

log_msg("Setting data layer to TPM values without additional normalization.")
obj <- tryCatch(
  SetAssayData(obj, assay = "RNA", layer = "data", new.data = expr),
  error = function(e) SetAssayData(obj, assay = "RNA", slot = "data", new.data = expr)
)

log_msg("Running variable feature selection, scaling, PCA, neighbors, resolution 1.0 clustering.")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 1.0, algorithm = 1, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)

cluster_col <- "seurat_clusters"
marker_genes <- intersect(c("CD68", "LYZ", "LST1", "AIF1", "C1QC", "C1QA", "SPP1", "THBS1", "IL1B", "S100A8", "S100A9"), rownames(obj))
cluster_avg <- AverageExpression(obj, assays = "RNA", features = marker_genes, slot = "data", verbose = FALSE)$RNA
cluster_marker_audit <- as.data.table(t(as.matrix(cluster_avg)), keep.rownames = "seurat_cluster")
cluster_counts <- as.data.table(table(obj[[cluster_col]][, 1]))
setnames(cluster_counts, c("seurat_cluster", "cells"))
cluster_marker_audit <- merge(cluster_marker_audit, cluster_counts, by = "seurat_cluster", all.x = TRUE)
fwrite(cluster_marker_audit, file.path(table_dir, "fig5e_full_cluster_marker_audit.tsv"), sep = "\t", quote = FALSE, na = "NA")

required_clusters <- c("8", "12", "15")
observed_clusters <- sort(unique(as.character(obj[[cluster_col]][, 1])))
missing_clusters <- setdiff(required_clusters, observed_clusters)
if (length(missing_clusters) > 0L) {
  stop("Required paper clusters missing at resolution 1.0: ", paste(missing_clusters, collapse = ", "),
       ". Upload fig5e_full_cluster_marker_audit.tsv for review.")
}

log_msg("Subsetting paper-stated myeloid clusters 8, 12, and 15.")
myeloid_cells <- colnames(obj)[as.character(obj[[cluster_col]][, 1]) %in% required_clusters]
if (length(myeloid_cells) == 0L) stop("No cells found in clusters 8, 12, and 15.")
myeloid_obj <- subset(obj, cells = myeloid_cells)

log_msg("Reclustering myeloid subset at resolution 0.4.")
myeloid_obj <- FindVariableFeatures(myeloid_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
myeloid_obj <- ScaleData(myeloid_obj, features = VariableFeatures(myeloid_obj), verbose = FALSE)
myeloid_obj <- RunPCA(myeloid_obj, features = VariableFeatures(myeloid_obj), npcs = 30, verbose = FALSE)
myeloid_obj <- FindNeighbors(myeloid_obj, dims = 1:20, verbose = FALSE)
myeloid_obj <- FindClusters(myeloid_obj, resolution = 0.4, algorithm = 1, verbose = FALSE)
myeloid_obj <- RunUMAP(myeloid_obj, dims = 1:20, verbose = FALSE)
myeloid_obj$myeloid_subcluster <- as.character(myeloid_obj$seurat_clusters)

fig2_ref_path <- file.path(input_dir, "fig2_monocyte_macrophage_cluster_average_expression_common_features.tsv.gz")
fig2_annot_path <- file.path(input_dir, "fig2_monocyte_macrophage_cluster_subtype_annotation_current_scope.tsv")
fig2_ref_raw <- fread(fig2_ref_path, header = FALSE, check.names = FALSE)
ref_headers <- as.character(unlist(fig2_ref_raw[1]))
fig2_ref <- fig2_ref_raw[-1]
setnames(fig2_ref, ref_headers)
fig2_ref[, gene := as.character(gene)]
ref_cols <- setdiff(names(fig2_ref), "gene")
for (ref_col in ref_cols) fig2_ref[, (ref_col) := as.numeric(get(ref_col))]
fig2_annot <- fread(fig2_annot_path)
assert_columns(fig2_annot, c("fig2_cluster", "fig2_subtype_current_scope"), "Fig. 2 annotation")
fig2_annot[, fig2_cluster := as.character(fig2_cluster)]

query_avg <- AverageExpression(myeloid_obj, assays = "RNA", slot = "data", group.by = "myeloid_subcluster", verbose = FALSE)$RNA
common_genes <- intersect(rownames(query_avg), fig2_ref$gene)
if (length(common_genes) < 50L) stop("Too few common genes for subcluster annotation: ", length(common_genes))
query_common <- as.matrix(query_avg[common_genes, , drop = FALSE])
ref_common <- as.matrix(fig2_ref[match(common_genes, gene), ..ref_cols])
rownames(ref_common) <- common_genes
storage.mode(query_common) <- "numeric"
storage.mode(ref_common) <- "numeric"

gene_means <- rowMeans(query_common, na.rm = TRUE)
gene_sds <- apply(query_common, 1, sd, na.rm = TRUE)
gene_sds[!is.finite(gene_sds) | gene_sds == 0] <- 1
query_scaled <- sweep(sweep(query_common, 1, gene_means, "-"), 1, gene_sds, "/")
ref_scaled <- sweep(sweep(ref_common, 1, gene_means, "-"), 1, gene_sds, "/")
norm_cols <- function(mat) {
  norms <- sqrt(colSums(mat * mat, na.rm = TRUE))
  norms[!is.finite(norms) | norms == 0] <- 1
  sweep(mat, 2, norms, "/")
}
cor_mat <- crossprod(norm_cols(query_scaled), norm_cols(ref_scaled))
best_ref_index <- max.col(cor_mat, ties.method = "first")
subcluster_annotation <- data.table(
  myeloid_subcluster = rownames(cor_mat),
  fig2_cluster = as.character(as.integer(ref_cols[best_ref_index])),
  transfer_score = cor_mat[cbind(seq_len(nrow(cor_mat)), best_ref_index)]
)
subcluster_annotation <- merge(
  subcluster_annotation,
  fig2_annot[, .(fig2_cluster, fig2_subtype_current_scope)],
  by = "fig2_cluster",
  all.x = TRUE,
  sort = FALSE
)
subcluster_annotation[, fig5e_lineage := fifelse(
  fig2_subtype_current_scope %in% c("Resting C1QC+ TAMs", "Activated C1QC+ TAMs"),
  "C1QC+ TAMs",
  fifelse(
    fig2_subtype_current_scope %in% c("THBS1+ MDSCs", "SPP1+ TAMs"),
    "THBS1+ MDSCs + SPP1+ TAMs",
    "Other myeloid"
  )
)]

subcluster_counts <- as.data.table(table(myeloid_obj$myeloid_subcluster))
setnames(subcluster_counts, c("myeloid_subcluster", "cells"))
subcluster_annotation <- merge(subcluster_annotation, subcluster_counts, by = "myeloid_subcluster", all.x = TRUE)
fwrite(subcluster_annotation, file.path(table_dir, "fig5e_full_myeloid_subcluster_annotation.tsv"), sep = "\t", quote = FALSE, na = "NA")

annotation_map <- subcluster_annotation[, .(myeloid_subcluster, fig2_subtype_current_scope, fig5e_lineage)]
myeloid_meta <- as.data.table(myeloid_obj@meta.data, keep.rownames = "cell_id")
myeloid_meta <- merge(myeloid_meta, annotation_map, by = "myeloid_subcluster", all.x = TRUE)

sample_den <- myeloid_meta[, .(denominator_cells = .N), by = .(sample_id, timepoint, response_group)]
sample_num <- myeloid_meta[
  fig5e_lineage %in% c("C1QC+ TAMs", "THBS1+ MDSCs + SPP1+ TAMs"),
  .N,
  by = .(sample_id, fig5e_lineage)
]
sample_grid <- CJ(
  sample_id = sample_den$sample_id,
  fig5e_lineage = c("C1QC+ TAMs", "THBS1+ MDSCs + SPP1+ TAMs"),
  unique = TRUE
)
sample_grid <- merge(sample_grid, sample_den, by = "sample_id", all.x = TRUE)
sample_props <- merge(sample_grid, sample_num, by = c("sample_id", "fig5e_lineage"), all.x = TRUE)
sample_props[is.na(N), N := 0L]
sample_props[, proportion := N / denominator_cells]
sample_props[, timepoint := factor(timepoint, levels = c("Pretreatment", "Posttreatment"))]
sample_props[, response_group := factor(response_group, levels = c("R", "NR"))]
sample_props[, fig5e_lineage := factor(fig5e_lineage, levels = c("C1QC+ TAMs", "THBS1+ MDSCs + SPP1+ TAMs"))]
fwrite(sample_props, file.path(table_dir, "fig5e_full_sample_proportions.tsv"), sep = "\t", quote = FALSE, na = "NA")

stats_dt <- sample_props[, {
  r_values <- proportion[response_group == "R"]
  nr_values <- proportion[response_group == "NR"]
  data.table(
    n_r = length(r_values),
    n_nr = length(nr_values),
    p_value = suppressWarnings(wilcox.test(r_values, nr_values, exact = FALSE)$p.value),
    median_r = median(r_values, na.rm = TRUE),
    median_nr = median(nr_values, na.rm = TRUE),
    mean_r = mean(r_values, na.rm = TRUE),
    mean_nr = mean(nr_values, na.rm = TRUE)
  )
}, by = .(timepoint, fig5e_lineage)]
stats_dt[, p_label := format_p_label(p_value)]
fwrite(stats_dt, file.path(table_dir, "fig5e_full_wilcoxon.tsv"), sep = "\t", quote = FALSE, na = "NA")

lineage_display <- c(
  "C1QC+ TAMs" = "C1QC+ TAMs",
  "THBS1+ MDSCs + SPP1+ TAMs" = "THBS1+ MDSCs\n+ SPP1+ TAMs"
)
sample_props[, fig5e_lineage_display := lineage_display[as.character(fig5e_lineage)]]
sample_props[, panel_label := paste(timepoint, fig5e_lineage_display, sep = "\n")]
panel_levels <- c(
  "Pretreatment\nC1QC+ TAMs",
  "Pretreatment\nTHBS1+ MDSCs\n+ SPP1+ TAMs",
  "Posttreatment\nC1QC+ TAMs",
  "Posttreatment\nTHBS1+ MDSCs\n+ SPP1+ TAMs"
)
sample_props[, panel_label := factor(panel_label, levels = panel_levels)]
stats_dt[, fig5e_lineage_display := lineage_display[as.character(fig5e_lineage)]]
stats_dt[, panel_label := factor(paste(timepoint, fig5e_lineage_display, sep = "\n"), levels = panel_levels)]
stats_dt[, `:=`(group1 = "R", group2 = "NR", y.position = 0.92, label = p_label)]

plot_e <- ggplot(sample_props, aes(response_group, proportion, fill = response_group)) +
  geom_boxplot(width = 0.58, outlier.size = 0.65, linewidth = 0.25) +
  ggpubr::stat_pvalue_manual(
    stats_dt,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.22,
    size = 1.9
  ) +
  facet_wrap(~ panel_label, nrow = 1) +
  scale_fill_manual(values = c("R" = "#D76B64", "NR" = "#2AA6A3")) +
  scale_y_continuous(limits = c(0, 1.0), breaks = c(0, 0.3, 0.6, 0.9), expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = "GSE120575 melanoma", x = NULL, y = "Proportion") +
  theme_classic(base_size = 7) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 5.3, lineheight = 0.95),
    axis.text.x = element_text(size = 6.2),
    axis.text.y = element_text(size = 6.2),
    axis.title.y = element_text(size = 7),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
    legend.position = "none",
    panel.spacing.x = unit(0.45, "lines"),
    plot.margin = margin(4, 4, 4, 4)
  )
ggsave(file.path(figure_dir, "fig5e_full_seurat_proportions.png"), plot_e, width = 5.2, height = 1.7, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "fig5e_full_seurat_proportions.pdf"), plot_e, width = 5.2, height = 1.7, bg = "white")

log_msg("Saving myeloid reclustered Seurat object.")
saveRDS(myeloid_obj, file.path(file_dir, "gse120575_myeloid_reclustered_seurat.rds"))
log_msg("Fig5E full Seurat workflow complete.")
