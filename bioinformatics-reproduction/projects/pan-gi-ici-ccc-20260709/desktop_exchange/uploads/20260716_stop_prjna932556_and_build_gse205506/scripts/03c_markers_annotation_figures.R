#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(presto)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 03c_markers_annotation_figures.R <gse205506_root> <output_dir> <seed>",
    call. = FALSE
  )
}

gse_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
seed <- as.integer(args[[3]])
analysis_dir <- file.path(gse_root, "r_analysis")
input_rds <- file.path(analysis_dir, "gse205506_seurat_integrated_clusters.rds")
final_rds <- file.path(analysis_dir, "gse205506_seurat_final_annotated.rds")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(
  input_rds,
  file.path(output_dir, "gse205506_qc_filter_summary.tsv"),
  file.path(output_dir, "gse205506_author_marker_reference_top10.tsv")
)
if (any(!file.exists(required_paths))) {
  stop("Missing required files: ", paste(required_paths[!file.exists(required_paths)], collapse = "; "), call. = FALSE)
}

write_tsv <- function(x, filename) {
  fwrite(x, file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

save_figure <- function(plot, stem, width, height) {
  plot <- plot + theme(plot.margin = margin(10, 14, 10, 10))
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, device = cairo_pdf, bg = "white", limitsize = FALSE
  )
}

object <- readRDS(input_rds)
if (ncol(object) != 238934L || length(unique(object$working_cluster)) != 28L) {
  stop("Integrated object does not match the verified 238934-cell, 28-cluster checkpoint", call. = FALSE)
}
DefaultAssay(object) <- "RNA"
Idents(object) <- "working_cluster"

message("Running full-cell presto Wilcoxon/AUC marker analysis")
marker_all <- as.data.table(wilcoxauc(
  object,
  group_by = "working_cluster",
  assay = "data",
  seurat_assay = "RNA"
))
setnames(marker_all, c("feature", "group"), c("gene", "cluster"))
marker_all[, cluster := as.character(cluster)]
marker_all[, positive_direction := auc > 0.5 & logFC > 0]
marker_all[, passes_primary_filter := padj <= 0.05 & positive_direction == TRUE & pct_in >= 5]
marker_filtered <- marker_all[passes_primary_filter == TRUE]
marker_rank_pool <- copy(marker_all)
marker_rank_pool[, cluster_numeric := as.integer(cluster)]
setorder(marker_rank_pool, cluster_numeric, -passes_primary_filter, -positive_direction, -auc, -logFC, -pct_in, gene)
marker_rank_pool[, rank := seq_len(.N), by = cluster]
marker_rank_pool[, cluster_numeric := NULL]
top50 <- marker_rank_pool[rank <= 50L]
top10 <- marker_rank_pool[rank <= 10L]
if (nrow(top50) != 28L * 50L || nrow(top10) != 28L * 10L) {
  stop("Top marker tables do not contain exactly 50 and 10 rows for each of 28 working clusters", call. = FALSE)
}
write_tsv(top10, "gse205506_cluster_top10_markers.tsv")
write_tsv(top50, "gse205506_cluster_top50_markers.tsv")

signature_genes <- list(
  T_NK = c("CD3D", "CD3E", "TRAC", "IL7R", "LTB", "NKG7", "GNLY", "KLRD1", "GZMB"),
  B_Plasma = c("CD79A", "MS4A1", "CD37", "CD74", "CD79B", "MZB1", "JCHAIN", "XBP1", "SDC1"),
  Myeloid = c("LST1", "TYROBP", "FCER1G", "CTSS", "LILRB1", "AIF1", "C1QA", "C1QB", "S100A8"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "MUC1", "CEACAM5", "TACSTD2"),
  Endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "CLDN5", "RAMP2", "PLVAP", "ENG"),
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A1", "LUM", "COL14A1", "PDGFRA"),
  Mast = c("TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2", "HDC", "HPGDS")
)
cycling_genes <- c("MKI67", "TOP2A", "UBE2C", "STMN1", "TYMS", "CENPF")
available_signatures <- lapply(signature_genes, intersect, y = rownames(object))
if (any(lengths(available_signatures) < 4L)) {
  stop("Fewer than four canonical genes were found for a broad cell-type signature", call. = FALSE)
}

cluster_levels <- as.character(sort(as.integer(as.character(unique(object$working_cluster)))))
canonical_genes <- unique(c(unlist(available_signatures), intersect(cycling_genes, rownames(object))))
rna_data <- GetAssayData(object, assay = "RNA", layer = "data")[canonical_genes, , drop = FALSE]
cluster_cells <- split(seq_len(ncol(object)), factor(object$working_cluster, levels = cluster_levels))
average_expression <- vapply(cluster_cells, function(index) {
  Matrix::rowMeans(rna_data[, index, drop = FALSE])
}, numeric(length(canonical_genes)))
rownames(average_expression) <- canonical_genes
colnames(average_expression) <- cluster_levels
gene_z <- t(scale(t(average_expression), center = TRUE, scale = TRUE))
gene_z[!is.finite(gene_z)] <- 0

score_rows <- rbindlist(lapply(names(available_signatures), function(signature) {
  genes <- available_signatures[[signature]]
  data.table(
    cluster = cluster_levels,
    signature = signature,
    signature_score = colMeans(gene_z[genes, , drop = FALSE]),
    mean_normalized_expression = colMeans(average_expression[genes, , drop = FALSE]),
    genes_used = paste(genes, collapse = ";")
  )
}))
score_rows[, cluster_numeric := as.integer(cluster)]
setorder(score_rows, cluster_numeric, -signature_score, signature)
score_rows[, signature_rank := seq_len(.N), by = cluster]
score_rows[, cluster_numeric := NULL]

annotation <- score_rows[signature_rank <= 2L, .(
  top_signature = signature[[1]],
  top_signature_score = signature_score[[1]],
  top_signature_mean_expression = mean_normalized_expression[[1]],
  second_signature = signature[[2]],
  second_signature_score = signature_score[[2]],
  score_margin = signature_score[[1]] - signature_score[[2]]
), by = cluster]
annotation[, cells := as.integer(table(factor(object$working_cluster, levels = cluster_levels)))]
annotation[, top10_markers := vapply(cluster, function(value) {
  paste(top10[cluster == value, gene], collapse = ";")
}, character(1))]
annotation[, broad_cell_type := top_signature]
annotation[, cycling_marker_count := vapply(cluster, function(value) {
  sum(intersect(top50[cluster == value, gene], cycling_genes) %in% cycling_genes)
}, integer(1))]
annotation[, cell_state := ifelse(cycling_marker_count >= 2L, "Cycling", "Non-cycling")]
annotation[, annotation_method := "highest cluster-standardized canonical signature; checked against cluster Top10/Top50 markers"]

broad_map <- setNames(annotation$broad_cell_type, annotation$cluster)
state_map <- setNames(annotation$cell_state, annotation$cluster)
object$cell_type_working <- unname(broad_map[as.character(object$working_cluster)])
object$cell_state_working <- unname(state_map[as.character(object$working_cluster)])
if (anyNA(object$cell_type_working) || anyNA(object$cell_state_working)) {
  missing_cluster_values <- unique(as.character(object$working_cluster)[is.na(object$cell_type_working)])
  stop(
    "Working annotation did not map to every cell; exact missing cluster values: ",
    paste(missing_cluster_values, collapse = ","),
    call. = FALSE
  )
}
write_tsv(annotation, "gse205506_annotation_evidence.tsv")
write_tsv(score_rows, "gse205506_annotation_signature_scores.tsv")

author_top10 <- fread(
  file.path(output_dir, "gse205506_author_marker_reference_top10.tsv"),
  sep = "\t", header = TRUE, na.strings = NULL
)
author_sets <- author_top10[, .(author_genes = list(unique(gene))), by = .(compartment_sheet, author_cluster = cluster)]
overlap_rows <- rbindlist(lapply(cluster_levels, function(cluster_value) {
  our_genes <- unique(top50[cluster == cluster_value, gene])
  author_sets[, .(
    cluster = cluster_value,
    compartment_sheet,
    author_cluster,
    our_top50_genes = length(our_genes),
    author_top10_genes = length(author_genes[[1]]),
    overlap_count = length(intersect(our_genes, author_genes[[1]])),
    overlap_genes = paste(intersect(our_genes, author_genes[[1]]), collapse = ";"),
    jaccard = length(intersect(our_genes, author_genes[[1]])) / length(union(our_genes, author_genes[[1]]))
  ), by = seq_len(nrow(author_sets))][, seq_len := NULL]
}))
overlap_rows[, cluster_numeric := as.integer(cluster)]
setorder(overlap_rows, cluster_numeric, -overlap_count, -jaccard, compartment_sheet, author_cluster)
overlap_rows[, overlap_rank := seq_len(.N), by = cluster]
overlap_rows[, cluster_numeric := NULL]
write_tsv(overlap_rows, "gse205506_author_marker_overlap_all.tsv")
write_tsv(overlap_rows[overlap_rank <= 3L], "gse205506_author_marker_overlap_top3.tsv")

metadata <- as.data.table(object@meta.data, keep.rownames = "cell_barcode")
patient_counts <- metadata[, .N, by = .(
  geo_subject,
  geo_accession,
  geo_genotype,
  geo_treatment,
  derived_timepoint_from_exact_geo_treatment,
  table_s1_response,
  working_cluster,
  cell_type_working,
  cell_state_working
)]
setnames(patient_counts, "N", "cells")
patient_counts[, cluster_numeric := as.integer(as.character(working_cluster))]
setorder(patient_counts, geo_subject, geo_accession, cluster_numeric)
patient_counts[, cluster_numeric := NULL]
write_tsv(patient_counts, "gse205506_patient_sample_cluster_cell_counts.tsv")

sample_retention <- fread(
  file.path(output_dir, "gse205506_qc_filter_summary.tsv"),
  sep = "\t", header = TRUE, na.strings = NULL
)
sample_retention[, geo_accession := factor(geo_accession, levels = geo_accession)]
p_qc <- ggplot(sample_retention, aes(geo_accession, final_retention_fraction, fill = timepoint)) +
  geom_col(width = 0.82) +
  scale_fill_manual(values = c("pre-treatment" = "#0072B2", "post-treatment" = "#D55E00")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Final singlet retention", fill = "Timepoint", title = "Per-sample QC retention") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
save_figure(p_qc, "gse205506_qc_sample_retention", 12, 6)

cluster_palette <- setNames(grDevices::hcl.colors(length(cluster_levels), "Dynamic"), cluster_levels)
broad_palette <- c(
  "T_NK" = "#0072B2",
  "B_Plasma" = "#E69F00",
  "Myeloid" = "#009E73",
  "Epithelial" = "#D55E00",
  "Endothelial" = "#56B4E9",
  "Fibroblast" = "#CC79A7",
  "Mast" = "#6F4E7C"
)
timepoint_palette <- c("pre-treatment" = "#0072B2", "post-treatment" = "#D55E00")
response_palette <- c("pCR" = "#009E73", "non-pCR" = "#CC79A7")

p_cluster <- DimPlot(
  object, reduction = "umap", group.by = "working_cluster",
  label = TRUE, repel = TRUE, label.size = 3.5, raster = TRUE,
  raster.dpi = c(300, 300), pt.size = 0.05, cols = cluster_palette
) + labs(title = "GSE205506 working clusters (resolution 0.6)", color = "Cluster") +
  theme_classic(base_size = 11) + theme(aspect.ratio = 1)
save_figure(p_cluster, "gse205506_umap_working_clusters", 9, 8)

p_cell_type <- DimPlot(
  object, reduction = "umap", group.by = "cell_type_working",
  label = TRUE, repel = TRUE, label.size = 4, raster = TRUE,
  raster.dpi = c(300, 300), pt.size = 0.05, cols = broad_palette
) + labs(title = "GSE205506 broad cell-type working annotation", color = "Cell type") +
  theme_classic(base_size = 11) + theme(aspect.ratio = 1)
save_figure(p_cell_type, "gse205506_umap_broad_cell_type", 9, 8)

p_timepoint <- DimPlot(
  object, reduction = "umap", group.by = "derived_timepoint_from_exact_geo_treatment",
  raster = TRUE, raster.dpi = c(300, 300), pt.size = 0.05, cols = timepoint_palette
) + labs(title = "Treatment timepoint", color = "Timepoint") +
  theme_classic(base_size = 11) + theme(aspect.ratio = 1)
save_figure(p_timepoint, "gse205506_umap_timepoint", 9, 8)

p_response <- DimPlot(
  object, reduction = "umap", group.by = "table_s1_response",
  raster = TRUE, raster.dpi = c(300, 300), pt.size = 0.05, cols = response_palette
) + labs(title = "Table S1 pathological response", color = "Response") +
  theme_classic(base_size = 11) + theme(aspect.ratio = 1)
save_figure(p_response, "gse205506_umap_response", 9, 8)

dot_features <- unique(unlist(signature_genes))
dot_features <- intersect(dot_features, rownames(object))
p_dot <- DotPlot(
  object,
  features = dot_features,
  group.by = "cell_type_working",
  assay = "RNA",
  dot.scale = 6
) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  labs(x = NULL, y = NULL, title = "Canonical marker evidence", color = "Scaled mean", size = "Percent expressed") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8))
save_figure(p_dot, "gse205506_marker_dotplot_broad_cell_type", 17, 6.5)

object@misc$marker_method <- list(
  package = paste0("presto ", as.character(packageVersion("presto"))),
  method = "full-cell Wilcoxon/AUC for each working cluster versus all other cells",
  marker_filter = "BH adjusted p <= 0.05, AUC > 0.5, positive presto logFC, pct_in >= 5",
  ranking = "primary-filter markers first, then positive-direction markers, then remaining genes by AUC descending"
)
object@misc$annotation_method <- "Working broad labels from cluster-standardized canonical signatures with Top10/Top50 marker evidence"
saveRDS(object, final_rds, compress = FALSE)

final_summary <- data.table(
  object_path = normalizePath(final_rds, winslash = "/", mustWork = TRUE),
  object_size_bytes = file.info(final_rds)$size,
  features = nrow(object),
  cells = ncol(object),
  working_clusters = uniqueN(object$working_cluster),
  broad_cell_types = uniqueN(object$cell_type_working),
  marker_rows_passing_filter = nrow(marker_filtered),
  top10_rows = nrow(top10),
  top50_rows = nrow(top50)
)
write_tsv(final_summary, "gse205506_final_object_summary.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_markers_annotation_figures.txt"))
message("Marker, annotation, figure, and final object workflow completed: ", final_rds)
