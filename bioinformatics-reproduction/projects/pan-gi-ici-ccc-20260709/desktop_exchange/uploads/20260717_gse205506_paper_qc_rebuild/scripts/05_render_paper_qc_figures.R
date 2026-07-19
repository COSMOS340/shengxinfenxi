#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

options(stringsAsFactors = FALSE, future.globals.maxSize = 32 * 1024^3)
future::plan("sequential")
set.seed(340)

gse_root <- "F:/pan-gi-ici-ccc-20260709/GSE205506"
project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
output_dir <- file.path(project_dir, "desktop_exchange/uploads/20260717_gse205506_paper_qc_rebuild")
figure_dir <- file.path(output_dir, "figures")
analysis_dir <- file.path(gse_root, "r_analysis_paper_qc_20260717")
final_path <- file.path(analysis_dir, "gse205506_formal_post_mt_rpca_atlas.rds")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

log_step <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
  flush.console()
}

save_figure <- function(plot, stem, width, height) {
  ggsave(
    filename = file.path(figure_dir, paste0(stem, ".png")),
    plot = plot,
    device = ragg::agg_png,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggsave(
    filename = file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}

base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8),
    axis.text = element_text(color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

broad_palette <- c(
  "Epithelial" = "#E7298A",
  "T/I/NK" = "#1F78B4",
  "B" = "#6A3D9A",
  "Myeloid" = "#D95F02",
  "Endothelial" = "#1B9E77",
  "Fibroblast" = "#A6761D",
  "Hybrid_or_unresolved_review" = "#666666"
)
response_palette <- c("pCR" = "#009E73", "non-pCR" = "#D55E00")
timepoint_palette <- c("pre-treatment" = "#0072B2", "post-treatment" = "#D55E00")

if (!file.exists(final_path)) stop("Final RPCA atlas is missing", call. = FALSE)
log_step("Loading final atlas for all-cell figures")
atlas <- readRDS(final_path)
DefaultAssay(atlas) <- "RNA"
required_metadata <- c(
  "paper_final_broad_compartment", "table_s1_response",
  "derived_timepoint_from_exact_geo_treatment", "geo_accession"
)
if (!all(required_metadata %in% colnames(atlas[[]]))) stop("Final atlas plotting metadata is incomplete", call. = FALSE)

set.seed(340)
p_broad <- DimPlot(
  atlas,
  reduction = "umap.rpca.final",
  group.by = "paper_final_broad_compartment",
  cols = broad_palette,
  shuffle = TRUE,
  seed = 340,
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 1,
  label = TRUE,
  repel = TRUE,
  label.box = FALSE
) +
  labs(title = "GSE205506 final broad atlas", color = "Broad compartment") +
  base_theme +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
save_figure(p_broad, "gse205506_paper_qc_umap_broad", 10, 7.5)

p_response <- DimPlot(
  atlas,
  reduction = "umap.rpca.final",
  group.by = "table_s1_response",
  shuffle = TRUE,
  seed = 340,
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 1
) +
  labs(title = "GSE205506 response", color = "Table S1 response") +
  base_theme +
  scale_color_manual(values = response_palette, na.value = "#BDBDBD") +
  guides(color = guide_legend(override.aes = list(size = 3)))
save_figure(p_response, "gse205506_paper_qc_umap_response", 9, 7.5)

p_timepoint <- DimPlot(
  atlas,
  reduction = "umap.rpca.final",
  group.by = "derived_timepoint_from_exact_geo_treatment",
  shuffle = TRUE,
  seed = 340,
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 1
) +
  labs(title = "GSE205506 treatment timepoint", color = "Timepoint") +
  base_theme +
  scale_color_manual(values = timepoint_palette, na.value = "#BDBDBD") +
  guides(color = guide_legend(override.aes = list(size = 3)))
save_figure(p_timepoint, "gse205506_paper_qc_umap_timepoint", 9, 7.5)

sample_levels <- sort(unique(atlas$geo_accession))
sample_palette <- setNames(grDevices::hcl.colors(length(sample_levels), "Dark 3"), sample_levels)
p_sample <- DimPlot(
  atlas,
  reduction = "umap.rpca.final",
  group.by = "geo_accession",
  cols = sample_palette,
  shuffle = TRUE,
  seed = 340,
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 1
) +
  labs(title = "GSE205506 sample distribution", color = "GEO sample") +
  base_theme +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 4, byrow = TRUE)) +
  theme(legend.position = "bottom")
save_figure(p_sample, "gse205506_paper_qc_umap_sample", 13, 10)

broad_features <- list(
  `T/I/NK` = c("CD3D", "CD3E", "TRAC", "TRBC1"),
  B = c("CD79A", "CD79B", "MS4A1", "TNFRSF17", "MZB1"),
  Myeloid = c("CD14", "CD68"),
  Epithelial = c("EPCAM", "CD24"),
  Fibroblast = c("COL1A2", "COL3A1", "MYH11", "ACTA2"),
  Endothelial = c("VWF", "PECAM1")
)
p_broad_dot <- DotPlot(
  atlas,
  features = broad_features,
  group.by = "paper_final_broad_compartment",
  assay = "RNA",
  cols = c("#F7F7F7", "#D55E00"),
  dot.scale = 6
) +
  RotatedAxis() +
  labs(title = "GSE205506 broad-compartment markers", x = NULL, y = "Broad compartment") +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
save_figure(p_broad_dot, "gse205506_paper_qc_marker_dotplot_broad", 13, 6.5)

rm(atlas)
invisible(gc())

log_step("Rendering QC before/after comparison")
cell_audit_path <- file.path(output_dir, "gse205506_paper_qc_cell_audit.tsv.gz")
cell_audit <- fread(cell_audit_path)
required_audit_columns <- c(
  "nFeature_RNA_before_gene_filter",
  "nCount_RNA_before_gene_filter",
  "percent_mt_before_gene_filter",
  "formal_final_inclusion"
)
if (!all(required_audit_columns %in% names(cell_audit))) stop("Cell audit plotting columns are incomplete", call. = FALSE)
set.seed(340)
before_index <- sample.int(nrow(cell_audit), min(80000L, nrow(cell_audit)))
after_pool <- which(cell_audit$formal_final_inclusion %in% TRUE)
after_index <- sample(after_pool, min(80000L, length(after_pool)))
qc_plot_data <- rbindlist(list(
  cell_audit[before_index, .(
    stage = "All raw barcodes",
    nFeature_RNA = nFeature_RNA_before_gene_filter,
    nCount_RNA = nCount_RNA_before_gene_filter,
    percent.mt = percent_mt_before_gene_filter
  )],
  cell_audit[after_index, .(
    stage = "Final paper-QC atlas",
    nFeature_RNA = nFeature_RNA_before_gene_filter,
    nCount_RNA = nCount_RNA_before_gene_filter,
    percent.mt = percent_mt_before_gene_filter
  )]
))
qc_long <- melt(
  qc_plot_data,
  id.vars = "stage",
  measure.vars = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  variable.name = "metric",
  value.name = "value"
)
qc_long[, metric := factor(
  metric,
  levels = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  labels = c("Detected genes", "UMI counts", "Mitochondrial percentage")
)]
qc_palette <- c("All raw barcodes" = "#BDBDBD", "Final paper-QC atlas" = "#0072B2")
p_qc <- ggplot(qc_long, aes(x = stage, y = value, fill = stage)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = qc_palette) +
  labs(title = "GSE205506 QC distributions before and after the formal route", x = NULL, y = NULL) +
  base_theme +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "none")
save_figure(p_qc, "gse205506_paper_qc_before_after", 13, 5.5)
rm(cell_audit, qc_plot_data, qc_long)
invisible(gc())

compartment_spec <- data.table(
  compartment = c("T/I/NK", "B", "Myeloid", "Endothelial", "Fibroblast"),
  slug = c("T_I_NK", "B", "Myeloid", "Endothelial", "Fibroblast")
)
for (index in seq_len(nrow(compartment_spec))) {
  compartment_name <- compartment_spec$compartment[[index]]
  slug <- compartment_spec$slug[[index]]
  log_step(paste("Rendering figures for", compartment_name))
  object_path <- file.path(analysis_dir, paste0("gse205506_compartment_", slug, ".rds"))
  if (!file.exists(object_path)) stop("Compartment object is missing: ", object_path, call. = FALSE)
  current <- readRDS(object_path)
  DefaultAssay(current) <- "RNA"

  cluster_levels <- sort(unique(current$compartment_working_cluster))
  cluster_palette <- setNames(grDevices::hcl.colors(length(cluster_levels), "Dark 3"), cluster_levels)
  reviewed_levels <- sort(unique(current$compartment_reviewed_label))
  reviewed_palette <- setNames(grDevices::hcl.colors(length(reviewed_levels), "Dark 3"), reviewed_levels)
  reviewed_palette[reviewed_levels == "Unresolved_review"] <- "#7A5195"

  p_cluster <- DimPlot(
    current,
    reduction = "umap.compartment",
    group.by = "compartment_working_cluster",
    cols = cluster_palette,
    shuffle = TRUE,
    seed = 340,
    raster = TRUE,
    raster.dpi = c(300, 300),
    pt.size = 0.8,
    label = TRUE,
    repel = TRUE,
    label.box = FALSE
  ) +
    labs(title = paste(compartment_name, "working clusters"), color = "Cluster") +
    base_theme +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))

  p_reviewed <- DimPlot(
    current,
    reduction = "umap.compartment",
    group.by = "compartment_reviewed_label",
    cols = reviewed_palette,
    shuffle = TRUE,
    seed = 340,
    raster = TRUE,
    raster.dpi = c(300, 300),
    pt.size = 0.8
  ) +
    labs(title = paste(compartment_name, "author-marker review"), color = "Reviewed label") +
    base_theme +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

  p_compartment <- p_cluster + p_reviewed + plot_layout(widths = c(1, 1.15))
  save_figure(p_compartment, paste0("gse205506_paper_qc_umap_", slug), 17, 7.5)

  marker_path <- file.path(output_dir, paste0("gse205506_paper_qc_", slug, "_top50_markers.tsv"))
  marker_table <- fread(marker_path)
  marker_table <- marker_table[!grepl("^MT-|^RPS[0-9]|^RPL[0-9]", gene)]
  marker_table[, cluster_numeric := as.integer(cluster)]
  setorder(marker_table, cluster_numeric, rank)
  dot_features <- marker_table[, head(gene, 2L), by = cluster]$V1
  dot_features <- unique(dot_features)
  p_dot <- DotPlot(
    current,
    features = dot_features,
    group.by = "compartment_working_cluster",
    assay = "RNA",
    cols = c("#F7F7F7", "#D55E00"),
    dot.scale = 5
  ) +
    RotatedAxis() +
    labs(title = paste(compartment_name, "cluster markers"), x = NULL, y = "Working cluster") +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7))
  dot_width <- max(12, min(28, 4 + 0.3 * length(dot_features)))
  dot_height <- max(5.5, min(12, 3.5 + 0.3 * length(cluster_levels)))
  save_figure(p_dot, paste0("gse205506_paper_qc_marker_dotplot_", slug), dot_width, dot_height)

  rm(current, marker_table, p_cluster, p_reviewed, p_compartment, p_dot)
  invisible(gc())
}

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_gse205506_paper_qc_figures.txt"))
log_step("All required paper-QC figures rendered")
