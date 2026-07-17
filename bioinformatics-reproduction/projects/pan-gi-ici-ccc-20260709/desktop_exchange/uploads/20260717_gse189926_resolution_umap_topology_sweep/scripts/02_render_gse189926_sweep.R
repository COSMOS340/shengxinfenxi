#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ggrastr)
  library(scales)
})

set.seed(340)
project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
output_dir <- file.path(project_dir, "desktop_exchange/uploads/20260717_gse189926_resolution_umap_topology_sweep")
object_path <- file.path(project_dir, "04_objects/20260717_gse189926_resolution_umap_topology_sweep/gse189926_resolution_umap_topology_sweep.rds")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

save_figure <- function(plot, stem, width, height) {
  plot <- plot + plot_annotation(theme = theme(plot.margin = margin(8, 10, 8, 8)))
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white", limitsize = FALSE)
}

if (!file.exists(object_path)) stop("Sweep object is missing", call. = FALSE)
object <- readRDS(object_path)
route_parameters <- fread(file.path(output_dir, "gse189926_umap_route_parameters.tsv"))
route_parameters[, `:=`(
  thread_count = 1L,
  thread_count_source = "Seurat 5.4.0 RunUMAP calls uwot with future::nbrOfWorkers(); this run used the default one-worker plan; the attempted direct n_threads argument was ignored and recorded in stderr"
)]
write_tsv(route_parameters, "gse189926_umap_route_parameters.tsv")
object@misc$gse189926_sweep$routes <- route_parameters
saveRDS(object, object_path, compress = FALSE)
object_summary <- fread(file.path(output_dir, "gse189926_sweep_object_summary.tsv"))
object_summary[, `:=`(
  size_bytes = as.numeric(file.info(object_path)$size),
  sha256 = digest::digest(file = object_path, algo = "sha256", serialize = FALSE)
)]
write_tsv(object_summary, "gse189926_sweep_object_summary.tsv")
resolution_selection <- fread(file.path(output_dir, "gse189926_resolution_selection_evidence.tsv"))
umap_selection <- fread(file.path(output_dir, "gse189926_umap_selection_evidence.tsv"))
linearity_values <- fread(file.path(output_dir, "gse189926_umap_local_linearity_cell_values.tsv.gz"))
filament_review <- fread(file.path(output_dir, "gse189926_filament_region_review.tsv"))
resolution_fields <- setNames(
  c("sweep_res_0_2", "sweep_res_0_3", "sweep_res_0_4", "sweep_res_0_5", "sweep_res_0_6", "sweep_res_0_8", "sweep_res_1_0"),
  c("0.2", "0.3", "0.4", "0.5", "0.6", "0.8", "1.0")
)

metadata <- as.data.table(object[[]], keep.rownames = "cell_id")
set.seed(340)
plot_order <- sample(seq_len(nrow(metadata)))
metadata <- metadata[plot_order]

broad_palette <- c(
  "Granulocyte_like" = "#E69F00",
  "T_cells" = "#0072B2",
  "CD8_T_NK" = "#56B4E9",
  "Plasma_cells" = "#CC79A7",
  "Inflammatory_myeloid" = "#D55E00",
  "Macrophages" = "#009E73",
  "Monocytes" = "#006D2C",
  "Mast_cells" = "#7A5195",
  "Epithelial_like" = "#B79F00",
  "Endothelial_like" = "#17BECF",
  "Stromal_perivascular_like" = "#8C564B",
  "Dendritic_cells" = "#E15759",
  "Unresolved" = "#666666"
)
broad_display_map <- c(
  "Granulocyte_like" = "Granulocyte",
  "T_cells" = "T",
  "CD8_T_NK" = "CD8/NK",
  "Plasma_cells" = "Plasma",
  "Inflammatory_myeloid" = "Inflamm. myeloid",
  "Macrophages" = "Macrophage",
  "Monocytes" = "Monocyte",
  "Mast_cells" = "Mast",
  "Epithelial_like" = "Epithelial",
  "Endothelial_like" = "Endothelial",
  "Stromal_perivascular_like" = "Stromal/perivascular",
  "Dendritic_cells" = "DC",
  "Unresolved" = "Unresolved"
)
broad_display_palette <- setNames(unname(broad_palette[names(broad_display_map)]), unname(broad_display_map))
outcome_palette <- c("PD" = "#D55E00", "SD" = "#0072B2", "PR" = "#009E73", "NE" = "#7A5195")
timepoint_palette <- c(
  "pre-treatment" = "#0072B2",
  "post-treatment" = "#D55E00",
  "post-treatment1" = "#D55E00",
  "post-treatment2" = "#009E73"
)

dynamic_palette <- function(values, chroma = 80, luminance = 55) {
  levels <- sort(unique(as.character(values)))
  setNames(grDevices::hcl(h = seq(15, 375, length.out = length(levels) + 1L)[seq_along(levels)], c = chroma, l = luminance), levels)
}

embedding_data <- function(reduction_name) {
  coords <- Embeddings(object, reduction = reduction_name)
  data.table(
    cell_id = metadata$cell_id,
    UMAP_1 = coords[metadata$cell_id, 1L],
    UMAP_2 = coords[metadata$cell_id, 2L]
  )
}

discrete_panel <- function(data, values, title, palette, show_labels = FALSE, show_legend = TRUE) {
  plot_data <- copy(data)
  plot_data[, plot_value := as.character(values)]
  p <- ggplot(plot_data, aes(UMAP_1, UMAP_2, color = plot_value)) +
    ggrastr::geom_point_rast(size = 0.18, alpha = 0.90, stroke = 0, raster.dpi = 300) +
    scale_color_manual(values = palette, na.value = "#BDBDBD", drop = FALSE) +
    coord_equal() +
    labs(title = title, color = NULL) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 6.5),
      legend.key.height = grid::unit(3, "mm"),
      plot.margin = margin(4, 4, 4, 4)
    )
  if (show_labels) {
    labels <- plot_data[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = plot_value]
    p <- p + ggrepel::geom_text_repel(
      data = labels,
      aes(label = plot_value),
      size = 2.4,
      color = "black",
      seed = 340,
      box.padding = 0.15,
      point.padding = 0.05,
      min.segment.length = 0,
      max.overlaps = Inf,
      show.legend = FALSE
    )
  }
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

continuous_panel <- function(data, values, title, limits, option = "C") {
  plot_data <- copy(data)
  plot_data[, plot_value := values]
  ggplot(plot_data, aes(UMAP_1, UMAP_2, color = plot_value)) +
    ggrastr::geom_point_rast(size = 0.18, alpha = 0.90, stroke = 0, raster.dpi = 300) +
    scale_color_viridis_c(option = option, limits = limits, oob = scales::squish, na.value = "#BDBDBD") +
    coord_equal() +
    labs(title = title, color = NULL) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 6.5),
      plot.margin = margin(4, 4, 4, 4)
    )
}

baseline_data <- embedding_data("umap.harmony")
resolution_panels <- lapply(names(resolution_fields), function(resolution_name) {
  field <- resolution_fields[[resolution_name]]
  values <- metadata[[field]]
  discrete_panel(
    baseline_data,
    values,
    paste0("Resolution ", resolution_name, " (fixed baseline UMAP)"),
    dynamic_palette(values),
    show_labels = TRUE,
    show_legend = FALSE
  )
})
resolution_figure <- wrap_plots(resolution_panels, ncol = 3) +
  plot_annotation(title = "GSE189926 fixed-graph cluster resolution sweep on one unchanged UMAP")
save_figure(resolution_figure, "gse189926_fixed_coordinate_seven_resolution_comparison", 18, 16)

route_data <- setNames(lapply(route_parameters$reduction_name, embedding_data), route_parameters$route)

broad_panels <- lapply(route_parameters$route, function(route) {
  display_values <- unname(broad_display_map[as.character(metadata$broad_label)])
  discrete_panel(route_data[[route]], display_values, route, broad_display_palette, show_labels = TRUE, show_legend = TRUE)
})
broad_figure <- wrap_plots(broad_panels, ncol = 3, guides = "collect") +
  plot_annotation(title = "GSE189926 seven fixed-Harmony UMAP routes: broad labels") &
  theme(legend.position = "bottom")
save_figure(broad_figure, "gse189926_broad_label_seven_route_umap_comparison", 18, 16)

fixed_cluster_values <- metadata$cluster_harmony_r0_5
fixed_cluster_palette <- dynamic_palette(fixed_cluster_values)
cluster_panels <- lapply(route_parameters$route, function(route) {
  discrete_panel(route_data[[route]], fixed_cluster_values, route, fixed_cluster_palette, show_labels = TRUE, show_legend = FALSE)
})
cluster_figure <- wrap_plots(cluster_panels, ncol = 3) +
  plot_annotation(title = "GSE189926 seven fixed-Harmony UMAP routes: fixed resolution 0.5 labels")
save_figure(cluster_figure, "gse189926_fixed_cluster_seven_route_umap_comparison", 18, 16)

sample_palette <- dynamic_palette(metadata$sample_accession, 85, 50)
sample_panels <- lapply(route_parameters$route, function(route) {
  discrete_panel(route_data[[route]], metadata$sample_accession, route, sample_palette, show_labels = FALSE, show_legend = TRUE)
})
sample_figure <- wrap_plots(sample_panels, ncol = 4, guides = "collect") +
  plot_annotation(title = "GSE189926 seven UMAP routes: sample accession") &
  theme(legend.position = "bottom") &
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.2, alpha = 1)))
save_figure(sample_figure, "gse189926_sample_seven_route_umap_comparison", 20, 13)

patient_palette <- dynamic_palette(metadata$patient_id, 85, 50)
patient_panels <- lapply(route_parameters$route, function(route) {
  discrete_panel(route_data[[route]], metadata$patient_id, route, patient_palette, show_labels = FALSE, show_legend = TRUE)
})
patient_figure <- wrap_plots(patient_panels, ncol = 4, guides = "collect") +
  plot_annotation(title = "GSE189926 seven UMAP routes: patient") &
  theme(legend.position = "bottom") &
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.2, alpha = 1)))
save_figure(patient_figure, "gse189926_patient_seven_route_umap_comparison", 20, 13)

outcome_values <- metadata[["characteristics_ch1::outcome"]]
outcome_panels <- lapply(route_parameters$route, function(route) {
  discrete_panel(route_data[[route]], outcome_values, route, outcome_palette, show_labels = FALSE, show_legend = TRUE)
})
outcome_figure <- wrap_plots(outcome_panels, ncol = 4, guides = "collect") +
  plot_annotation(title = "GSE189926 seven UMAP routes: raw outcome") &
  theme(legend.position = "bottom")
save_figure(outcome_figure, "gse189926_raw_outcome_seven_route_umap_comparison", 20, 11)

timepoint_panels <- lapply(route_parameters$route, function(route) {
  discrete_panel(route_data[[route]], metadata$timepoint, route, timepoint_palette, show_labels = FALSE, show_legend = TRUE)
})
timepoint_figure <- wrap_plots(timepoint_panels, ncol = 4, guides = "collect") +
  plot_annotation(title = "GSE189926 seven UMAP routes: timepoint") &
  theme(legend.position = "bottom")
save_figure(timepoint_figure, "gse189926_timepoint_seven_route_umap_comparison", 20, 11)

qc_metrics <- c("nFeature_RNA", "nCount_RNA", "pct_mt", "pct_ribo", "pct_hb")
for (metric in qc_metrics) {
  values <- metadata[[metric]]
  limits <- as.numeric(quantile(values, c(0.01, 0.99), na.rm = TRUE, type = 8))
  panels <- lapply(route_parameters$route, function(route) {
    continuous_panel(route_data[[route]], values, route, limits)
  })
  figure <- wrap_plots(panels, ncol = 4, guides = "collect") +
    plot_annotation(title = paste0("GSE189926 seven UMAP routes: ", metric)) &
    theme(legend.position = "bottom")
  save_figure(figure, paste0("gse189926_", metric, "_seven_route_umap_comparison"), 20, 11)
}

linearity_limits <- as.numeric(quantile(linearity_values$local_linearity, c(0.01, 0.99), na.rm = TRUE, type = 8))
linearity_panels <- lapply(route_parameters$route, function(route_name) {
  values <- linearity_values[route == route_name]
  values <- values[sample.int(nrow(values))]
  continuous_panel(
    values[, .(cell_id, UMAP_1, UMAP_2)],
    values$local_linearity,
    route_name,
    linearity_limits,
    option = "B"
  )
})
linearity_figure <- wrap_plots(linearity_panels, ncol = 4, guides = "collect") +
  plot_annotation(title = "GSE189926 seven UMAP routes: local linearity on the shared stratified sample") &
  theme(legend.position = "bottom")
save_figure(linearity_figure, "gse189926_local_linearity_seven_route_umap_comparison", 20, 11)

recommended_resolution <- resolution_selection[evidence_rank == 1L]
recommended_umap <- umap_selection[evidence_rank == 1L]
status_lines <- c(
  "# GSE189926 resolution and UMAP topology sweep",
  "",
  "Status: `COMPLETED_AWAITING_MAC_REVIEW`",
  "",
  "- Exact input path, size, SHA-256, cell count, assays, layers, reductions, commands, and metadata fields were verified.",
  "- The saved object had no graph; `harmony_snn_sweep` was rebuilt from Harmony dimensions 1:30 with `k.param=15`.",
  "- Seven resolutions were calculated on the same graph and displayed on one unchanged baseline UMAP.",
  "- Seven UMAP routes were evaluated separately on the fixed Harmony reduction.",
  "- All 89,587 accepted QC-pass cells were retained; local linearity was used only as a topology diagnostic.",
  sprintf("- Evidence-ranked resolution for Mac review: `%s`; formal object unchanged.", recommended_resolution$resolution),
  sprintf("- Evidence-ranked UMAP route for Mac review: `%s`; formal object unchanged.", recommended_umap$route),
  sprintf("- Baseline q99 filament audit contains %d cells across %d resolution-0.5 clusters.", sum(filament_review$cells), nrow(filament_review)),
  "- Every PNG/PDF comparison sheet uses a white background, fixed seed 340 plotting order, alpha 0.90, and the approved bright palette.",
  "- The sweep object is local-only and does not overwrite the prior object."
)
writeLines(status_lines, file.path(output_dir, "STATUS.md"))

review_lines <- c(
  "# Mac review request: GSE189926 resolution and UMAP topology sweep",
  "",
  "## Separate decisions",
  "",
  sprintf("- Resolution evidence rank 1: `%s` (clusters=%s, mean adjacent ARI=%s, broad-marker conflict flags=%s).", recommended_resolution$resolution, recommended_resolution$clusters, signif(recommended_resolution$mean_adjacent_ari, 4), recommended_resolution$broad_marker_conflict_clusters),
  sprintf("- UMAP evidence rank 1: `%s` (q90 local linearity=%s, mean same-sample neighbor fraction=%s, centroid/dispersion ratio=%s).", recommended_umap$route, signif(recommended_umap$q90_local_linearity, 4), signif(recommended_umap$mean_same_sample_neighbor_fraction, 4), signif(recommended_umap$centroid_separation_to_dispersion_ratio, 4)),
  "- These ranks use separate, fully exported evidence tables. No choice was made from visual compactness alone.",
  "",
  "## Filament review",
  "",
  sprintf("- Baseline all-cell q99 local-linearity rule flagged %d cells; no cells were removed.", sum(filament_review$cells)),
  "- `gse189926_filament_region_review.tsv` records sample/patient dominance, QC associations, broad-label mixing, and affected cluster Top50 markers.",
  "- The local-linearity metric is a geometry diagnostic and is not used as a biological quality or exclusion score.",
  "",
  "## Requested decision",
  "",
  sprintf("1. Review whether resolution `%s` should replace the current formal resolution after inspecting marker support and transitions.", recommended_resolution$resolution),
  sprintf("2. Review whether UMAP route `%s` should replace the current display embedding after inspecting biology, sample mixing, QC overlays, and filament evidence.", recommended_umap$route),
  "3. Do not update the formal GSE189926 object until both decisions are approved."
)
writeLines(review_lines, file.path(output_dir, "MAC_REVIEW_REQUEST.md"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_gse189926_figures.txt"))
cat("GSE189926 sweep figures and review files completed\n")
