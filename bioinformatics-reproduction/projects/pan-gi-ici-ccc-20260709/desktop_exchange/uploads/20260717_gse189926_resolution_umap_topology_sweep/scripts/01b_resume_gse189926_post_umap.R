#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
output_dir <- file.path(project_dir, "desktop_exchange/uploads/20260717_gse189926_resolution_umap_topology_sweep")
object_dir <- file.path(project_dir, "04_objects/20260717_gse189926_resolution_umap_topology_sweep")
checkpoint <- file.path(object_dir, "gse189926_sweep_post_umap_checkpoint.rds")
output_object <- file.path(object_dir, "gse189926_resolution_umap_topology_sweep.rds")

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

if (!file.exists(checkpoint)) stop("Post-UMAP checkpoint is missing", call. = FALSE)
object <- readRDS(checkpoint)
metadata <- as.data.table(object[[]], keep.rownames = "cell_id")
fixed_cluster <- "cluster_harmony_r0_5"
route_parameters <- fread(file.path(output_dir, "gse189926_umap_route_parameters.tsv"))
top50 <- fread(file.path(output_dir, "gse189926_resolution_top50_markers.tsv"))
filament_cells <- fread(file.path(output_dir, "gse189926_filament_region_cell_audit.tsv.gz"))
linearity_summary <- fread(file.path(output_dir, "gse189926_umap_local_linearity_summary.tsv"))
mixing_metrics <- fread(file.path(output_dir, "gse189926_umap_route_mixing_metrics.tsv"))
geometry_summary <- fread(file.path(output_dir, "gse189926_umap_fixed_cluster_geometry.tsv"))
resolution_summary <- fread(file.path(output_dir, "gse189926_resolution_sweep_summary.tsv"))

overall_q95 <- metadata[, .(
  pct_mt_q95 = as.numeric(quantile(pct_mt, 0.95, type = 8)),
  pct_ribo_q95 = as.numeric(quantile(pct_ribo, 0.95, type = 8)),
  pct_hb_q95 = as.numeric(quantile(pct_hb, 0.95, type = 8))
)]
filament_review <- filament_cells[, {
  sample_table <- sort(table(sample_accession), decreasing = TRUE)
  patient_table <- sort(table(patient_id), decreasing = TRUE)
  broad_table <- sort(table(broad_label), decreasing = TRUE)
  cluster_value <- get(fixed_cluster)[[1L]]
  evidence <- c()
  if (sample_table[[1L]] / .N >= 0.5) evidence <- c(evidence, "one_sample_dominance")
  if (patient_table[[1L]] / .N >= 0.5) evidence <- c(evidence, "one_patient_dominance")
  if (median(pct_mt) >= overall_q95$pct_mt_q95) evidence <- c(evidence, "high_mitochondrial_metric")
  if (median(pct_ribo) >= overall_q95$pct_ribo_q95) evidence <- c(evidence, "high_ribosomal_metric")
  if (median(pct_hb) >= overall_q95$pct_hb_q95) evidence <- c(evidence, "high_hemoglobin_metric")
  if (broad_table[[1L]] / .N < 0.6) evidence <- c(evidence, "mixed_broad_labels")
  if (!length(evidence)) evidence <- "geometry_only_review"
  data.table(
    cells = .N,
    maximum_sample = names(sample_table)[[1L]],
    maximum_sample_fraction = as.numeric(sample_table[[1L]] / .N),
    maximum_patient = names(patient_table)[[1L]],
    maximum_patient_fraction = as.numeric(patient_table[[1L]] / .N),
    dominant_broad_label = names(broad_table)[[1L]],
    dominant_broad_fraction = as.numeric(broad_table[[1L]] / .N),
    median_nFeature_RNA = as.numeric(median(nFeature_RNA)),
    median_nCount_RNA = as.numeric(median(nCount_RNA)),
    median_pct_mt = as.numeric(median(pct_mt)),
    median_pct_ribo = as.numeric(median(pct_ribo)),
    median_pct_hb = as.numeric(median(pct_hb)),
    median_local_linearity = as.numeric(median(local_linearity)),
    evidence_labels = paste(evidence, collapse = ";"),
    action = "retain cells and request review; no automatic rename or deletion",
    top50_markers = paste(top50[resolution == "0.5" & cluster == cluster_value][order(rank), gene], collapse = ";")
  )
}, by = .(cluster_harmony_r0_5)]
write_tsv(filament_review, "gse189926_filament_region_review.tsv")

resolution_selection <- copy(resolution_summary)
resolution_selection[, `:=`(
  rank_marker_support = frank(broad_marker_conflict_clusters, ties.method = "min"),
  rank_small_clusters = frank(clusters_below_100, ties.method = "min"),
  rank_stability = frank(-mean_adjacent_ari, ties.method = "min"),
  rank_sample_representation = frank(median_max_sample_fraction, ties.method = "min"),
  rank_patient_representation = frank(median_max_patient_fraction, ties.method = "min")
)]
resolution_selection[, evidence_rank_sum := rank_marker_support + rank_small_clusters + rank_stability + rank_sample_representation + rank_patient_representation]
resolution_selection[, resolution_numeric := as.numeric(resolution)]
setorder(resolution_selection, evidence_rank_sum, resolution_numeric)
resolution_selection[, evidence_rank := seq_len(.N)]
resolution_selection[, resolution_numeric := NULL]
resolution_selection[, recommendation := fifelse(evidence_rank == 1L, "recommended for Mac review; formal object unchanged", "retain as evaluated alternative")]
write_tsv(resolution_selection, "gse189926_resolution_selection_evidence.tsv")

umap_selection <- merge(linearity_summary, mixing_metrics, by = c("route", "sampled_cells"))
umap_selection <- merge(umap_selection, geometry_summary, by = "route")
baseline_patient <- umap_selection[route == "baseline", mean_same_patient_neighbor_fraction]
umap_selection[, patient_structure_ratio_to_baseline := mean_same_patient_neighbor_fraction / baseline_patient]
umap_selection[, `:=`(
  rank_linearity = frank(q90_local_linearity, ties.method = "min"),
  rank_sample_mixing = frank(mean_same_sample_neighbor_fraction, ties.method = "min"),
  rank_patient_structure_retention = frank(abs(log(patient_structure_ratio_to_baseline)), ties.method = "min"),
  rank_cluster_geometry = frank(-centroid_separation_to_dispersion_ratio, ties.method = "min")
)]
umap_selection[, evidence_rank_sum := rank_linearity + rank_sample_mixing + rank_patient_structure_retention + rank_cluster_geometry]
setorder(umap_selection, evidence_rank_sum, route)
umap_selection[, evidence_rank := seq_len(.N)]
umap_selection[, recommendation := fifelse(evidence_rank == 1L, "recommended for Mac review; formal object unchanged", "retain as evaluated alternative")]
write_tsv(umap_selection, "gse189926_umap_selection_evidence.tsv")

object@misc$gse189926_sweep <- list(
  request_date = "2026-07-17",
  input_sha256 = "bf112efde4b1816348c8939ab2127f00ae1e959f4132adda17a74d8673052abb",
  graph = "harmony_snn_sweep rebuilt from harmony dimensions 1:30 with k.param=15",
  resolutions = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0),
  routes = route_parameters,
  local_linearity_definition = "lambda1 / max(lambda2, .Machine$double.eps) from 30 neighbors in 2D UMAP coordinates",
  no_cells_removed = TRUE
)
saveRDS(object, output_object, compress = FALSE)
object_summary <- data.table(
  object_path = normalizePath(output_object, winslash = "/", mustWork = TRUE),
  size_bytes = as.numeric(file.info(output_object)$size),
  sha256 = digest::digest(file = output_object, algo = "sha256", serialize = FALSE),
  cells = ncol(object),
  features = nrow(object),
  resolution_fields = paste(c("sweep_res_0_2", "sweep_res_0_3", "sweep_res_0_4", "sweep_res_0_5", "sweep_res_0_6", "sweep_res_0_8", "sweep_res_1_0"), collapse = ";"),
  umap_reductions = paste(route_parameters$reduction_name, collapse = ";"),
  delivery = "local_only"
)
write_tsv(object_summary, "gse189926_sweep_object_summary.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_gse189926_sweep.txt"))
cat("GSE189926 post-UMAP resume completed\n")
