#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(RANN)
  library(presto)
  library(mclust)
})

options(stringsAsFactors = FALSE, future.globals.maxSize = 16 * 1024^3)
set.seed(340)

project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
input_path <- file.path(project_dir, "04_objects/20260710_gse189926_r_only_rerun/gse189926_r_qc_pass_annotated.rds")
output_dir <- file.path(project_dir, "desktop_exchange/uploads/20260717_gse189926_resolution_umap_topology_sweep")
object_dir <- file.path(project_dir, "04_objects/20260717_gse189926_resolution_umap_topology_sweep")
output_object <- file.path(object_dir, "gse189926_resolution_umap_topology_sweep.rds")
pre_umap_checkpoint <- file.path(object_dir, "gse189926_sweep_pre_umap_checkpoint.rds")
post_umap_checkpoint <- file.path(object_dir, "gse189926_sweep_post_umap_checkpoint.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

log_step <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
  flush.console()
}

resolution_values <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
resolution_names <- c("0.2", "0.3", "0.4", "0.5", "0.6", "0.8", "1.0")
resolution_fields <- setNames(paste0("sweep_res_", gsub("\\.", "_", resolution_names)), resolution_names)

route_parameters <- data.table(
  route = c("baseline", "neighbor30", "neighbor50", "compact30", "diffuse30", "diffuse50", "dims20_neighbor30"),
  reduction_name = c(
    "umap.harmony", "umap.sweep.neighbor30", "umap.sweep.neighbor50",
    "umap.sweep.compact30", "umap.sweep.diffuse30", "umap.sweep.diffuse50",
    "umap.sweep.dims20_neighbor30"
  ),
  harmony_dims = c(30L, 30L, 30L, 30L, 30L, 30L, 20L),
  n_neighbors = c(15L, 30L, 50L, 30L, 30L, 50L, 30L),
  min_dist = c(0.3, 0.3, 0.3, 0.1, 0.5, 0.5, 0.3),
  seed = 340L,
  method = "uwot",
  metric = "cosine",
  thread_count = 1L,
  thread_count_source = "Seurat 5.4.0 RunUMAP calls uwot with future::nbrOfWorkers(); default plan reported one worker; direct n_threads argument is not accepted",
  execution = c("verified existing baseline", rep("calculated in this sweep", 6L)),
  uwot_version = as.character(packageVersion("uwot")),
  Seurat_version = as.character(packageVersion("Seurat"))
)

expected_size <- 3739942614
expected_sha <- "bf112efde4b1816348c8939ab2127f00ae1e959f4132adda17a74d8673052abb"
if (!file.exists(input_path)) stop("Exact requested input object is missing", call. = FALSE)
actual_size <- as.numeric(file.info(input_path)$size)
if (actual_size != expected_size) stop("Input object size mismatch: ", actual_size, call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("R package digest is required", call. = FALSE)
actual_sha <- digest::digest(file = input_path, algo = "sha256", serialize = FALSE)
if (!identical(tolower(actual_sha), expected_sha)) stop("Input object SHA-256 mismatch", call. = FALSE)

log_step("Reading verified GSE189926 object")
object <- readRDS(input_path)
required_metadata <- c(
  "sample_accession", "patient_id", "timepoint", "characteristics_ch1::outcome",
  "nFeature_RNA", "nCount_RNA", "pct_mt", "pct_ribo", "pct_hb",
  "cluster_harmony_r0_5", "broad_label", "refined_label"
)
missing_metadata <- setdiff(required_metadata, colnames(object[[]]))
if (length(missing_metadata)) stop("Missing exact metadata fields: ", paste(missing_metadata, collapse = "; "), call. = FALSE)
if (ncol(object) != 89587L) stop("Expected 89587 cells, observed ", ncol(object), call. = FALSE)
if (!"RNA" %in% Assays(object)) stop("RNA assay is missing", call. = FALSE)
if (!all(c("counts", "data") %in% Layers(object[["RNA"]]))) stop("Required RNA layers are missing", call. = FALSE)
if (!all(c("harmony", "umap.harmony") %in% Reductions(object))) stop("Required reductions are missing", call. = FALSE)
if (ncol(Embeddings(object, "harmony")) < 30L) stop("Harmony has fewer than 30 dimensions", call. = FALSE)

baseline_umap_command <- object@commands[["RunUMAP.RNA.harmony"]]@params
baseline_neighbor_command <- object@commands[["FindNeighbors.RNA.harmony"]]@params
baseline_audit <- data.table(
  metric = c(
    "input_path", "input_size_bytes", "input_sha256", "cells", "features", "assays", "RNA_layers",
    "reductions", "graphs_observed_before_sweep", "metadata_fields", "harmony_dimensions",
    "harmony_grouping_source", "baseline_umap_n_neighbors", "baseline_umap_min_dist",
    "baseline_umap_seed", "baseline_umap_metric", "baseline_neighbor_k", "baseline_neighbor_dims",
    "graph_action"
  ),
  value = c(
    normalizePath(input_path, winslash = "/", mustWork = TRUE), as.character(actual_size), actual_sha,
    as.character(ncol(object)), as.character(nrow(object)), paste(Assays(object), collapse = ";"),
    paste(Layers(object[["RNA"]]), collapse = ";"), paste(Reductions(object), collapse = ";"),
    paste(Graphs(object), collapse = ";"), paste(colnames(object[[]]), collapse = ";"),
    as.character(ncol(Embeddings(object, "harmony"))),
    "sample_accession recorded by the accepted R-only build",
    as.character(baseline_umap_command$n.neighbors), as.character(baseline_umap_command$min.dist),
    as.character(baseline_umap_command$seed.use), as.character(baseline_umap_command$metric),
    as.character(baseline_neighbor_command$k.param), paste(range(baseline_neighbor_command$dims), collapse = ":"),
    "harmony_snn was absent from the saved object and was rebuilt from harmony dimensions 1:30 with k.param=15"
  ),
  pass = c(rep(TRUE, 18L), TRUE)
)
write_tsv(baseline_audit, "gse189926_sweep_input_object_audit.tsv")

log_step("Rebuilding Harmony graph from exact requested parameters")
object <- FindNeighbors(
  object,
  reduction = "harmony",
  dims = 1:30,
  k.param = 15,
  graph.name = c("harmony_nn_sweep", "harmony_snn_sweep"),
  verbose = FALSE
)

log_step("Running fixed-graph resolution sweep")
for (i in seq_along(resolution_values)) {
  object <- FindClusters(
    object,
    graph.name = "harmony_snn_sweep",
    resolution = resolution_values[[i]],
    algorithm = 1,
    random.seed = 340,
    cluster.name = resolution_fields[[resolution_names[[i]]]],
    verbose = FALSE
  )
}

metadata <- as.data.table(object[[]], keep.rownames = "cell_id")
for (field in resolution_fields) metadata[, (field) := as.character(get(field))]

resolution_cluster_counts <- rbindlist(lapply(resolution_names, function(resolution_name) {
  field <- resolution_fields[[resolution_name]]
  out <- metadata[, .(cells = .N), by = .(cluster = get(field))]
  out[, resolution := resolution_name]
  out[]
}))
resolution_cluster_counts[, cluster_numeric := as.integer(cluster)]
resolution_cluster_counts[, resolution_numeric := as.numeric(resolution)]
setorder(resolution_cluster_counts, resolution_numeric, cluster_numeric)
resolution_cluster_counts[, c("cluster_numeric", "resolution_numeric") := NULL]
write_tsv(resolution_cluster_counts, "gse189926_resolution_cluster_counts.tsv")

resolution_dominance <- rbindlist(lapply(resolution_names, function(resolution_name) {
  field <- resolution_fields[[resolution_name]]
  base <- metadata[, .(
    cells = .N,
    samples = uniqueN(sample_accession),
    patients = uniqueN(patient_id)
  ), by = .(cluster = get(field))]
  sample_counts <- metadata[, .N, by = .(cluster = get(field), sample_accession)]
  patient_counts <- metadata[, .N, by = .(cluster = get(field), patient_id)]
  sample_max <- sample_counts[, .SD[which.max(N)], by = cluster]
  patient_max <- patient_counts[, .SD[which.max(N)], by = cluster]
  setnames(sample_max, c("sample_accession", "N"), c("max_sample", "max_sample_cells"))
  setnames(patient_max, c("patient_id", "N"), c("max_patient", "max_patient_cells"))
  out <- merge(base, sample_max, by = "cluster", sort = FALSE)
  out <- merge(out, patient_max, by = "cluster", sort = FALSE)
  out[, `:=`(
    resolution = resolution_name,
    max_sample_fraction = max_sample_cells / cells,
    max_patient_fraction = max_patient_cells / cells
  )]
  out[]
}))
resolution_dominance[, cluster_numeric := as.integer(cluster)]
resolution_dominance[, resolution_numeric := as.numeric(resolution)]
setorder(resolution_dominance, resolution_numeric, cluster_numeric)
resolution_dominance[, c("cluster_numeric", "resolution_numeric") := NULL]
write_tsv(resolution_dominance, "gse189926_resolution_sample_patient_dominance.tsv")

resolution_transitions <- rbindlist(lapply(seq_len(length(resolution_names) - 1L), function(i) {
  from_name <- resolution_names[[i]]
  to_name <- resolution_names[[i + 1L]]
  out <- metadata[, .N, by = .(
    from_cluster = get(resolution_fields[[from_name]]),
    to_cluster = get(resolution_fields[[to_name]])
  )]
  out[, `:=`(from_resolution = from_name, to_resolution = to_name)]
  out[]
}))
setnames(resolution_transitions, "N", "cells")
write_tsv(resolution_transitions, "gse189926_resolution_transition_tables.tsv")

adjacent_ari <- rbindlist(lapply(seq_len(length(resolution_names) - 1L), function(i) {
  from_name <- resolution_names[[i]]
  to_name <- resolution_names[[i + 1L]]
  data.table(
    from_resolution = from_name,
    to_resolution = to_name,
    adjusted_rand_index = adjustedRandIndex(
      metadata[[resolution_fields[[from_name]]]],
      metadata[[resolution_fields[[to_name]]]]
    ),
    package = paste0("mclust ", packageVersion("mclust"))
  )
}))
write_tsv(adjacent_ari, "gse189926_resolution_adjacent_ari.tsv")

log_step("Calculating full-cell RNA markers for all resolutions")
DefaultAssay(object) <- "RNA"
marker_results <- lapply(resolution_names, function(resolution_name) {
  field <- resolution_fields[[resolution_name]]
  log_step(paste("Marker analysis resolution", resolution_name))
  marker_all <- as.data.table(wilcoxauc(object, group_by = field, assay = "data", seurat_assay = "RNA"))
  setnames(marker_all, c("feature", "group"), c("gene", "cluster"))
  marker_all[, `:=`(
    resolution = resolution_name,
    cluster = as.character(cluster),
    positive_direction = auc > 0.5 & logFC > 0,
    passes_primary_filter = padj <= 0.05 & auc > 0.5 & logFC > 0 & pct_in >= 5
  )]
  marker_all[, cluster_numeric := as.integer(cluster)]
  setorder(marker_all, cluster_numeric, -passes_primary_filter, -positive_direction, -auc, -logFC, -pct_in, gene)
  marker_all[, rank := seq_len(.N), by = cluster]
  marker_all[, cluster_numeric := NULL]
  list(top10 = marker_all[rank <= 10L], top50 = marker_all[rank <= 50L])
})
top10 <- rbindlist(lapply(marker_results, `[[`, "top10"), use.names = TRUE, fill = TRUE)
top50 <- rbindlist(lapply(marker_results, `[[`, "top50"), use.names = TRUE, fill = TRUE)
write_tsv(top10, "gse189926_resolution_top10_markers.tsv")
write_tsv(top50, "gse189926_resolution_top50_markers.tsv")

marker_sets <- list(
  Granulocyte = c("FCGR3B", "CSF3R", "S100A8", "S100A9"),
  T_cell = c("CD3D", "CD3E", "TRAC", "IL7R", "CCR7"),
  CD8_T_NK = c("CD3D", "CD8A", "NKG7", "GNLY", "CCL5"),
  Plasma = c("JCHAIN", "MZB1", "XBP1", "SDC1"),
  Inflammatory_myeloid = c("LYZ", "S100A8", "S100A9", "FCGR3A"),
  Macrophage = c("LST1", "C1QA", "C1QB", "C1QC", "APOE"),
  Monocyte = c("LYZ", "LST1", "S100A8", "S100A9", "FCN1"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  Endothelial = c("PECAM1", "VWF", "PLVAP", "EMCN"),
  Stromal_perivascular = c("COL1A1", "COL1A2", "COL3A1", "RGS5", "ACTA2"),
  Dendritic = c("CD1C", "FCER1A", "LAMP3", "FSCN1")
)

broad_marker_evidence <- rbindlist(lapply(resolution_names, function(resolution_name) {
  clusters <- sort(unique(top50[resolution == resolution_name, cluster]))
  rbindlist(lapply(clusters, function(cluster_value) {
    genes <- top50[resolution == resolution_name & cluster == cluster_value, gene]
    overlaps <- vapply(marker_sets, function(markers) length(intersect(genes, markers)), integer(1))
    ordering <- order(-overlaps, names(overlaps))
    top_name <- names(overlaps)[ordering[[1L]]]
    second_name <- names(overlaps)[ordering[[2L]]]
    top_count <- overlaps[[ordering[[1L]]]]
    second_count <- overlaps[[ordering[[2L]]]]
    data.table(
      resolution = resolution_name,
      cluster = cluster_value,
      strongest_marker_program = top_name,
      strongest_overlap_count = top_count,
      strongest_overlap_genes = paste(intersect(genes, marker_sets[[top_name]]), collapse = ";"),
      second_marker_program = second_name,
      second_overlap_count = second_count,
      second_overlap_genes = paste(intersect(genes, marker_sets[[second_name]]), collapse = ";"),
      broad_marker_conflict_review = top_count == 0L || (second_count >= 2L && second_count >= 0.75 * top_count),
      conflict_rule = "review if no canonical overlap or second program has >=2 overlaps and >=75% of strongest overlap; review flag only, no cell deletion"
    )
  }))
}))
write_tsv(broad_marker_evidence, "gse189926_resolution_broad_marker_evidence.tsv")

resolution_summary <- rbindlist(lapply(resolution_names, function(resolution_name) {
  counts <- resolution_cluster_counts[resolution == resolution_name]
  dominance <- resolution_dominance[resolution == resolution_name]
  evidence <- broad_marker_evidence[resolution == resolution_name]
  ari_values <- c(
    adjacent_ari[to_resolution == resolution_name, adjusted_rand_index],
    adjacent_ari[from_resolution == resolution_name, adjusted_rand_index]
  )
  data.table(
    resolution = resolution_name,
    clusters = nrow(counts),
    clusters_below_50 = sum(counts$cells < 50L),
    clusters_below_100 = sum(counts$cells < 100L),
    clusters_below_200 = sum(counts$cells < 200L),
    smallest_cluster_cells = min(counts$cells),
    median_cluster_cells = median(counts$cells),
    median_max_sample_fraction = median(dominance$max_sample_fraction),
    maximum_sample_fraction = max(dominance$max_sample_fraction),
    median_max_patient_fraction = median(dominance$max_patient_fraction),
    maximum_patient_fraction = max(dominance$max_patient_fraction),
    broad_marker_conflict_clusters = sum(evidence$broad_marker_conflict_review),
    mean_adjacent_ari = mean(ari_values)
  )
}))
write_tsv(resolution_summary, "gse189926_resolution_sweep_summary.tsv")
saveRDS(object, pre_umap_checkpoint, compress = FALSE)

log_step("Calculating six new fixed-Harmony UMAP routes")
for (i in 2:nrow(route_parameters)) {
  row <- route_parameters[i]
  object <- RunUMAP(
    object,
    reduction = "harmony",
    dims = seq_len(row$harmony_dims),
    n.neighbors = row$n_neighbors,
    min.dist = row$min_dist,
    metric = row$metric,
    seed.use = row$seed,
    reduction.name = row$reduction_name,
    reduction.key = paste0("SWEEP", i, "_"),
    verbose = FALSE
  )
}
write_tsv(route_parameters, "gse189926_umap_route_parameters.tsv")
saveRDS(object, post_umap_checkpoint, compress = FALSE)

stratified_cell_sample <- function(meta, group_field, maximum_cells = 20000L, seed = 340L) {
  if (nrow(meta) <= maximum_cells) return(meta$cell_id)
  group_sizes <- meta[, .N, by = group_field]
  group_sizes[, allocation := floor(maximum_cells * N / sum(N))]
  group_sizes[allocation < 1L, allocation := 1L]
  remaining <- maximum_cells - sum(group_sizes$allocation)
  if (remaining > 0L) {
    group_sizes[, fractional := maximum_cells * N / sum(N) - floor(maximum_cells * N / sum(N))]
    setorderv(group_sizes, c("fractional", group_field), c(-1L, 1L))
    group_sizes[seq_len(remaining), allocation := allocation + 1L]
  }
  if (remaining < 0L) {
    setorderv(group_sizes, c("allocation", group_field), c(1L, 1L))
    reducible <- which(group_sizes$allocation > 1L)
    group_sizes[reducible[seq_len(-remaining)], allocation := allocation - 1L]
  }
  set.seed(seed)
  unlist(lapply(seq_len(nrow(group_sizes)), function(i) {
    values <- meta[get(group_field) == group_sizes[[group_field]][[i]], cell_id]
    sample(values, group_sizes$allocation[[i]], replace = FALSE)
  }), use.names = FALSE)
}

local_linearity <- function(coords, k = 30L) {
  neighbors <- RANN::nn2(coords, k = k + 1L)$nn.idx[, -1L, drop = FALSE]
  x <- matrix(coords[neighbors, 1L], nrow = nrow(neighbors))
  y <- matrix(coords[neighbors, 2L], nrow = nrow(neighbors))
  mean_x <- rowMeans(x)
  mean_y <- rowMeans(y)
  var_x <- rowMeans((x - mean_x)^2)
  var_y <- rowMeans((y - mean_y)^2)
  cov_xy <- rowMeans((x - mean_x) * (y - mean_y))
  root <- sqrt((var_x - var_y)^2 + 4 * cov_xy^2)
  lambda1 <- (var_x + var_y + root) / 2
  lambda2 <- (var_x + var_y - root) / 2
  lambda1 / pmax(lambda2, .Machine$double.eps)
}

fixed_cluster <- "cluster_harmony_r0_5"
sampled_cells <- stratified_cell_sample(metadata, fixed_cluster, 20000L, 340L)
sampled_meta <- metadata[match(sampled_cells, cell_id)]

route_cell_values <- list()
mixing_rows <- list()
linearity_summary_rows <- list()
linearity_cluster_rows <- list()
geometry_rows <- list()
for (i in seq_len(nrow(route_parameters))) {
  route <- route_parameters$route[[i]]
  reduction_name <- route_parameters$reduction_name[[i]]
  coords_all <- Embeddings(object, reduction = reduction_name)
  coords <- coords_all[sampled_cells, , drop = FALSE]
  linearity <- local_linearity(coords, 30L)
  neighbors <- RANN::nn2(coords, k = 31L)$nn.idx[, -1L, drop = FALSE]
  same_sample <- rowMeans(matrix(sampled_meta$sample_accession[neighbors], nrow = nrow(neighbors)) == sampled_meta$sample_accession)
  same_patient <- rowMeans(matrix(sampled_meta$patient_id[neighbors], nrow = nrow(neighbors)) == sampled_meta$patient_id)
  route_values <- data.table(
    cell_id = sampled_cells,
    route = route,
    UMAP_1 = coords[, 1L],
    UMAP_2 = coords[, 2L],
    local_linearity = linearity,
    cluster_harmony_r0_5 = sampled_meta[[fixed_cluster]],
    sample_accession = sampled_meta$sample_accession,
    patient_id = sampled_meta$patient_id
  )
  route_cell_values[[route]] <- route_values
  mixing_rows[[route]] <- data.table(
    route = route,
    sampled_cells = length(sampled_cells),
    neighbors = 30L,
    mean_same_sample_neighbor_fraction = mean(same_sample),
    median_same_sample_neighbor_fraction = median(same_sample),
    mean_same_patient_neighbor_fraction = mean(same_patient),
    median_same_patient_neighbor_fraction = median(same_patient),
    sample_pure_connected_regions = NA_integer_,
    sample_pure_connected_regions_note = "not computed because no existing project diagnostic method was found; no substitute implementation used"
  )
  linearity_summary_rows[[route]] <- data.table(
    route = route,
    sampled_cells = length(linearity),
    median_local_linearity = median(linearity),
    q90_local_linearity = as.numeric(quantile(linearity, 0.90, type = 8)),
    q99_local_linearity = as.numeric(quantile(linearity, 0.99, type = 8)),
    fraction_above_10 = mean(linearity > 10),
    fraction_above_20 = mean(linearity > 20),
    fraction_above_50 = mean(linearity > 50)
  )
  linearity_cluster_rows[[route]] <- route_values[, .(
    cells = .N,
    median_local_linearity = median(local_linearity),
    q90_local_linearity = as.numeric(quantile(local_linearity, 0.90, type = 8)),
    q99_local_linearity = as.numeric(quantile(local_linearity, 0.99, type = 8)),
    fraction_above_10 = mean(local_linearity > 10),
    fraction_above_20 = mean(local_linearity > 20),
    fraction_above_50 = mean(local_linearity > 50)
  ), by = .(route, cluster_harmony_r0_5)]

  cluster_geometry <- route_values[, .(
    centroid_x = mean(UMAP_1),
    centroid_y = mean(UMAP_2),
    within_cluster_dispersion = mean(sqrt((UMAP_1 - mean(UMAP_1))^2 + (UMAP_2 - mean(UMAP_2))^2))
  ), by = cluster_harmony_r0_5]
  centroid_matrix <- as.matrix(cluster_geometry[, .(centroid_x, centroid_y)])
  geometry_rows[[route]] <- data.table(
    route = route,
    fixed_cluster_field = fixed_cluster,
    cluster_centroid_mean_pairwise_separation = mean(as.numeric(dist(centroid_matrix))),
    median_within_cluster_dispersion = median(cluster_geometry$within_cluster_dispersion),
    centroid_separation_to_dispersion_ratio = mean(as.numeric(dist(centroid_matrix))) / median(cluster_geometry$within_cluster_dispersion)
  )
}

route_cell_values_dt <- rbindlist(route_cell_values)
fwrite(route_cell_values_dt, file.path(output_dir, "gse189926_umap_local_linearity_cell_values.tsv.gz"), sep = "\t", quote = FALSE, na = "", compress = "gzip")
mixing_metrics <- rbindlist(mixing_rows)
linearity_summary <- rbindlist(linearity_summary_rows)
linearity_by_cluster <- rbindlist(linearity_cluster_rows)
geometry_summary <- rbindlist(geometry_rows)
write_tsv(mixing_metrics, "gse189926_umap_route_mixing_metrics.tsv")
write_tsv(linearity_summary, "gse189926_umap_local_linearity_summary.tsv")
write_tsv(linearity_by_cluster, "gse189926_umap_local_linearity_by_cluster.tsv")
write_tsv(geometry_summary, "gse189926_umap_fixed_cluster_geometry.tsv")

log_step("Auditing baseline filament-like cells at the all-cell level")
baseline_coords <- Embeddings(object, reduction = "umap.harmony")
baseline_linearity_all <- local_linearity(baseline_coords, 30L)
filament_threshold <- as.numeric(quantile(baseline_linearity_all, 0.99, type = 8))
filament_index <- which(baseline_linearity_all >= filament_threshold)
doublet_fields <- grep("doublet|scDblFinder", colnames(metadata), ignore.case = TRUE, value = TRUE)
filament_columns <- unique(c(
  "cell_id", fixed_cluster, "broad_label", "refined_label", "sample_accession", "patient_id",
  "characteristics_ch1::outcome", "timepoint", "nFeature_RNA", "nCount_RNA", "pct_mt", "pct_ribo", "pct_hb",
  doublet_fields
))
filament_cells <- metadata[filament_index, ..filament_columns]
filament_cells[, `:=`(
  baseline_umap_1 = baseline_coords[filament_index, 1L],
  baseline_umap_2 = baseline_coords[filament_index, 2L],
  local_linearity = baseline_linearity_all[filament_index],
  selection_rule = paste0("baseline local_linearity >= all-cell q99 = ", signif(filament_threshold, 8))
)]
fwrite(filament_cells, file.path(output_dir, "gse189926_filament_region_cell_audit.tsv.gz"), sep = "\t", quote = FALSE, na = "", compress = "gzip")

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
    maximum_sample_fraction = sample_table[[1L]] / .N,
    maximum_patient = names(patient_table)[[1L]],
    maximum_patient_fraction = patient_table[[1L]] / .N,
    dominant_broad_label = names(broad_table)[[1L]],
    dominant_broad_fraction = broad_table[[1L]] / .N,
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
  input_sha256 = actual_sha,
  graph = "harmony_snn_sweep rebuilt from harmony dimensions 1:30 with k.param=15",
  resolutions = resolution_values,
  routes = route_parameters,
  stratified_sample_cells = sampled_cells,
  local_linearity_definition = "lambda1 / max(lambda2, .Machine$double.eps) from 30 neighbors in 2D UMAP coordinates",
  no_cells_removed = TRUE
)

log_step("Saving sweep object without overwriting the accepted input")
saveRDS(object, output_object, compress = FALSE)
object_summary <- data.table(
  object_path = normalizePath(output_object, winslash = "/", mustWork = TRUE),
  size_bytes = as.numeric(file.info(output_object)$size),
  sha256 = digest::digest(file = output_object, algo = "sha256", serialize = FALSE),
  cells = ncol(object),
  features = nrow(object),
  resolution_fields = paste(unname(resolution_fields), collapse = ";"),
  umap_reductions = paste(route_parameters$reduction_name, collapse = ";"),
  delivery = "local_only"
)
write_tsv(object_summary, "gse189926_sweep_object_summary.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_gse189926_sweep.txt"))
log_step("GSE189926 sweep analysis complete")
