#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(harmony)
  library(RANN)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 03b_normalize_integrate_cluster.R <gse205506_root> <output_dir> <seed>",
    call. = FALSE
  )
}

gse_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
seed <- as.integer(args[[3]])
analysis_dir <- file.path(gse_root, "r_analysis")
input_rds <- file.path(analysis_dir, "gse205506_seurat_qc_singlets.rds")
pca_rds <- file.path(analysis_dir, "gse205506_seurat_normalized_pca.rds")
harmony_rds <- file.path(analysis_dir, "gse205506_seurat_harmony.rds")
cluster_rds <- file.path(analysis_dir, "gse205506_seurat_integrated_clusters.rds")

if (!file.exists(input_rds)) stop("Missing QC singlet object: ", input_rds, call. = FALSE)

write_tsv <- function(x, filename) {
  fwrite(x, file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

object <- readRDS(input_rds)
if (ncol(object) != 238934L || length(unique(object$geo_accession)) != 40L ||
    length(unique(object$geo_subject)) != 19L) {
  stop("QC singlet object dimensions or metadata do not match the verified checkpoint", call. = FALSE)
}

DefaultAssay(object) <- "RNA"
set.seed(seed)
message("Joining 40 verified RNA count layers")
object <- JoinLayers(object, assay = "RNA")
object <- NormalizeData(
  object,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)
object <- FindVariableFeatures(
  object,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)
pca_features <- VariableFeatures(object)
object <- ScaleData(object, assay = "RNA", features = pca_features, verbose = TRUE)
object <- RunPCA(
  object,
  assay = "RNA",
  features = pca_features,
  npcs = 50,
  approx = TRUE,
  seed.use = seed,
  verbose = TRUE
)
saveRDS(object, pca_rds, compress = FALSE)
message("PCA checkpoint saved: ", pca_rds)

object <- RunHarmony(
  object,
  group.by.vars = "geo_accession",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  project.dim = FALSE,
  max.iter.harmony = 20,
  verbose = TRUE
)
saveRDS(object, harmony_rds, compress = FALSE)
message("Harmony checkpoint saved: ", harmony_rds)

eta_squared <- function(embedding, group, reduction_name) {
  group <- factor(group)
  result <- lapply(seq_len(ncol(embedding)), function(i) {
    values <- embedding[, i]
    overall <- mean(values)
    group_n <- as.numeric(table(group))
    group_means <- as.numeric(tapply(values, group, mean))
    between_ss <- sum(group_n * (group_means - overall)^2)
    total_ss <- sum((values - overall)^2)
    data.table(
      reduction = reduction_name,
      dimension = i,
      sample_eta_squared = if (total_ss > 0) between_ss / total_ss else NA_real_
    )
  })
  rbindlist(result)
}

knn_mixing <- function(embedding, sample_id, reduction_name, seed_value) {
  set.seed(seed_value)
  selected <- if (nrow(embedding) > 20000L) sample.int(nrow(embedding), 20000L) else seq_len(nrow(embedding))
  selected_embedding <- embedding[selected, 1:30, drop = FALSE]
  selected_sample <- sample_id[selected]
  nn <- nn2(selected_embedding, k = 31)$nn.idx[, -1, drop = FALSE]
  same_sample <- vapply(seq_len(nrow(nn)), function(i) {
    mean(selected_sample[nn[i, ]] == selected_sample[[i]])
  }, numeric(1))
  data.table(
    reduction = reduction_name,
    cells_evaluated = length(selected),
    neighbors_per_cell = 30L,
    mean_same_sample_neighbor_fraction = mean(same_sample),
    median_same_sample_neighbor_fraction = median(same_sample)
  )
}

pca_embedding <- Embeddings(object, "pca")[, 1:30, drop = FALSE]
harmony_embedding <- Embeddings(object, "harmony")[, 1:30, drop = FALSE]
batch_eta <- rbindlist(list(
  eta_squared(pca_embedding, object$geo_accession, "pca_uncorrected"),
  eta_squared(harmony_embedding, object$geo_accession, "harmony_by_geo_accession")
))
batch_knn <- rbindlist(list(
  knn_mixing(pca_embedding, object$geo_accession, "pca_uncorrected", seed),
  knn_mixing(harmony_embedding, object$geo_accession, "harmony_by_geo_accession", seed)
))
write_tsv(batch_eta, "gse205506_batch_structure_pc_eta_squared.tsv")
write_tsv(batch_knn, "gse205506_batch_structure_knn_audit.tsv")

set.seed(seed)
object <- RunUMAP(
  object,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  seed.use = seed,
  n_threads = 8,
  verbose = TRUE
)
object <- FindNeighbors(
  object,
  reduction = "harmony",
  dims = 1:30,
  k.param = 20,
  graph.name = c("harmony_nn", "harmony_snn"),
  verbose = TRUE
)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
cluster_rows <- vector("list", length(resolutions))
for (i in seq_along(resolutions)) {
  resolution <- resolutions[[i]]
  cluster_name <- paste0("cluster_res_", format(resolution, nsmall = 1))
  object <- FindClusters(
    object,
    graph.name = "harmony_snn",
    cluster.name = cluster_name,
    resolution = resolution,
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    random.seed = seed,
    verbose = TRUE
  )
  sizes <- table(object[[cluster_name, drop = TRUE]])
  cluster_rows[[i]] <- data.table(
    resolution = resolution,
    cluster_column = cluster_name,
    clusters = length(sizes),
    smallest_cluster_cells = min(as.integer(sizes)),
    median_cluster_cells = median(as.integer(sizes)),
    largest_cluster_cells = max(as.integer(sizes))
  )
}

object$working_cluster <- object$cluster_res_0.6
Idents(object) <- "working_cluster"
object@misc$integration_method <- list(
  uncorrected_reduction = "PCA on 3000 variable features",
  clustering_reduction = "Harmony dimensions 1:30 grouped by exact geo_accession",
  expression_values = "RNA LogNormalize data; Harmony is not used as an expression matrix",
  caution = "Sample is partially confounded with treatment, tissue state, and subject; uncorrected PCA audit is retained"
)
object@misc$working_cluster_resolution <- 0.6

write_tsv(rbindlist(cluster_rows), "gse205506_clustering_resolution_summary.tsv")
parameters <- data.table(
  parameter = c(
    "seed", "normalization_method", "normalization_scale_factor", "variable_features",
    "pca_components", "harmony_group", "harmony_dimensions", "umap_neighbors",
    "umap_min_dist", "neighbor_k", "cluster_resolutions", "working_resolution"
  ),
  value = c(
    seed, "LogNormalize", 10000, 3000, 50, "geo_accession", "1:30", 30,
    0.3, 20, paste(resolutions, collapse = ","), 0.6
  )
)
write_tsv(parameters, "gse205506_analysis_parameters.tsv")

saveRDS(object, cluster_rds, compress = FALSE)
object_summary <- data.table(
  object_path = normalizePath(cluster_rds, winslash = "/", mustWork = TRUE),
  object_size_bytes = file.info(cluster_rds)$size,
  features = nrow(object),
  cells = ncol(object),
  samples = uniqueN(object$geo_accession),
  subjects = uniqueN(object$geo_subject),
  working_resolution = 0.6,
  working_clusters = length(unique(object$working_cluster))
)
write_tsv(object_summary, "gse205506_integrated_object_summary.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_normalize_integrate_cluster.txt"))
message("Integrated clustered Seurat object completed: ", cluster_rds)
