#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BPCells)
  library(data.table)
  library(harmony)
  library(Matrix)
  library(MatrixGenerics)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: 10_GSE318850_full_harmony_all_cells.R ",
    "COUNTS_H5 SOURCE_METADATA_RDS QC_DOUBLET_METADATA_RDS OUTPUT_DIR"
  )
}

counts_h5 <- normalizePath(args[[1]], mustWork = TRUE)
source_metadata_rds <- normalizePath(args[[2]], mustWork = TRUE)
qc_metadata_rds <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260814)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 14 * 1024^3)

message("Opening the full BPCells count matrix and exact metadata")
counts <- open_matrix_10x_hdf5(counts_h5)
source_metadata <- readRDS(source_metadata_rds)
qc_metadata <- as.data.frame(readRDS(qc_metadata_rds))
stopifnot(
  identical(colnames(counts), source_metadata$cell),
  identical(source_metadata$cell, qc_metadata$cell),
  all(qc_metadata$retained_singlet %in% c(TRUE, FALSE)),
  all(c("orig.ident", "Patient", "Sample") %in% colnames(qc_metadata))
)

singlet_index <- which(qc_metadata$retained_singlet)
singlet_cells <- qc_metadata$cell[singlet_index]
singlet_metadata <- qc_metadata[singlet_index, , drop = FALSE]
rownames(singlet_metadata) <- singlet_cells
singlet_counts <- counts[, singlet_index, drop = FALSE]
stopifnot(
  ncol(singlet_counts) == 649163L,
  identical(colnames(singlet_counts), rownames(singlet_metadata))
)

message("Creating the all-cell singlet Seurat object")
object <- CreateSeuratObject(
  counts = singlet_counts,
  assay = "RNA",
  project = "GSE318850_full_harmony",
  meta.data = singlet_metadata,
  min.cells = 0,
  min.features = 0
)

message("Normalizing all 649,163 singlets and selecting 2,500 variable genes")
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
  nfeatures = 2500,
  verbose = TRUE
)
hvg <- VariableFeatures(object, assay = "RNA")
stopifnot(length(hvg) == 2500L, all(hvg %in% rownames(object)))
writeLines(hvg, file.path(output_dir, "GSE318850_full_2500_HVGs.txt"))

message("Building a lazy centered and scaled all-cell matrix with BPCells")
normalized <- LayerData(object, assay = "RNA", layer = "data")
normalized_hvg <- normalized[hvg, , drop = FALSE]
feature_means <- rowMeans(normalized_hvg)
feature_sd <- sqrt(rowVars(normalized_hvg))
valid_features <- is.finite(feature_sd) & feature_sd > 0
if (!all(valid_features)) {
  hvg <- hvg[valid_features]
  normalized_hvg <- normalized_hvg[valid_features, , drop = FALSE]
  feature_means <- feature_means[valid_features]
  feature_sd <- feature_sd[valid_features]
}
scaled_hvg <- multiply_rows(
  add_rows(normalized_hvg, -feature_means),
  1 / feature_sd
)

message("Running exact all-cell truncated SVD for 40 principal components")
pca_fit <- BPCells::svds(
  t(scaled_hvg),
  k = 40L,
  opts = list(tol = 1e-4, maxitr = 2000L)
)
pca_embeddings <- sweep(pca_fit$u, 2, pca_fit$d, `*`)
rownames(pca_embeddings) <- singlet_cells
colnames(pca_embeddings) <- paste0("PC_", seq_len(ncol(pca_embeddings)))
pca_loadings <- pca_fit$v
rownames(pca_loadings) <- hvg
colnames(pca_loadings) <- colnames(pca_embeddings)
pca_stdev <- pca_fit$d / sqrt(nrow(pca_embeddings) - 1)

object[["pca.full"]] <- CreateDimReducObject(
  embeddings = pca_embeddings,
  loadings = pca_loadings,
  stdev = pca_stdev,
  key = "PCFULL_",
  assay = "RNA"
)
fwrite(
  data.table(
    component = seq_along(pca_stdev),
    stdev = pca_stdev,
    variance = pca_stdev^2,
    variance_fraction = pca_stdev^2 / sum(pca_stdev^2),
    cumulative_variance_fraction = cumsum(pca_stdev^2 / sum(pca_stdev^2))
  ),
  file.path(output_dir, "GSE318850_full_PCA_variance.tsv"),
  sep = "\t", quote = FALSE
)
rm(pca_fit, pca_embeddings, pca_loadings, normalized_hvg, scaled_hvg)
invisible(gc())

message("Saving the all-cell PCA checkpoint")
saveRDS(
  object,
  file.path(output_dir, "GSE318850_full_pca_checkpoint.rds"),
  compress = FALSE
)

message("Running Harmony on all 649,163 singlets")
object <- RunHarmony(
  object = object,
  group.by.vars = "orig.ident",
  reduction.use = "pca.full",
  dims.use = 1:30,
  reduction.save = "harmony.full",
  project.dim = FALSE,
  theta = 2,
  lambda = 1,
  max_iter = 20L,
  block.size = 0.05,
  verbose = TRUE
)

message("Computing all-cell centroid separation diagnostics before and after Harmony")
centroid_separation <- function(embedding, groups, stage, block_size = 10000L) {
  groups <- factor(groups)
  group_index <- as.integer(groups)
  group_names <- levels(groups)
  centroids <- rowsum(embedding, group = groups, reorder = FALSE) /
    as.numeric(table(groups))
  centroid_norm <- rowSums(centroids^2)
  score <- numeric(nrow(embedding))
  for (block_start in seq.int(1L, nrow(embedding), by = block_size)) {
    block_end <- min(block_start + block_size - 1L, nrow(embedding))
    block_index <- block_start:block_end
    block <- embedding[block_index, , drop = FALSE]
    distance_squared <- outer(rowSums(block^2), centroid_norm, `+`) -
      2 * tcrossprod(block, centroids)
    distance_squared[distance_squared < 0] <- 0
    own_distance <- sqrt(distance_squared[cbind(
      seq_len(nrow(distance_squared)),
      group_index[block_index]
    )])
    distance_squared[cbind(
      seq_len(nrow(distance_squared)),
      group_index[block_index]
    )] <- Inf
    nearest_other_distance <- sqrt(apply(distance_squared, 1, min))
    score[block_index] <-
      (nearest_other_distance - own_distance) /
      pmax(nearest_other_distance, own_distance, .Machine$double.eps)
  }
  data.table(
    cell = rownames(embedding),
    stage = stage,
    group = as.character(groups),
    separation_score = score
  )
}

pca_for_diagnostic <- Embeddings(object, "pca.full")[, 1:30, drop = FALSE]
harmony_for_diagnostic <- Embeddings(object, "harmony.full")[, 1:30, drop = FALSE]
batch_diagnostic <- rbindlist(list(
  centroid_separation(pca_for_diagnostic, object$orig.ident, "Before Harmony"),
  centroid_separation(harmony_for_diagnostic, object$orig.ident, "After Harmony")
))
fwrite(
  batch_diagnostic,
  file.path(output_dir, "GSE318850_full_library_centroid_separation.tsv.gz"),
  sep = "\t", quote = FALSE
)
fwrite(
  batch_diagnostic[, .(
    n_cells = .N,
    mean_score = mean(separation_score),
    median_score = median(separation_score),
    q25 = unname(quantile(separation_score, 0.25)),
    q75 = unname(quantile(separation_score, 0.75))
  ), by = .(stage)],
  file.path(output_dir, "GSE318850_full_library_centroid_separation_summary.tsv"),
  sep = "\t", quote = FALSE
)
rm(pca_for_diagnostic, harmony_for_diagnostic, batch_diagnostic)
invisible(gc())

message("Saving the all-cell PCA and Harmony checkpoint")
saveRDS(
  object,
  file.path(output_dir, "GSE318850_full_harmony_checkpoint.rds"),
  compress = FALSE
)

message("Building one all-cell Harmony neighbor index for both SNN and UMAP")
harmony_for_neighbors <- Embeddings(object, "harmony.full")[, 1:25, drop = FALSE]
full_neighbor <- FindNeighbors(
  harmony_for_neighbors,
  k.param = 30,
  nn.method = "annoy",
  n.trees = 50,
  return.neighbor = TRUE,
  verbose = TRUE
)
object[["full_nn"]] <- full_neighbor

neighbor_indices <- Indices(full_neighbor)
neighbor_j <- as.numeric(t(neighbor_indices))
neighbor_i <- ((seq_along(neighbor_j) - 1L) %/% ncol(neighbor_indices)) + 1L
full_nn_graph <- sparseMatrix(
  i = neighbor_i,
  j = neighbor_j,
  x = 1,
  dims = c(nrow(harmony_for_neighbors), nrow(harmony_for_neighbors)),
  dimnames = list(rownames(harmony_for_neighbors), rownames(harmony_for_neighbors))
)
object[["full_nn_graph"]] <- as.Graph(full_nn_graph)
full_snn <- Seurat:::ComputeSNN(
  nn_ranked = neighbor_indices,
  prune = 1 / 15
)
rownames(full_snn) <- rownames(harmony_for_neighbors)
colnames(full_snn) <- rownames(harmony_for_neighbors)
object[["full_snn"]] <- as.Graph(full_snn)
rm(
  harmony_for_neighbors, full_neighbor, neighbor_indices,
  neighbor_i, neighbor_j, full_nn_graph, full_snn
)
invisible(gc())

message("Scanning all-cell clustering resolutions 0.2 to 0.6")
resolutions <- c(0.2, 0.3, 0.4, 0.5, 0.6)
for (resolution in resolutions) {
  object <- FindClusters(
    object,
    graph.name = "full_snn",
    resolution = resolution,
    cluster.name = paste0("full_res_", resolution),
    random.seed = 20260814,
    verbose = TRUE
  )
}
object$harmony_cluster <- object$full_res_0.4
Idents(object) <- "harmony_cluster"

resolution_table <- rbindlist(lapply(resolutions, function(resolution) {
  cluster_sizes <- table(object[[paste0("full_res_", resolution), drop = TRUE]])
  data.table(
    resolution = resolution,
    n_clusters = length(cluster_sizes),
    min_cluster_size = min(cluster_sizes),
    median_cluster_size = as.numeric(median(cluster_sizes)),
    max_cluster_size = max(cluster_sizes),
    n_clusters_lt_100 = sum(cluster_sizes < 100)
  )
}))
fwrite(
  resolution_table,
  file.path(output_dir, "GSE318850_full_resolution_summary.tsv"),
  sep = "\t", quote = FALSE
)

message("Running UMAP on all 649,163 singlets using the full neighbor index")
object <- RunUMAP(
  object,
  nn.name = "full_nn",
  n.neighbors = 30,
  min.dist = 0.3,
  spread = 1,
  metric = "cosine",
  reduction.name = "umap.full.harmony",
  seed.use = 20260814,
  verbose = TRUE
)

full_metadata <- object[[]]
full_metadata$cell <- rownames(full_metadata)
fwrite(
  as.data.table(full_metadata),
  file.path(output_dir, "GSE318850_full_harmony_metadata.tsv.gz"),
  sep = "\t", quote = FALSE
)
fwrite(
  data.table(
    cell = rownames(Embeddings(object, "umap.full.harmony")),
    Embeddings(object, "umap.full.harmony")
  ),
  file.path(output_dir, "GSE318850_full_harmony_umap_coordinates.tsv.gz"),
  sep = "\t", quote = FALSE
)
fwrite(
  as.data.table(as.data.frame(table(full_metadata$harmony_cluster)))[
    , setNames(.SD, c("harmony_cluster", "n_cells"))
  ],
  file.path(output_dir, "GSE318850_full_harmony_cluster_sizes.tsv"),
  sep = "\t", quote = FALSE
)

message("Saving the all-cell UMAP stage before marker annotation")
saveRDS(
  object,
  file.path(output_dir, "GSE318850_full_harmony_stage1.rds"),
  compress = FALSE
)

message("Calculating manual-annotation markers from all 649,163 singlets")
row_major_dir <- file.path(output_dir, "RNA_normalized_row_major")
if (!dir.exists(row_major_dir)) {
  normalized_row_major <- transpose_storage_order(
    normalized,
    outdir = row_major_dir,
    tmpdir = output_dir,
    load_bytes = 16 * 1024^2,
    sort_bytes = 512 * 1024^2
  )
} else {
  normalized_row_major <- open_matrix_dir(row_major_dir)
}
marker_table <- as.data.table(BPCells::marker_features(
  normalized_row_major,
  groups = object$harmony_cluster,
  method = "wilcoxon"
))
marker_table[, p_val_adj := p.adjust(p_val_raw, method = "BH"), by = foreground]
marker_table[, avg_log_expression_difference := foreground_mean - background_mean]

group_factor <- factor(object$harmony_cluster)
membership <- sparseMatrix(
  i = seq_along(group_factor),
  j = as.integer(group_factor),
  x = 1,
  dims = c(length(group_factor), nlevels(group_factor)),
  dimnames = list(NULL, levels(group_factor))
)
detection_counts <- as.matrix((singlet_counts > 0) %*% membership)
group_sizes <- as.numeric(table(group_factor))
foreground_pct <- sweep(detection_counts, 2, group_sizes, `/`)
background_pct <- sweep(
  matrix(rowSums(detection_counts), nrow = nrow(detection_counts), ncol = ncol(detection_counts)) -
    detection_counts,
  2,
  length(group_factor) - group_sizes,
  `/`
)

cluster_index <- match(marker_table$foreground, levels(group_factor))
feature_index <- match(marker_table$feature, rownames(singlet_counts))
marker_table[, pct_in_cluster := foreground_pct[cbind(feature_index, cluster_index)]]
marker_table[, pct_outside_cluster := background_pct[cbind(feature_index, cluster_index)]]
setorder(
  marker_table,
  foreground,
  p_val_adj,
  -avg_log_expression_difference,
  -pct_in_cluster
)
fwrite(
  marker_table,
  file.path(output_dir, "GSE318850_full_harmony_cluster_markers.tsv.gz"),
  sep = "\t", quote = FALSE
)
fwrite(
  marker_table[, head(.SD, 50), by = foreground],
  file.path(output_dir, "GSE318850_full_harmony_cluster_top50_markers.tsv"),
  sep = "\t", quote = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "GSE318850_full_harmony_sessionInfo.txt")
)
message("All-cell Harmony, clustering, UMAP and marker analysis complete")
