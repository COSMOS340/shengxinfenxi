#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(presto)
  library(RANN)
})

options(stringsAsFactors = FALSE, future.globals.maxSize = 32 * 1024^3)
future::plan("sequential")
set.seed(340)

gse_root <- "F:/pan-gi-ici-ccc-20260709/GSE205506"
project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
output_dir <- file.path(project_dir, "desktop_exchange/uploads/20260717_gse205506_paper_qc_rebuild")
analysis_dir <- file.path(gse_root, "r_analysis_paper_qc_20260717")
input_path <- file.path(analysis_dir, "gse205506_post_compartment_mt_qc_raw.rds")
final_path <- file.path(analysis_dir, "gse205506_formal_post_mt_rpca_atlas.rds")
anchor_checkpoint_path <- file.path(analysis_dir, "gse205506_final_rpca_anchors_light_checkpoint.rds")
light_checkpoint_path <- file.path(analysis_dir, "gse205506_final_rpca_light_checkpoint.rds")
integration_output_feature_count <- 1000L

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

log_step <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
  flush.console()
}

eta_squared <- function(embedding, group, reduction_name) {
  group <- factor(group)
  rbindlist(lapply(seq_len(ncol(embedding)), function(index) {
    values <- embedding[, index]
    total_ss <- sum((values - mean(values))^2)
    group_means <- tapply(values, group, mean)
    group_sizes <- table(group)
    between_ss <- sum(as.numeric(group_sizes) * (group_means - mean(values))^2)
    data.table(
      reduction = reduction_name,
      dimension = index,
      sample_eta_squared = if (total_ss > 0) between_ss / total_ss else NA_real_
    )
  }))
}

knn_audit <- function(embedding, group, reduction_name, maximum_cells = 20000L, neighbors = 30L) {
  set.seed(340)
  selected <- if (nrow(embedding) > maximum_cells) sample(seq_len(nrow(embedding)), maximum_cells) else seq_len(nrow(embedding))
  selected_embedding <- embedding[selected, , drop = FALSE]
  selected_group <- as.character(group[selected])
  nearest <- RANN::nn2(selected_embedding, k = neighbors + 1L)$nn.idx[, -1L, drop = FALSE]
  same_sample <- vapply(seq_len(nrow(nearest)), function(index) {
    mean(selected_group[nearest[index, ]] == selected_group[[index]])
  }, numeric(1))
  data.table(
    reduction = reduction_name,
    cells_evaluated = length(selected),
    neighbors_per_cell = neighbors,
    mean_same_sample_neighbor_fraction = mean(same_sample),
    median_same_sample_neighbor_fraction = median(same_sample),
    sampling_seed = 340L
  )
}

build_sequential_sample_tree <- function(number_of_samples) {
  if (number_of_samples < 2L) stop("Sequential sample tree requires at least two samples", call. = FALSE)
  sample_tree <- matrix(NA_integer_, nrow = number_of_samples - 1L, ncol = 2L)
  sample_tree[1L, ] <- c(-1L, -2L)
  if (number_of_samples > 2L) {
    for (step in 2L:(number_of_samples - 1L)) {
      sample_tree[step, ] <- c(step - 1L, -(step + 1L))
    }
  }
  sample_tree
}

run_rpca <- function(input_path, prefix, anchor_checkpoint_path) {
  input_object <- readRDS(input_path)
  if (uniqueN(input_object$geo_accession) != 40L || uniqueN(input_object$geo_subject) != 19L) {
    stop("Post-mitochondrial-QC object lost sample or subject traceability", call. = FALSE)
  }
  DefaultAssay(input_object) <- "RNA"
  object_list <- SplitObject(input_object, split.by = "geo_accession")
  if (length(object_list) != 40L) stop("RPCA split did not produce 40 sample objects", call. = FALSE)
  sample_names <- names(object_list)
  object_list <- lapply(sample_names, function(sample_name) {
    log_step(paste(prefix, "normalizing sample", sample_name))
    current <- object_list[[sample_name]]
    current <- NormalizeData(current, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    current <- FindVariableFeatures(current, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    current
  })
  names(object_list) <- sample_names
  integration_features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 2000)
  correction_features <- head(integration_features, integration_output_feature_count)
  object_list <- lapply(object_list, function(current) {
    subset(current, features = integration_features)
  })
  object_sizes <- vapply(object_list, ncol, numeric(1))
  sample_order <- names(sort(object_sizes, decreasing = TRUE))
  object_list <- object_list[sample_order]
  sample_tree <- build_sequential_sample_tree(length(object_list))
  write_tsv(data.table(
    integration_order = seq_along(sample_order),
    geo_accession = sample_order,
    cells = object_sizes[sample_order],
    sample_tree_dataset_index = -seq_along(sample_order)
  ), paste0("gse205506_paper_qc_", prefix, "_rpca_sequential_integration_order.tsv"))
  rm(input_object)
  invisible(gc())
  object_list <- lapply(object_list, function(current) {
    current <- ScaleData(current, features = integration_features, vars.to.regress = "nCount_RNA", verbose = FALSE)
    current <- RunPCA(current, features = integration_features, npcs = 30, approx = TRUE, seed.use = 340, verbose = FALSE)
    current
  })
  if (file.exists(anchor_checkpoint_path)) {
    log_step(paste(prefix, "loading light RPCA anchor checkpoint"))
    anchor_checkpoint <- readRDS(anchor_checkpoint_path)
    anchors <- anchor_checkpoint$anchors
    if (!identical(anchor_checkpoint$features, integration_features) || !identical(anchor_checkpoint$sample_order, sample_order)) {
      stop("RPCA anchor checkpoint does not match the current integration features or sample order", call. = FALSE)
    }
  } else {
    log_step(paste(prefix, "finding RPCA anchors"))
    anchors <- FindIntegrationAnchors(
      object.list = object_list,
      anchor.features = integration_features,
      reduction = "rpca",
      dims = 1:20,
      k.anchor = 5,
      verbose = TRUE
    )
    saveRDS(
      list(anchors = anchors, features = integration_features, sample_order = sample_order, sample_tree = sample_tree),
      anchor_checkpoint_path,
      compress = FALSE
    )
  }
  rm(object_list)
  invisible(gc())
  log_step(paste(prefix, "integrating data"))
  integrated <- IntegrateData(
    anchorset = anchors,
    dims = 1:20,
    features.to.integrate = correction_features,
    sample.tree = sample_tree,
    preserve.order = TRUE,
    verbose = TRUE
  )
  rm(anchors)
  invisible(gc())
  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated, vars.to.regress = "nCount_RNA", verbose = TRUE)
  integrated <- RunPCA(
    integrated,
    assay = "integrated",
    npcs = 30,
    approx = TRUE,
    seed.use = 340,
    reduction.name = paste0("pca.rpca.", prefix),
    reduction.key = paste0("RPCAPCA", toupper(prefix), "_"),
    verbose = TRUE
  )
  integrated <- RunUMAP(
    integrated,
    reduction = paste0("pca.rpca.", prefix),
    dims = 1:20,
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    seed.use = 340,
    reduction.name = paste0("umap.rpca.", prefix),
    reduction.key = paste0("RPCAUMAP", toupper(prefix), "_"),
    verbose = TRUE
  )
  integrated <- FindNeighbors(
    integrated,
    reduction = paste0("pca.rpca.", prefix),
    dims = 1:20,
    k.param = 20,
    graph.name = c(paste0("rpca_", prefix, "_nn"), paste0("rpca_", prefix, "_snn")),
    verbose = TRUE
  )
  if (file.exists(anchor_checkpoint_path)) unlink(anchor_checkpoint_path)
  list(
    object = integrated,
    features = integration_features,
    correction_features = correction_features,
    sample_order = sample_order,
    sample_tree = sample_tree
  )
}

if (!file.exists(input_path)) stop("Post-mitochondrial-QC raw object is missing", call. = FALSE)
object <- readRDS(input_path)
if (uniqueN(object$geo_accession) != 40L || uniqueN(object$geo_subject) != 19L) {
  stop("Post-mitochondrial-QC object lost sample or subject traceability", call. = FALSE)
}

log_step("Computing unintegrated RNA PCA for batch audit")
uncorrected <- object
DefaultAssay(uncorrected) <- "RNA"
uncorrected <- NormalizeData(uncorrected, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
uncorrected <- FindVariableFeatures(uncorrected, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
uncorrected <- ScaleData(
  uncorrected,
  features = VariableFeatures(uncorrected),
  vars.to.regress = "nCount_RNA",
  verbose = TRUE
)
uncorrected <- RunPCA(
  uncorrected,
  features = VariableFeatures(uncorrected),
  npcs = 20,
  approx = TRUE,
  seed.use = 340,
  reduction.name = "pca.uncorrected.final",
  reduction.key = "UNCPCA_",
  verbose = TRUE
)
uncorrected_embedding <- Embeddings(uncorrected, "pca.uncorrected.final")[, 1:20, drop = FALSE]
uncorrected_groups <- uncorrected$geo_accession
rm(uncorrected)
rm(object)
invisible(gc())

if (file.exists(light_checkpoint_path)) {
  log_step("Loading final light RPCA checkpoint")
  light_checkpoint <- readRDS(light_checkpoint_path)
  integrated_light <- light_checkpoint$object
  integration_features <- light_checkpoint$features
  correction_features <- light_checkpoint$correction_features
} else {
  log_step("Running final post-mitochondrial-QC memory-bounded RPCA integration")
  rpca_result <- run_rpca(input_path, "final", anchor_checkpoint_path)
  integrated_light <- rpca_result$object
  integration_features <- rpca_result$features
  correction_features <- rpca_result$correction_features
  integrated_light <- FindClusters(
    integrated_light,
    graph.name = "rpca_final_snn",
    resolution = 1.2,
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    random.seed = 340,
    cluster.name = "paper_final_cluster_r1_2",
    verbose = TRUE
  )
  integrated_light$paper_final_cluster_r1_2 <- as.character(integrated_light$paper_final_cluster_r1_2)
  saveRDS(
    list(object = integrated_light, features = integration_features, correction_features = correction_features),
    light_checkpoint_path,
    compress = FALSE
  )
}

log_step("Restoring all RNA genes after final light RPCA integration")
final_atlas <- readRDS(input_path)
DefaultAssay(final_atlas) <- "RNA"
final_atlas <- NormalizeData(final_atlas, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
target_cells <- colnames(final_atlas)
if (!setequal(target_cells, colnames(integrated_light))) {
  stop("Post-mitochondrial-QC and light RPCA objects do not contain the same cells", call. = FALSE)
}
final_atlas[["integrated"]] <- integrated_light[["integrated"]]
for (reduction_name in c("pca.rpca.final", "umap.rpca.final")) {
  final_atlas[[reduction_name]] <- integrated_light[[reduction_name]]
}
for (graph_name in c("rpca_final_nn", "rpca_final_snn")) {
  reordered_graph <- integrated_light[[graph_name]][target_cells, target_cells, drop = FALSE]
  final_atlas[[graph_name]] <- as(reordered_graph, "Graph")
  if (!identical(colnames(final_atlas[[graph_name]]), target_cells)) {
    stop("Final graph cell order does not match the full RNA object", call. = FALSE)
  }
}
final_atlas$paper_final_cluster_r1_2 <- integrated_light$paper_final_cluster_r1_2[match(colnames(final_atlas), colnames(integrated_light))]
if (anyNA(final_atlas$paper_final_cluster_r1_2)) stop("Final light RPCA cluster mapping failed", call. = FALSE)
rm(integrated_light)
invisible(gc())
final_atlas <- JoinLayers(final_atlas, assay = "RNA")
DefaultAssay(final_atlas) <- "RNA"
Idents(final_atlas) <- "paper_final_cluster_r1_2"
final_atlas[["pca.uncorrected.audit"]] <- CreateDimReducObject(
  embeddings = uncorrected_embedding[colnames(final_atlas), , drop = FALSE],
  key = "UNCPCA_",
  assay = "RNA"
)
saveRDS(final_atlas, final_path, compress = FALSE)
write_tsv(data.table(
  integration_feature = integration_features,
  anchor_feature = TRUE,
  integrated_correction_feature = integration_features %in% correction_features
), "gse205506_paper_qc_final_integration_features.tsv")
log_step("Saved final RPCA checkpoint before marker analysis")

log_step("Calculating final all-cell RNA markers")
marker_all <- as.data.table(wilcoxauc(
  final_atlas,
  group_by = "paper_final_cluster_r1_2",
  assay = "data",
  seurat_assay = "RNA"
))
setnames(marker_all, c("feature", "group"), c("gene", "cluster"))
marker_all[, `:=`(
  cluster = as.character(cluster),
  positive_direction = auc > 0.5 & logFC > 0,
  passes_primary_filter = padj <= 0.05 & auc > 0.5 & logFC > 0 & pct_in >= 5
)]
marker_all[, cluster_numeric := as.integer(cluster)]
setorder(marker_all, cluster_numeric, -passes_primary_filter, -positive_direction, -auc, -logFC, -pct_in, gene)
marker_all[, rank := seq_len(.N), by = cluster]
marker_all[, cluster_numeric := NULL]
top10 <- marker_all[rank <= 10L]
top50 <- marker_all[rank <= 50L]
write_tsv(top10, "gse205506_paper_qc_final_all_cell_top10_markers.tsv")
write_tsv(top50, "gse205506_paper_qc_final_all_cell_top50_markers.tsv")

marker_sets <- list(
  `T/I/NK` = c("CD3D", "CD3E", "TRAC", "TRBC1"),
  B = c("CD79A", "CD79B", "MS4A1", "TNFRSF17", "MZB1"),
  Myeloid = c("CD14", "CD68"),
  Epithelial = c("EPCAM", "CD24"),
  Fibroblast = c("COL1A2", "COL3A1", "MYH11", "ACTA2"),
  Endothelial = c("VWF", "PECAM1")
)
canonical_genes <- unique(unlist(marker_sets))
missing_markers <- setdiff(canonical_genes, rownames(final_atlas))
if (length(missing_markers)) stop("Required broad markers are absent: ", paste(missing_markers, collapse = "; "), call. = FALSE)
cluster_levels <- sort(unique(final_atlas$paper_final_cluster_r1_2))
rna_data <- LayerData(final_atlas, assay = "RNA", layer = "data")[canonical_genes, , drop = FALSE]
cluster_cells <- split(seq_len(ncol(final_atlas)), factor(final_atlas$paper_final_cluster_r1_2, levels = cluster_levels))
average_expression <- vapply(cluster_cells, function(index) Matrix::rowMeans(rna_data[, index, drop = FALSE]), numeric(length(canonical_genes)))
detection_fraction <- vapply(cluster_cells, function(index) Matrix::rowMeans(rna_data[, index, drop = FALSE] > 0), numeric(length(canonical_genes)))
rownames(average_expression) <- canonical_genes
rownames(detection_fraction) <- canonical_genes
colnames(average_expression) <- cluster_levels
colnames(detection_fraction) <- cluster_levels
gene_z <- t(scale(t(average_expression), center = TRUE, scale = TRUE))
gene_z[!is.finite(gene_z)] <- 0

score_rows <- rbindlist(lapply(names(marker_sets), function(compartment) {
  genes <- marker_sets[[compartment]]
  data.table(
    cluster = cluster_levels,
    compartment = compartment,
    signature_score = colMeans(gene_z[genes, , drop = FALSE]),
    mean_normalized_expression = colMeans(average_expression[genes, , drop = FALSE]),
    mean_detection_fraction = colMeans(detection_fraction[genes, , drop = FALSE]),
    marker_genes = paste(genes, collapse = ";")
  )
}))
score_rows[, cluster_numeric := as.integer(cluster)]
setorder(score_rows, cluster_numeric, -signature_score, compartment)
score_rows[, signature_rank := seq_len(.N), by = cluster]
score_rows[, cluster_numeric := NULL]
write_tsv(score_rows, "gse205506_paper_qc_final_broad_signature_scores.tsv")

metadata <- as.data.table(final_atlas[[]], keep.rownames = "cell_id")
cluster_distribution <- metadata[, .(
  cells = .N,
  samples = uniqueN(geo_accession),
  subjects = uniqueN(geo_subject)
), by = .(cluster = paper_final_cluster_r1_2)]
sample_counts <- metadata[, .N, by = .(cluster = paper_final_cluster_r1_2, geo_accession)]
subject_counts <- metadata[, .N, by = .(cluster = paper_final_cluster_r1_2, geo_subject)]
sample_max <- sample_counts[, .SD[which.max(N)], by = cluster]
subject_max <- subject_counts[, .SD[which.max(N)], by = cluster]
setnames(sample_max, c("geo_accession", "N"), c("maximum_sample", "maximum_sample_cells"))
setnames(subject_max, c("geo_subject", "N"), c("maximum_subject", "maximum_subject_cells"))
cluster_distribution <- merge(cluster_distribution, sample_max, by = "cluster", sort = FALSE)
cluster_distribution <- merge(cluster_distribution, subject_max, by = "cluster", sort = FALSE)
cluster_distribution[, `:=`(
  maximum_sample_fraction = maximum_sample_cells / cells,
  maximum_subject_fraction = maximum_subject_cells / cells
)]

annotation <- score_rows[signature_rank <= 2L, .(
  assigned_broad_compartment = compartment[[1L]],
  assigned_signature_score = signature_score[[1L]],
  second_broad_compartment = compartment[[2L]],
  second_signature_score = signature_score[[2L]],
  signature_score_margin = signature_score[[1L]] - signature_score[[2L]]
), by = cluster]
annotation <- merge(annotation, cluster_distribution, by = "cluster", sort = FALSE)
annotation[, `:=`(
  top10_markers = vapply(cluster, function(value) paste(top10[cluster == value][order(rank), gene], collapse = ";"), character(1)),
  top50_markers = vapply(cluster, function(value) paste(top50[cluster == value][order(rank), gene], collapse = ";"), character(1))
)]
annotation[, positive_marker_evidence := vapply(seq_len(nrow(annotation)), function(index) {
  genes <- marker_sets[[annotation$assigned_broad_compartment[[index]]]]
  paste(intersect(top50[cluster == annotation$cluster[[index]], gene], genes), collapse = ";")
}, character(1))]
annotation[, negative_or_conflicting_marker_evidence := vapply(seq_len(nrow(annotation)), function(index) {
  genes <- marker_sets[[annotation$second_broad_compartment[[index]]]]
  paste(intersect(top50[cluster == annotation$cluster[[index]], gene], genes), collapse = ";")
}, character(1))]
coexpression_rows <- rbindlist(lapply(seq_len(nrow(annotation)), function(index) {
  cluster_value <- annotation$cluster[[index]]
  cell_index <- which(final_atlas$paper_final_cluster_r1_2 == cluster_value)
  assigned_genes <- marker_sets[[annotation$assigned_broad_compartment[[index]]]]
  second_genes <- marker_sets[[annotation$second_broad_compartment[[index]]]]
  assigned_count <- Matrix::colSums(rna_data[assigned_genes, cell_index, drop = FALSE] > 0)
  second_count <- Matrix::colSums(rna_data[second_genes, cell_index, drop = FALSE] > 0)
  data.table(
    cluster = cluster_value,
    assigned_marker_positive_fraction = mean(assigned_count >= 1L),
    assigned_marker_coexpression_fraction_ge_2 = mean(assigned_count >= 2L),
    second_marker_positive_fraction = mean(second_count >= 1L),
    second_marker_coexpression_fraction_ge_2 = mean(second_count >= 2L)
  )
}))
annotation <- merge(annotation, coexpression_rows, by = "cluster", sort = FALSE)
annotation[, cluster_numeric := as.integer(cluster)]
setorder(annotation, cluster_numeric)
annotation[, cluster_numeric := NULL]
annotation[, annotation_method := "highest cluster-standardized source-paper broad marker signature with Top10/Top50, co-expression, conflicting evidence, sample and subject support"]
write_tsv(annotation, "gse205506_paper_qc_final_broad_cluster_evidence.tsv")

broad_map <- setNames(annotation$assigned_broad_compartment, annotation$cluster)
final_atlas$paper_final_broad_compartment <- unname(broad_map[final_atlas$paper_final_cluster_r1_2])
if (anyNA(final_atlas$paper_final_broad_compartment)) stop("Final broad annotation did not map to every cell", call. = FALSE)

integrated_embedding <- Embeddings(final_atlas, "pca.rpca.final")[, 1:20, drop = FALSE]
integrated_groups <- final_atlas$geo_accession
batch_knn <- rbindlist(list(
  knn_audit(uncorrected_embedding, uncorrected_groups, "pca_uncorrected_post_mt"),
  knn_audit(integrated_embedding, integrated_groups, "rpca_integrated_post_mt")
))
batch_eta <- rbindlist(list(
  eta_squared(uncorrected_embedding, uncorrected_groups, "pca_uncorrected_post_mt"),
  eta_squared(integrated_embedding, integrated_groups, "rpca_integrated_post_mt")
))
write_tsv(batch_knn, "gse205506_paper_qc_batch_knn_audit.tsv")
write_tsv(batch_eta, "gse205506_paper_qc_batch_pc_eta_squared.tsv")

metadata <- as.data.table(final_atlas[[]], keep.rownames = "cell_id")
count_dimensions <- list(
  project = character(),
  sample = "geo_accession",
  subject = "geo_subject",
  response = "table_s1_response",
  timepoint = "derived_timepoint_from_exact_geo_treatment",
  tissue = "geo_tissue",
  broad_compartment = "paper_final_broad_compartment"
)
final_counts <- rbindlist(lapply(names(count_dimensions), function(dimension_name) {
  column_name <- count_dimensions[[dimension_name]]
  if (!length(column_name)) {
    data.table(dimension = dimension_name, value = "GSE205506", cells = nrow(metadata))
  } else {
    metadata[, .(cells = .N), by = .(value = as.character(get(column_name)))][, dimension := dimension_name][]
  }
}), use.names = TRUE, fill = TRUE)
setcolorder(final_counts, c("dimension", "value", "cells"))
write_tsv(final_counts, "gse205506_paper_qc_final_counts.tsv")

cell_audit_path <- file.path(output_dir, "gse205506_paper_qc_cell_audit.tsv.gz")
cell_audit <- fread(cell_audit_path)
final_index <- match(metadata$cell_id, cell_audit$cell_id)
if (anyNA(final_index) || anyDuplicated(final_index)) stop("Cell audit join failed for final atlas", call. = FALSE)
cell_audit[final_index, `:=`(
  final_cluster_r1_2 = metadata$paper_final_cluster_r1_2,
  final_broad_compartment = metadata$paper_final_broad_compartment
)]
fwrite(cell_audit, cell_audit_path, sep = "\t", quote = FALSE, na = "", compress = "gzip")

final_atlas@misc$paper_final_atlas <- list(
  integration = "Seurat RPCA FindIntegrationAnchors reduction=rpca with 2000-feature light objects and dims 1:20; IntegrateData top-ranked 1000 features.to.integrate, dims 1:20, sequential sample.tree, preserve.order=TRUE; all RNA genes restored before markers",
  scaling_regression = "nCount_RNA",
  clustering = "integrated assay PCA 1:20, resolution 1.2",
  broad_compartments = names(marker_sets),
  batch_audit = "uncorrected RNA PCA versus RPCA-integrated PCA; KNN sample fraction and per-PC sample eta-squared"
)
saveRDS(final_atlas, final_path, compress = FALSE)

if (!requireNamespace("digest", quietly = TRUE)) stop("R package digest is required", call. = FALSE)
final_sha256 <- digest::digest(file = final_path, algo = "sha256", serialize = FALSE)
summary_table <- data.table(
  object_path = normalizePath(final_path, winslash = "/", mustWork = TRUE),
  size_bytes = as.numeric(file.info(final_path)$size),
  sha256 = final_sha256,
  features = nrow(final_atlas),
  cells = ncol(final_atlas),
  samples = uniqueN(final_atlas$geo_accession),
  subjects = uniqueN(final_atlas$geo_subject),
  final_clusters = uniqueN(final_atlas$paper_final_cluster_r1_2),
  broad_compartments = uniqueN(final_atlas$paper_final_broad_compartment),
  article_reported_cells = 155397L,
  difference_from_article = ncol(final_atlas) - 155397L,
  delivery = "local_only"
)
write_tsv(summary_table, "gse205506_paper_qc_final_object_summary.tsv")

package_names <- c("R", "Seurat", "SeuratObject", "Matrix", "SingleCellExperiment", "scDblFinder", "ggplot2", "presto", "RANN")
package_versions <- data.table(
  package = package_names,
  version = vapply(package_names, function(package_name) {
    if (package_name == "R") return(as.character(getRversion()))
    if (!requireNamespace(package_name, quietly = TRUE)) return(NA_character_)
    as.character(packageVersion(package_name))
  }, character(1))
)
write_tsv(package_versions, "gse205506_paper_qc_package_versions.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_gse205506_final_rpca_atlas.txt"))
if (file.exists(light_checkpoint_path)) unlink(light_checkpoint_path)
log_step("Final post-mitochondrial-QC RPCA atlas complete")
