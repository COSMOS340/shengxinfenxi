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
input_path <- file.path(analysis_dir, "gse205506_article_basic_qc_primary.rds")
first_round_path <- file.path(analysis_dir, "gse205506_first_round_rpca_broad_classification.rds")
post_mt_raw_path <- file.path(analysis_dir, "gse205506_post_compartment_mt_qc_raw.rds")
anchor_checkpoint_path <- file.path(analysis_dir, "gse205506_first_round_rpca_anchors_light_checkpoint.rds")
light_checkpoint_path <- file.path(analysis_dir, "gse205506_first_round_rpca_light_checkpoint.rds")
integration_output_feature_count <- 1000L

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

log_step <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
  flush.console()
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

if (!file.exists(input_path)) stop("Paper basic-QC object is missing", call. = FALSE)

run_rpca <- function(input_path, prefix, anchor_checkpoint_path) {
  input_object <- readRDS(input_path)
  if (ncol(input_object) != 200663L || uniqueN(input_object$geo_accession) != 40L || uniqueN(input_object$geo_subject) != 19L) {
    stop("Paper basic-QC object does not match the verified checkpoint", call. = FALSE)
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

if (file.exists(light_checkpoint_path)) {
  log_step("Loading first-round light RPCA checkpoint")
  light_checkpoint <- readRDS(light_checkpoint_path)
  integrated_light <- light_checkpoint$object
  integration_features <- light_checkpoint$features
  correction_features <- light_checkpoint$correction_features
} else {
  log_step("Running first-round memory-bounded RPCA integration")
  rpca_result <- run_rpca(input_path, "first", anchor_checkpoint_path)
  integrated_light <- rpca_result$object
  integration_features <- rpca_result$features
  correction_features <- rpca_result$correction_features
  integrated_light <- FindClusters(
    integrated_light,
    graph.name = "rpca_first_snn",
    resolution = 1.2,
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    random.seed = 340,
    cluster.name = "paper_first_round_cluster_r1_2",
    verbose = TRUE
  )
  integrated_light$paper_first_round_cluster_r1_2 <- as.character(integrated_light$paper_first_round_cluster_r1_2)
  saveRDS(
    list(object = integrated_light, features = integration_features, correction_features = correction_features),
    light_checkpoint_path,
    compress = FALSE
  )
}

log_step("Restoring all RNA genes after light RPCA integration")
first_round <- readRDS(input_path)
DefaultAssay(first_round) <- "RNA"
first_round <- NormalizeData(first_round, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
target_cells <- colnames(first_round)
if (!setequal(target_cells, colnames(integrated_light))) {
  stop("Full RNA and light RPCA objects do not contain the same cells", call. = FALSE)
}
first_round[["integrated"]] <- integrated_light[["integrated"]]
for (reduction_name in c("pca.rpca.first", "umap.rpca.first")) {
  first_round[[reduction_name]] <- integrated_light[[reduction_name]]
}
for (graph_name in c("rpca_first_nn", "rpca_first_snn")) {
  reordered_graph <- integrated_light[[graph_name]][target_cells, target_cells, drop = FALSE]
  first_round[[graph_name]] <- as(reordered_graph, "Graph")
  if (!identical(colnames(first_round[[graph_name]]), target_cells)) {
    stop("First-round graph cell order does not match the full RNA object", call. = FALSE)
  }
}
first_round$paper_first_round_cluster_r1_2 <- integrated_light$paper_first_round_cluster_r1_2[match(colnames(first_round), colnames(integrated_light))]
if (anyNA(first_round$paper_first_round_cluster_r1_2)) stop("Light RPCA cluster mapping failed", call. = FALSE)
rm(integrated_light)
invisible(gc())
first_round <- JoinLayers(first_round, assay = "RNA")
DefaultAssay(first_round) <- "RNA"
Idents(first_round) <- "paper_first_round_cluster_r1_2"
saveRDS(first_round, first_round_path, compress = FALSE)
write_tsv(data.table(
  integration_feature = integration_features,
  anchor_feature = TRUE,
  integrated_correction_feature = integration_features %in% correction_features
), "gse205506_paper_qc_first_round_integration_features.tsv")
log_step("Saved first-round RPCA checkpoint before marker analysis")

log_step("Calculating first-round full-cell RNA markers")
marker_all <- as.data.table(wilcoxauc(
  first_round,
  group_by = "paper_first_round_cluster_r1_2",
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
write_tsv(top10, "gse205506_paper_qc_first_round_top10_markers.tsv")
write_tsv(top50, "gse205506_paper_qc_first_round_top50_markers.tsv")

marker_sets <- list(
  `T/I/NK` = c("CD3D", "CD3E", "TRAC", "TRBC1"),
  B = c("CD79A", "CD79B", "MS4A1", "TNFRSF17", "MZB1"),
  Myeloid = c("CD14", "CD68"),
  Epithelial = c("EPCAM", "CD24"),
  Fibroblast = c("COL1A2", "COL3A1", "MYH11", "ACTA2"),
  Endothelial = c("VWF", "PECAM1")
)
missing_markers <- setdiff(unique(unlist(marker_sets)), rownames(first_round))
if (length(missing_markers)) stop("Required broad markers are absent: ", paste(missing_markers, collapse = "; "), call. = FALSE)

cluster_levels <- sort(unique(first_round$paper_first_round_cluster_r1_2))
canonical_genes <- unique(unlist(marker_sets))
rna_data <- LayerData(first_round, assay = "RNA", layer = "data")[canonical_genes, , drop = FALSE]
cluster_cells <- split(seq_len(ncol(first_round)), factor(first_round$paper_first_round_cluster_r1_2, levels = cluster_levels))
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
write_tsv(score_rows, "gse205506_paper_qc_first_round_broad_signature_scores.tsv")

metadata <- as.data.table(first_round[[]], keep.rownames = "cell_id")
cluster_distribution <- metadata[, .(
  cells = .N,
  samples = uniqueN(geo_accession),
  subjects = uniqueN(geo_subject)
), by = .(cluster = paper_first_round_cluster_r1_2)]
sample_counts <- metadata[, .N, by = .(cluster = paper_first_round_cluster_r1_2, geo_accession)]
subject_counts <- metadata[, .N, by = .(cluster = paper_first_round_cluster_r1_2, geo_subject)]
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
  assigned_mean_expression = mean_normalized_expression[[1L]],
  assigned_mean_detection_fraction = mean_detection_fraction[[1L]],
  second_broad_compartment = compartment[[2L]],
  second_signature_score = signature_score[[2L]],
  second_mean_expression = mean_normalized_expression[[2L]],
  second_mean_detection_fraction = mean_detection_fraction[[2L]],
  signature_score_margin = signature_score[[1L]] - signature_score[[2L]]
), by = cluster]
annotation <- merge(annotation, cluster_distribution, by = "cluster", sort = FALSE)
annotation[, `:=`(
  top10_markers = vapply(cluster, function(value) paste(top10[cluster == value][order(rank), gene], collapse = ";"), character(1)),
  top50_markers = vapply(cluster, function(value) paste(top50[cluster == value][order(rank), gene], collapse = ";"), character(1)),
  positive_marker_evidence = vapply(seq_len(.N), function(i) {
    genes <- marker_sets[[assigned_broad_compartment[[i]]]]
    paste(intersect(top50[cluster == annotation$cluster[[i]], gene], genes), collapse = ";")
  }, character(1)),
  negative_or_conflicting_marker_evidence = vapply(seq_len(.N), function(i) {
    genes <- marker_sets[[second_broad_compartment[[i]]]]
    paste(intersect(top50[cluster == annotation$cluster[[i]], gene], genes), collapse = ";")
  }, character(1))
)]

broad_map <- setNames(annotation$assigned_broad_compartment, annotation$cluster)
first_round$paper_first_round_broad_compartment <- unname(broad_map[first_round$paper_first_round_cluster_r1_2])
if (anyNA(first_round$paper_first_round_broad_compartment)) stop("Broad annotation did not map to every cell", call. = FALSE)

coexpression_rows <- rbindlist(lapply(seq_len(nrow(annotation)), function(i) {
  cluster_value <- annotation$cluster[[i]]
  cell_index <- which(first_round$paper_first_round_cluster_r1_2 == cluster_value)
  assigned_genes <- marker_sets[[annotation$assigned_broad_compartment[[i]]]]
  second_genes <- marker_sets[[annotation$second_broad_compartment[[i]]]]
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
annotation[, annotation_method := "highest cluster-standardized source-paper broad marker signature; Top10/Top50, co-expression, negative evidence, sample and subject support exported"]
write_tsv(annotation, "gse205506_paper_qc_broad_cluster_evidence.tsv")

hybrid_review <- annotation[, .(
  cluster,
  cells,
  samples,
  subjects,
  assigned_broad_compartment,
  second_broad_compartment,
  signature_score_margin,
  assigned_marker_coexpression_fraction_ge_2,
  second_marker_coexpression_fraction_ge_2,
  maximum_sample,
  maximum_sample_fraction,
  maximum_subject,
  maximum_subject_fraction,
  marker_conflict_review = signature_score_margin < 0.25 | second_marker_coexpression_fraction_ge_2 >= assigned_marker_coexpression_fraction_ge_2,
  abundance_rank_smallest_first = frank(cells, ties.method = "min"),
  abundance_review = "paper provides no numerical low-abundance exclusion threshold; rank exported for Mac review only",
  action = "retain cluster; no whole-cluster deletion in this request",
  top10_markers,
  positive_marker_evidence,
  negative_or_conflicting_marker_evidence
)]
write_tsv(hybrid_review, "gse205506_paper_qc_hybrid_review.tsv")

metadata <- as.data.table(first_round[[]], keep.rownames = "cell_id")
metadata[, paper_first_round_broad_compartment := as.character(paper_first_round_broad_compartment)]
metadata[, mt_model_group := fifelse(
  paper_first_round_broad_compartment %in% c("T/I/NK", "B"), "Lymphoid",
  paper_first_round_broad_compartment
)]
valid_model_groups <- c("Lymphoid", "Myeloid", "Fibroblast", "Endothelial", "Epithelial")
if (!setequal(unique(metadata$mt_model_group), valid_model_groups)) {
  stop("Unexpected mitochondrial model groups: ", paste(unique(metadata$mt_model_group), collapse = ";"), call. = FALSE)
}

metadata[, `:=`(
  mt_model_center = NA_real_,
  mt_model_sigma = NA_real_,
  mt_model_p_upper = NA_real_,
  mt_model_p_bonferroni = NA_real_,
  compartment_mt_qc_pass = FALSE
)]
model_rows <- list()
for (group_name in c("Lymphoid", "Myeloid", "Fibroblast", "Endothelial")) {
  index <- which(metadata$mt_model_group == group_name)
  center <- median(metadata$percent.mt[index], na.rm = TRUE)
  sigma <- mad(metadata$percent.mt[index], center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(sigma) || sigma == 0) stop("Non-finite or zero mitochondrial sigma for ", group_name, call. = FALSE)
  p_upper <- pnorm(metadata$percent.mt[index], mean = center, sd = sigma, lower.tail = FALSE)
  p_bonferroni <- p.adjust(p_upper, method = "bonferroni")
  pass <- p_bonferroni >= 0.05
  metadata[index, `:=`(
    mt_model_center = center,
    mt_model_sigma = sigma,
    mt_model_p_upper = p_upper,
    mt_model_p_bonferroni = p_bonferroni,
    compartment_mt_qc_pass = pass
  )]
  model_rows[[group_name]] <- data.table(
    model_group = group_name,
    broad_compartments = if (group_name == "Lymphoid") "T/I/NK;B" else group_name,
    cells = length(index),
    center = center,
    sigma = sigma,
    removed_cells = sum(!pass),
    removal_fraction = mean(!pass),
    rule = "one-sided normal upper tail using median center, MAD constant 1.4826, Bonferroni-adjusted p < 0.05"
  )
}
epithelial_index <- which(metadata$mt_model_group == "Epithelial")
epithelial_pass <- metadata$percent.mt[epithelial_index] <= 75
metadata[epithelial_index, `:=`(
  mt_model_center = 75,
  mt_model_sigma = NA_real_,
  mt_model_p_upper = NA_real_,
  mt_model_p_bonferroni = NA_real_,
  compartment_mt_qc_pass = epithelial_pass
)]
model_rows[["Epithelial"]] <- data.table(
  model_group = "Epithelial",
  broad_compartments = "Epithelial",
  cells = length(epithelial_index),
  center = 75,
  sigma = NA_real_,
  removed_cells = sum(!epithelial_pass),
  removal_fraction = mean(!epithelial_pass),
  rule = "retain percent.mt <= 75; remove percent.mt > 75"
)
mt_model <- rbindlist(model_rows, use.names = TRUE, fill = TRUE)
write_tsv(mt_model, "gse205506_paper_qc_compartment_mt_model.tsv")

reported_fraction <- c(Lymphoid = 0.0920, Myeloid = 0.1284, Fibroblast = 0.0811, Endothelial = 0.0850, Epithelial = 0.2975)
removal_comparison <- copy(mt_model)
removal_comparison[, article_reported_removal_fraction := unname(reported_fraction[model_group])]
removal_comparison[, absolute_difference := removal_fraction - article_reported_removal_fraction]
removal_comparison[, comparison_note := "comparison target only; no parameter tuning or downsampling used"]
write_tsv(removal_comparison, "gse205506_paper_qc_removal_fraction_comparison.tsv")

cell_audit_path <- file.path(output_dir, "gse205506_paper_qc_cell_audit.tsv.gz")
cell_audit <- fread(cell_audit_path)
audit_derived_columns <- c(
  "first_round_broad_compartment",
  "mt_model_group",
  "mt_model_center",
  "mt_model_sigma",
  "mt_model_p_upper",
  "mt_model_p_bonferroni",
  "compartment_mt_qc_pass",
  "formal_final_inclusion"
)
existing_audit_derived_columns <- intersect(audit_derived_columns, names(cell_audit))
if (length(existing_audit_derived_columns)) {
  cell_audit[, (existing_audit_derived_columns) := NULL]
}
cell_audit[, `:=`(
  first_round_broad_compartment = NA_character_,
  mt_model_group = NA_character_,
  mt_model_center = NA_real_,
  mt_model_sigma = NA_real_,
  mt_model_p_upper = NA_real_,
  mt_model_p_bonferroni = NA_real_,
  compartment_mt_qc_pass = NA,
  formal_final_inclusion = FALSE
)]
basic_index <- match(metadata$cell_id, cell_audit$cell_id)
if (anyNA(basic_index) || anyDuplicated(basic_index)) stop("Cell audit join failed for first-round object", call. = FALSE)
cell_audit[basic_index, `:=`(
  first_round_broad_compartment = metadata$paper_first_round_broad_compartment,
  mt_model_group = metadata$mt_model_group,
  mt_model_center = metadata$mt_model_center,
  mt_model_sigma = metadata$mt_model_sigma,
  mt_model_p_upper = metadata$mt_model_p_upper,
  mt_model_p_bonferroni = metadata$mt_model_p_bonferroni,
  compartment_mt_qc_pass = metadata$compartment_mt_qc_pass,
  formal_final_inclusion = metadata$compartment_mt_qc_pass
)]
fwrite(cell_audit, cell_audit_path, sep = "\t", quote = FALSE, na = "", compress = "gzip")

final_cell_ids <- metadata[compartment_mt_qc_pass == TRUE, .(
  cell_id,
  geo_accession,
  geo_subject,
  geo_tissue,
  geo_treatment,
  derived_timepoint_from_exact_geo_treatment,
  table_s1_response,
  paper_first_round_cluster_r1_2,
  paper_first_round_broad_compartment,
  percent.mt,
  mt_model_group,
  mt_model_p_bonferroni
)]
fwrite(final_cell_ids, file.path(output_dir, "gse205506_paper_qc_final_cell_ids.tsv.gz"), sep = "\t", quote = FALSE, na = "", compress = "gzip")

sample_counts_path <- file.path(output_dir, "gse205506_paper_qc_sample_counts.tsv")
sample_counts_table <- fread(sample_counts_path)
mt_sample_columns <- c(
  "cells_entered_compartment_mt_qc",
  "cells_pass_compartment_mt_qc",
  "cells_fail_compartment_mt_qc",
  "compartment_mt_retention_fraction"
)
existing_mt_sample_columns <- intersect(mt_sample_columns, names(sample_counts_table))
if (length(existing_mt_sample_columns)) {
  sample_counts_table[, (existing_mt_sample_columns) := NULL]
}
mt_sample <- metadata[, .(
  cells_entered_compartment_mt_qc = .N,
  cells_pass_compartment_mt_qc = sum(compartment_mt_qc_pass),
  cells_fail_compartment_mt_qc = sum(!compartment_mt_qc_pass)
), by = geo_accession]
sample_counts_table <- merge(sample_counts_table, mt_sample, by = "geo_accession", all.x = TRUE, sort = FALSE)
sample_counts_table[, compartment_mt_retention_fraction := cells_pass_compartment_mt_qc / cells_entered_compartment_mt_qc]
write_tsv(sample_counts_table, "gse205506_paper_qc_sample_counts.tsv")

stage_counts_path <- file.path(output_dir, "gse205506_paper_qc_stage_counts.tsv")
stage_counts <- fread(stage_counts_path)
stage_counts <- stage_counts[!stage %in% c("compartment_mt_qc_pass", "compartment_mt_qc_fail")]
stage_counts <- rbind(
  stage_counts,
  data.table(
    stage = c("compartment_mt_qc_pass", "compartment_mt_qc_fail"),
    cells = c(sum(metadata$compartment_mt_qc_pass), sum(!metadata$compartment_mt_qc_pass)),
    primary_route = TRUE,
    notes = c("paper basic-QC pass cells retained by compartment-aware mitochondrial QC", "paper basic-QC pass cells excluded by compartment-aware mitochondrial QC")
  ),
  fill = TRUE
)
write_tsv(stage_counts, "gse205506_paper_qc_stage_counts.tsv")

doublet_path <- file.path(output_dir, "gse205506_paper_qc_doublet_sensitivity.tsv")
doublet_sample <- fread(doublet_path)
if ("record_type" %in% names(doublet_sample)) {
  doublet_sample <- doublet_sample[record_type == "per_sample_call_summary"]
}
doublet_sample[, record_type := "per_sample_call_summary"]
doublet_composition <- rbindlist(lapply(c("all_primary_cells", "exclude_scDblFinder_doublets_sensitivity"), function(filter_state) {
  current <- if (filter_state == "all_primary_cells") metadata else metadata[scDblFinder.class.sensitivity == "singlet"]
  composition <- current[, .(cells = .N), by = .(
    paper_first_round_cluster_r1_2,
    paper_first_round_broad_compartment
  )]
  composition[, `:=`(
    record_type = "broad_cluster_composition",
    filter_state = filter_state
  )]
  composition[, proportion_within_filter_state := cells / sum(cells)]
  composition
}), use.names = TRUE, fill = TRUE)
doublet_output <- rbindlist(list(doublet_sample, doublet_composition), use.names = TRUE, fill = TRUE)
write_tsv(doublet_output, "gse205506_paper_qc_doublet_sensitivity.tsv")

first_round$paper_mt_model_group <- metadata$mt_model_group[match(colnames(first_round), metadata$cell_id)]
first_round$paper_mt_model_center <- metadata$mt_model_center[match(colnames(first_round), metadata$cell_id)]
first_round$paper_mt_model_sigma <- metadata$mt_model_sigma[match(colnames(first_round), metadata$cell_id)]
first_round$paper_mt_p_upper <- metadata$mt_model_p_upper[match(colnames(first_round), metadata$cell_id)]
first_round$paper_mt_p_bonferroni <- metadata$mt_model_p_bonferroni[match(colnames(first_round), metadata$cell_id)]
first_round$paper_compartment_mt_qc_pass <- metadata$compartment_mt_qc_pass[match(colnames(first_round), metadata$cell_id)]
first_round@misc$paper_first_round <- list(
  integration = "Seurat RPCA FindIntegrationAnchors reduction=rpca with 2000-feature light objects and dims 1:20; IntegrateData top-ranked 1000 features.to.integrate, dims 1:20, sequential sample.tree, preserve.order=TRUE; all 27499 RNA genes restored before markers",
  scaling_regression = "nCount_RNA",
  clustering = "integrated assay PCA 1:20, resolution 1.2",
  broad_compartments = names(marker_sets),
  mitochondrial_qc = "epithelial percent.mt<=75; combined lymphoid, myeloid, fibroblast, endothelial one-sided median/MAD normal model with Bonferroni p>=0.05"
)
saveRDS(first_round, first_round_path, compress = FALSE)

object <- readRDS(input_path)
post_mt_raw <- subset(object, cells = final_cell_ids$cell_id)
post_mt_raw$paper_first_round_cluster_r1_2 <- first_round$paper_first_round_cluster_r1_2[match(colnames(post_mt_raw), colnames(first_round))]
post_mt_raw$paper_first_round_broad_compartment <- first_round$paper_first_round_broad_compartment[match(colnames(post_mt_raw), colnames(first_round))]
post_mt_raw$paper_mt_model_group <- first_round$paper_mt_model_group[match(colnames(post_mt_raw), colnames(first_round))]
post_mt_raw$paper_mt_p_bonferroni <- first_round$paper_mt_p_bonferroni[match(colnames(post_mt_raw), colnames(first_round))]
saveRDS(post_mt_raw, post_mt_raw_path, compress = FALSE)
rm(object)

if (!requireNamespace("digest", quietly = TRUE)) stop("R package digest is required", call. = FALSE)
object_inventory <- rbindlist(lapply(c(first_round_path, post_mt_raw_path), function(path) {
  data.table(
    object_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    size_bytes = as.numeric(file.info(path)$size),
    sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE),
    cells = if (path == first_round_path) ncol(first_round) else ncol(post_mt_raw),
    delivery = "local_only"
  )
}))
write_tsv(object_inventory, "gse205506_paper_qc_first_round_object_inventory.tsv")

parameters_path <- file.path(output_dir, "gse205506_paper_qc_parameters.tsv")
parameters <- fread(parameters_path)
parameters <- rbind(
  parameters,
  data.table(
    parameter = c(
      "rpca_variable_features", "rpca_anchor_reduction", "rpca_anchor_dims", "integrated_scale_regress",
      "first_round_pca_dims", "first_round_resolution", "first_round_umap_neighbors", "first_round_umap_min_dist",
      "rpca_sample_tree", "rpca_preserve_order", "rpca_light_object_features",
      "rpca_features_to_integrate",
      "epithelial_mt_threshold", "non_epithelial_mt_model", "mt_adjustment"
    ),
    value = c(2000, "rpca", "1:20", "nCount_RNA", "1:20", 1.2, 30, 0.3, "sequential accumulated object plus one sample per step", TRUE, 2000, 1000, 75, "median center; MAD constant 1.4826 normal upper tail", "bonferroni p>=0.05 retained")
  ),
  fill = TRUE
)
write_tsv(parameters, "gse205506_paper_qc_parameters.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_gse205506_first_round_rpca_mt_qc.txt"))
if (file.exists(light_checkpoint_path)) unlink(light_checkpoint_path)
log_step("First-round RPCA broad classification and compartment mitochondrial QC complete")
