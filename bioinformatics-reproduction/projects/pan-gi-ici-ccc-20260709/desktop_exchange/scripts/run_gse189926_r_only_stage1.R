suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(harmony)
  library(RANN)
})

options(stringsAsFactors = FALSE, future.globals.maxSize = 16 * 1024^3)
set.seed(340)

project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
exchange_dir <- file.path(project_dir, "desktop_exchange")
matrix_inventory_path <- file.path(exchange_dir, "uploads/20260709_matrix_download/matrix_file_inventory.tsv")
sample_qc_path <- file.path(exchange_dir, "uploads/20260709_preprocess_priority_response/gse189926_sample_qc.tsv")
sample_parse_path <- file.path(exchange_dir, "uploads/20260709_preprocess_priority_response/gse189926_sample_parse_audit.tsv")
out_dir <- file.path(exchange_dir, "uploads/20260710_gse189926_r_only_rerun")
object_dir <- file.path(project_dir, "04_objects/20260710_gse189926_r_only_rerun")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "")
}

log_step <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
  flush.console()
}

read_sparse_text_matrix <- function(path, sample_accession, chunk_lines = 200L) {
  log_step(sprintf("Importing %s", basename(path)))
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L, warn = FALSE)
  if (length(header) != 1L) stop(sprintf("Missing header in %s", path))
  header_fields <- strsplit(header, "\t", fixed = TRUE)[[1L]]
  if (length(header_fields) < 2L || nzchar(header_fields[[1L]])) {
    stop(sprintf("Unexpected first header field in %s", path))
  }
  raw_cell_ids <- header_fields[-1L]
  n_cells <- length(raw_cell_ids)
  genes <- list()
  i_parts <- list()
  j_parts <- list()
  x_parts <- list()
  row_offset <- 0L
  part <- 0L

  repeat {
    lines <- readLines(con, n = chunk_lines, warn = FALSE)
    if (!length(lines)) break
    part <- part + 1L
    dt <- fread(
      text = paste(lines, collapse = "\n"),
      sep = "\t",
      header = FALSE,
      data.table = FALSE,
      showProgress = FALSE,
      check.names = FALSE
    )
    if (ncol(dt) != n_cells + 1L) {
      stop(sprintf(
        "Column count mismatch in %s at rows %d-%d: expected %d, observed %d",
        path, row_offset + 1L, row_offset + nrow(dt), n_cells + 1L, ncol(dt)
      ))
    }
    chunk_genes <- as.character(dt[[1L]])
    if (anyNA(chunk_genes) || any(!nzchar(chunk_genes))) {
      stop(sprintf("Missing feature name in %s at rows %d-%d", path, row_offset + 1L, row_offset + nrow(dt)))
    }
    values <- as.matrix(dt[-1L])
    storage.mode(values) <- "double"
    if (anyNA(values)) {
      stop(sprintf("Non-numeric expression value in %s at rows %d-%d", path, row_offset + 1L, row_offset + nrow(dt)))
    }
    nz <- which(values != 0, arr.ind = TRUE)
    genes[[part]] <- chunk_genes
    if (nrow(nz)) {
      i_parts[[part]] <- as.integer(row_offset + nz[, 1L])
      j_parts[[part]] <- as.integer(nz[, 2L])
      x_parts[[part]] <- values[nz]
    } else {
      i_parts[[part]] <- integer()
      j_parts[[part]] <- integer()
      x_parts[[part]] <- numeric()
    }
    row_offset <- row_offset + nrow(dt)
    rm(dt, values, nz)
  }

  all_genes <- unlist(genes, use.names = FALSE)
  unique_genes <- unique(all_genes)
  duplicate_feature_rows <- length(all_genes) - length(unique_genes)
  row_map <- match(all_genes, unique_genes)
  original_i <- unlist(i_parts, use.names = FALSE)
  matrix_i <- row_map[original_i]
  matrix_j <- unlist(j_parts, use.names = FALSE)
  matrix_x <- unlist(x_parts, use.names = FALSE)
  cell_ids <- paste(sample_accession, raw_cell_ids, sep = "_")
  mat <- sparseMatrix(
    i = matrix_i,
    j = matrix_j,
    x = matrix_x,
    dims = c(length(unique_genes), n_cells),
    dimnames = list(unique_genes, cell_ids),
    giveCsparse = TRUE
  )
  storage.mode(mat@x) <- "double"
  list(
    matrix = mat,
    observed_feature_rows = length(all_genes),
    unique_features = length(unique_genes),
    duplicate_feature_rows_collapsed = duplicate_feature_rows,
    cells = n_cells,
    nonzero_values = length(mat@x)
  )
}

align_to_union <- function(mat, union_features) {
  feature_map <- match(rownames(mat), union_features)
  if (anyNA(feature_map)) stop("Feature union alignment failed")
  out <- sparseMatrix(
    i = feature_map[mat@i + 1L],
    j = rep.int(seq_len(ncol(mat)), diff(mat@p)),
    x = mat@x,
    dims = c(length(union_features), ncol(mat)),
    dimnames = list(union_features, colnames(mat)),
    giveCsparse = TRUE
  )
  storage.mode(out@x) <- "double"
  out
}

technical_feature_flag <- function(features) {
  grepl("^MT-", features) |
    grepl("^RP[SL]", features) |
    grepl("^HB[ABDEGMQZ][0-9]*$", features) |
    grepl("^MALAT1($|-)", features) |
    grepl("^(IGHA|IGHG|IGHM|IGHD|IGHE|IGKC|IGLC)", features)
}

percent_for_features <- function(counts, feature_index) {
  totals <- Matrix::colSums(counts)
  numerator <- if (any(feature_index)) Matrix::colSums(counts[feature_index, , drop = FALSE]) else rep(0, ncol(counts))
  ifelse(totals > 0, 100 * numerator / totals, 0)
}

run_log_normalized_route <- function(object, reduction_prefix, run_clusters = FALSE) {
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  hvg_initial <- VariableFeatures(object)
  hvg_technical <- technical_feature_flag(hvg_initial)
  pca_features <- hvg_initial[!hvg_technical]
  if (length(pca_features) < 500L) stop("Fewer than 500 non-technical variable features remained for PCA")
  VariableFeatures(object) <- pca_features
  object <- ScaleData(object, features = pca_features, verbose = FALSE)
  pca_name <- paste0("pca.", reduction_prefix)
  umap_name <- paste0("umap.", reduction_prefix)
  nn_name <- paste0(reduction_prefix, "_nn")
  snn_name <- paste0(reduction_prefix, "_snn")
  object <- RunPCA(object, features = pca_features, npcs = 40, reduction.name = pca_name, verbose = FALSE)
  object <- FindNeighbors(
    object,
    reduction = pca_name,
    dims = 1:30,
    k.param = 15,
    graph.name = c(nn_name, snn_name),
    verbose = FALSE
  )
  object <- RunUMAP(
    object,
    reduction = pca_name,
    dims = 1:30,
    n.neighbors = 15,
    min.dist = 0.3,
    seed.use = 340,
    reduction.name = umap_name,
    verbose = FALSE
  )
  if (run_clusters) {
    for (res in c(0.3, 0.5, 0.8)) {
      object <- FindClusters(
        object,
        graph.name = snn_name,
        resolution = res,
        algorithm = 1,
        random.seed = 340,
        verbose = FALSE
      )
      field <- switch(as.character(res), "0.3" = "cluster_qc_r0_3", "0.5" = "cluster_qc_r0_5", "0.8" = "cluster_qc_r0_8")
      object[[field]] <- as.character(Idents(object))
    }
  }
  attr(object, "route_audit") <- list(
    hvg_initial = length(hvg_initial),
    technical_hvg_excluded = sum(hvg_technical),
    pca_features = length(pca_features),
    pca_feature_names = pca_features,
    pca_name = pca_name,
    umap_name = umap_name,
    nn_name = nn_name,
    snn_name = snn_name
  )
  object
}

mixing_metrics <- function(coords, metadata, route_name, k = 15L) {
  nn <- RANN::nn2(coords, k = k + 1L)$nn.idx[, -1L, drop = FALSE]
  same_sample <- rowMeans(matrix(metadata$sample_accession[nn], nrow = nrow(nn)) == metadata$sample_accession)
  same_patient <- rowMeans(matrix(metadata$patient_id[nn], nrow = nrow(nn)) == metadata$patient_id)
  data.frame(
    route_name = route_name,
    k_neighbors = k,
    cells_used = nrow(coords),
    mean_same_sample_neighbor_fraction = mean(same_sample),
    median_same_sample_neighbor_fraction = median(same_sample),
    mean_same_patient_neighbor_fraction = mean(same_patient),
    median_same_patient_neighbor_fraction = median(same_patient),
    interpretation = "Lower same-sample fraction indicates less sample-dominated local structure; patient structure is retained as a separate diagnostic.",
    stringsAsFactors = FALSE
  )
}

marker_tables <- function(object, cluster_field) {
  data_mat <- LayerData(object, assay = "RNA", layer = "data")
  clusters <- as.character(object[[cluster_field, drop = TRUE]])
  cluster_levels <- sort(unique(clusters), na.last = TRUE)
  cluster_index <- match(clusters, cluster_levels)
  membership <- sparseMatrix(
    i = seq_along(cluster_index),
    j = cluster_index,
    x = 1,
    dims = c(length(cluster_index), length(cluster_levels))
  )
  cluster_sizes <- as.numeric(table(factor(clusters, levels = cluster_levels)))
  expression_sums <- data_mat %*% membership
  detection_sums <- (data_mat > 0) %*% membership
  total_expression <- Matrix::rowSums(data_mat)
  total_detection <- Matrix::rowSums(data_mat > 0)
  total_cells <- ncol(data_mat)
  genes <- rownames(data_mat)
  technical <- technical_feature_flag(genes)
  raw_rows <- list()
  clean_rows <- list()

  for (idx in seq_along(cluster_levels)) {
    n_in <- cluster_sizes[[idx]]
    n_out <- total_cells - n_in
    avg_in <- as.numeric(expression_sums[, idx]) / n_in
    avg_out <- (total_expression - as.numeric(expression_sums[, idx])) / n_out
    pct_in <- as.numeric(detection_sums[, idx]) / n_in
    pct_out <- (total_detection - as.numeric(detection_sums[, idx])) / n_out
    difference <- avg_in - avg_out
    base <- data.frame(
      cluster = cluster_levels[[idx]],
      gene = genes,
      avg_log_normalized_expression_in = avg_in,
      avg_log_normalized_expression_out = avg_out,
      avg_log_expression_difference = difference,
      pct_expressing_in = pct_in,
      pct_expressing_out = pct_out,
      cluster_cell_count = n_in,
      technical_gene = technical,
      stringsAsFactors = FALSE
    )
    eligible <- base$pct_expressing_in >= 0.05 & base$avg_log_expression_difference > 0
    raw <- base[eligible, , drop = FALSE]
    raw <- raw[order(-raw$avg_log_expression_difference, -raw$pct_expressing_in, raw$gene), , drop = FALSE]
    raw <- head(raw, 50L)
    raw$rank <- seq_len(nrow(raw))
    clean <- base[eligible & !base$technical_gene, , drop = FALSE]
    clean <- clean[order(-clean$avg_log_expression_difference, -clean$pct_expressing_in, clean$gene), , drop = FALSE]
    clean <- head(clean, 50L)
    clean$rank <- seq_len(nrow(clean))
    raw_rows[[idx]] <- raw[, c("cluster", "rank", setdiff(names(raw), c("cluster", "rank")))]
    clean_rows[[idx]] <- clean[, c("cluster", "rank", setdiff(names(clean), c("cluster", "rank")))]
  }
  list(raw = rbindlist(raw_rows), no_technical = rbindlist(clean_rows))
}

inventory <- fread(matrix_inventory_path, sep = "\t", data.table = FALSE, check.names = FALSE)
inventory <- inventory[inventory$dataset_id == "GSE189926", , drop = FALSE]
sample_qc <- fread(sample_qc_path, sep = "\t", data.table = FALSE, check.names = FALSE)
sample_parse <- fread(sample_parse_path, sep = "\t", data.table = FALSE, check.names = FALSE)

if (nrow(inventory) != 22L) stop(sprintf("Expected 22 GSE189926 matrix files, observed %d", nrow(inventory)))
if (!identical(sort(inventory$file_name), sort(sample_qc$file_name))) stop("Matrix inventory and sample metadata file names do not match exactly")
if (!all(sample_parse$parse_status == "parsed")) stop("Sample parse audit contains a non-parsed row")
if (!all(file.exists(inventory$local_path))) stop("One or more matrix inventory paths do not exist")

required_metadata_fields <- c(
  "sample_accession", "sample_title", "patient_id", "timepoint",
  "characteristics_ch1::outcome", "characteristics_ch1::treatment",
  "characteristics_ch1::mmr", "characteristics_ch1::cell type",
  "characteristics_ch1::tissue"
)
missing_metadata_fields <- setdiff(required_metadata_fields, names(sample_qc))
if (length(missing_metadata_fields)) stop(sprintf("Missing metadata fields: %s", paste(missing_metadata_fields, collapse = ", ")))

matrix_list <- vector("list", nrow(inventory))
import_rows <- vector("list", nrow(inventory))
for (idx in seq_len(nrow(inventory))) {
  file_row <- inventory[idx, , drop = FALSE]
  meta_idx <- match(file_row$file_name, sample_qc$file_name)
  sample_accession <- sample_qc$sample_accession[[meta_idx]]
  imported <- read_sparse_text_matrix(file_row$local_path, sample_accession)
  matrix_list[[idx]] <- imported$matrix
  import_rows[[idx]] <- data.frame(
    sample_accession = sample_accession,
    file_name = file_row$file_name,
    local_path = file_row$local_path,
    file_size_bytes = as.numeric(file_row$file_size_bytes),
    sha256_from_inventory = file_row$sha256,
    import_status = "completed",
    import_method = "R chunked gzip text reader to Matrix dgCMatrix; no full matrix dense expansion",
    observed_feature_rows = imported$observed_feature_rows,
    unique_features = imported$unique_features,
    duplicate_feature_rows_collapsed = imported$duplicate_feature_rows_collapsed,
    cells_imported = imported$cells,
    nonzero_values = imported$nonzero_values,
    stringsAsFactors = FALSE
  )
  rm(imported)
  gc(verbose = FALSE)
}

union_features <- unique(unlist(lapply(matrix_list, rownames), use.names = FALSE))
feature_alignment <- rbindlist(lapply(seq_along(matrix_list), function(idx) {
  data.frame(
    sample_accession = import_rows[[idx]]$sample_accession,
    observed_unique_features = nrow(matrix_list[[idx]]),
    union_features = length(union_features),
    missing_features_filled_with_zero = length(union_features) - nrow(matrix_list[[idx]]),
    alignment_status = "completed",
    alignment_method = "union of exact feature names; missing rows represented as sparse zeros",
    stringsAsFactors = FALSE
  )
}))
write_tsv(rbindlist(import_rows), file.path(out_dir, "gse189926_r_import_audit.tsv"))
write_tsv(feature_alignment, file.path(out_dir, "gse189926_r_feature_alignment_audit.tsv"))

log_step(sprintf("Aligning %d features across %d samples", length(union_features), length(matrix_list)))
aligned_list <- lapply(matrix_list, align_to_union, union_features = union_features)
counts <- do.call(cbind, aligned_list)
rm(matrix_list, aligned_list)
gc(verbose = FALSE)

cell_sample <- sub("_.*$", "", colnames(counts))
meta_index <- match(cell_sample, sample_qc$sample_accession)
if (anyNA(meta_index)) stop("At least one prefixed cell ID does not map to sample_accession")
metadata <- sample_qc[meta_index, required_metadata_fields, drop = FALSE]
rownames(metadata) <- colnames(counts)
metadata$timepoint_simple <- metadata$timepoint
metadata$timepoint_simple[metadata$timepoint %in% c("post-treatment1", "post-treatment2")] <- "post-treatment"
metadata$binary_response_group <- "not_in_binary_response"
metadata$binary_response_group[metadata[["characteristics_ch1::outcome"]] == "PR"] <- "responder"
metadata$binary_response_group[metadata[["characteristics_ch1::outcome"]] %in% c("SD", "PD")] <- "non_responder"

log_step(sprintf("Creating Seurat object with %d cells and %d features", ncol(counts), nrow(counts)))
object <- CreateSeuratObject(counts = counts, project = "GSE189926", assay = "RNA", meta.data = metadata, min.cells = 0, min.features = 0)
rm(counts, metadata)
gc(verbose = FALSE)

raw_counts <- LayerData(object, assay = "RNA", layer = "counts")
features <- rownames(raw_counts)
object$pct_mt <- percent_for_features(raw_counts, grepl("^MT-", features))
object$pct_ribo <- percent_for_features(raw_counts, grepl("^RP[SL]", features))
object$pct_hb <- percent_for_features(raw_counts, grepl("^HB[ABDEGMQZ][0-9]*$", features))

metric_quantiles <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  c(q1 = q[[1L]], q3 = q[[2L]], iqr = q[[2L]] - q[[1L]])
}
q_gene <- metric_quantiles(object$nFeature_RNA)
q_count <- metric_quantiles(object$nCount_RNA)
q_mt <- metric_quantiles(object$pct_mt)
q_ribo <- metric_quantiles(object$pct_ribo)
q_hb <- metric_quantiles(object$pct_hb)

thresholds <- c(
  low_detected_genes = max(200, q_gene[["q1"]] - 1.5 * q_gene[["iqr"]]),
  high_detected_genes = q_gene[["q3"]] + 3 * q_gene[["iqr"]],
  high_total_counts = q_count[["q3"]] + 3 * q_count[["iqr"]],
  high_mitochondrial_percentage = min(25, max(10, q_mt[["q3"]] + 1.5 * q_mt[["iqr"]])),
  high_ribosomal_percentage = min(50, max(20, q_ribo[["q3"]] + 3 * q_ribo[["iqr"]])),
  high_hemoglobin_percentage = max(5, q_hb[["q3"]] + 3 * q_hb[["iqr"]])
)

object$qc_low_detected_genes <- object$nFeature_RNA < thresholds[["low_detected_genes"]]
object$qc_high_detected_genes <- object$nFeature_RNA > thresholds[["high_detected_genes"]]
object$qc_high_total_counts <- object$nCount_RNA > thresholds[["high_total_counts"]]
object$qc_high_mitochondrial_percentage <- object$pct_mt > thresholds[["high_mitochondrial_percentage"]]
object$qc_high_ribosomal_percentage <- object$pct_ribo > thresholds[["high_ribosomal_percentage"]]
object$qc_high_hemoglobin_percentage <- object$pct_hb > thresholds[["high_hemoglobin_percentage"]]
flag_fields <- c(
  "qc_low_detected_genes", "qc_high_detected_genes", "qc_high_total_counts",
  "qc_high_mitochondrial_percentage", "qc_high_ribosomal_percentage", "qc_high_hemoglobin_percentage"
)
object$qc_pass <- rowSums(object[[]][, flag_fields, drop = FALSE]) == 0L

threshold_rule <- c(
  low_detected_genes = "nFeature_RNA < max(200, Q1 - 1.5*IQR)",
  high_detected_genes = "nFeature_RNA > Q3 + 3*IQR",
  high_total_counts = "nCount_RNA > Q3 + 3*IQR",
  high_mitochondrial_percentage = "pct_mt > min(25, max(10, Q3 + 1.5*IQR))",
  high_ribosomal_percentage = "pct_ribo > min(50, max(20, Q3 + 3*IQR))",
  high_hemoglobin_percentage = "pct_hb > max(5, Q3 + 3*IQR)"
)
metric_name <- c(
  low_detected_genes = "nFeature_RNA",
  high_detected_genes = "nFeature_RNA",
  high_total_counts = "nCount_RNA",
  high_mitochondrial_percentage = "pct_mt",
  high_ribosomal_percentage = "pct_ribo",
  high_hemoglobin_percentage = "pct_hb"
)
quantile_lookup <- list(
  low_detected_genes = q_gene, high_detected_genes = q_gene, high_total_counts = q_count,
  high_mitochondrial_percentage = q_mt, high_ribosomal_percentage = q_ribo,
  high_hemoglobin_percentage = q_hb
)
flag_lookup <- c(
  low_detected_genes = "qc_low_detected_genes",
  high_detected_genes = "qc_high_detected_genes",
  high_total_counts = "qc_high_total_counts",
  high_mitochondrial_percentage = "qc_high_mitochondrial_percentage",
  high_ribosomal_percentage = "qc_high_ribosomal_percentage",
  high_hemoglobin_percentage = "qc_high_hemoglobin_percentage"
)
qc_summary <- rbindlist(lapply(names(thresholds), function(name) {
  q <- quantile_lookup[[name]]
  flag <- object[[flag_lookup[[name]], drop = TRUE]]
  data.frame(
    qc_flag = name,
    metric = metric_name[[name]],
    threshold = thresholds[[name]],
    threshold_rule = threshold_rule[[name]],
    q1 = q[["q1"]], q3 = q[["q3"]], iqr = q[["iqr"]],
    affected_cells = sum(flag), total_cells = length(flag), affected_fraction = mean(flag),
    stringsAsFactors = FALSE
  )
}))
qc_summary <- rbind(qc_summary, data.frame(
  qc_flag = "qc_pass", metric = "all_flags", threshold = NA_real_,
  threshold_rule = "TRUE when none of the six QC flags is TRUE",
  q1 = NA_real_, q3 = NA_real_, iqr = NA_real_,
  affected_cells = sum(object$qc_pass), total_cells = ncol(object), affected_fraction = mean(object$qc_pass)
))
write_tsv(qc_summary, file.path(out_dir, "gse189926_r_qc_summary.tsv"))

qc_by_sample <- rbindlist(lapply(split(object[[]], object$sample_accession), function(d) {
  data.frame(
    sample_accession = d$sample_accession[[1L]],
    patient_id = d$patient_id[[1L]],
    timepoint = d$timepoint[[1L]],
    raw_outcome = d[["characteristics_ch1::outcome"]][[1L]],
    cells_input = nrow(d),
    cells_qc_pass = sum(d$qc_pass),
    qc_pass_fraction = mean(d$qc_pass),
    median_nCount_RNA = median(d$nCount_RNA),
    median_nFeature_RNA = median(d$nFeature_RNA),
    median_pct_mt = median(d$pct_mt),
    median_pct_ribo = median(d$pct_ribo),
    median_pct_hb = median(d$pct_hb),
    low_detected_genes = sum(d$qc_low_detected_genes),
    high_detected_genes = sum(d$qc_high_detected_genes),
    high_total_counts = sum(d$qc_high_total_counts),
    high_mitochondrial_percentage = sum(d$qc_high_mitochondrial_percentage),
    high_ribosomal_percentage = sum(d$qc_high_ribosomal_percentage),
    high_hemoglobin_percentage = sum(d$qc_high_hemoglobin_percentage),
    stringsAsFactors = FALSE
  )
}))
write_tsv(qc_by_sample, file.path(out_dir, "gse189926_r_qc_flag_by_sample.tsv"))

raw_object_path <- file.path(object_dir, "gse189926_r_unfiltered_raw.rds")
log_step(sprintf("Saving unfiltered raw object: %s", raw_object_path))
saveRDS(object, raw_object_path, compress = FALSE)

log_step("Running unfiltered unintegrated Seurat route")
unfiltered <- run_log_normalized_route(object, reduction_prefix = "unfiltered", run_clusters = FALSE)
unfiltered_audit <- attr(unfiltered, "route_audit")
unfiltered_coords <- Embeddings(unfiltered, reduction = "umap.unfiltered")
unfiltered_metadata <- unfiltered[[]][
  rownames(unfiltered_coords),
  c(required_metadata_fields, "timepoint_simple", "binary_response_group", "qc_pass"),
  drop = FALSE
]
unfiltered_plot_data <- cbind(
  cell_id = rownames(unfiltered_coords),
  UMAP_1 = unfiltered_coords[, 1L], UMAP_2 = unfiltered_coords[, 2L],
  unfiltered_metadata
)
saveRDS(unfiltered_plot_data, file.path(object_dir, "gse189926_r_unfiltered_umap_plot_data.rds"), compress = FALSE)
unfiltered_mix <- mixing_metrics(unfiltered_coords, unfiltered[[]], "unfiltered_unintegrated_umap", k = 15L)
rm(unfiltered, unfiltered_coords, unfiltered_plot_data, unfiltered_metadata)
gc(verbose = FALSE)

log_step("Reloading raw object and creating QC-pass subset")
object <- readRDS(raw_object_path)
qc_object <- subset(object, subset = qc_pass)
rm(object)
gc(verbose = FALSE)
log_step(sprintf("QC-pass subset has %d cells", ncol(qc_object)))

log_step("Running QC-pass unintegrated Seurat route")
qc_object <- run_log_normalized_route(qc_object, reduction_prefix = "qc", run_clusters = FALSE)
qc_route_audit <- attr(qc_object, "route_audit")
qc_mix <- mixing_metrics(Embeddings(qc_object, "umap.qc"), qc_object[[]], "qc_pass_unintegrated_umap", k = 15L)

log_step("Running Harmony correction by sample_accession")
qc_object <- harmony::RunHarmony(
  qc_object,
  group.by.vars = "sample_accession",
  reduction.use = "pca.qc",
  dims.use = 1:30,
  reduction.save = "harmony",
  plot_convergence = FALSE,
  verbose = FALSE
)
qc_object <- FindNeighbors(
  qc_object,
  reduction = "harmony",
  dims = 1:30,
  k.param = 15,
  graph.name = c("harmony_nn", "harmony_snn"),
  verbose = FALSE
)
for (res in c(0.3, 0.5, 0.8)) {
  qc_object <- FindClusters(
    qc_object,
    graph.name = "harmony_snn",
    resolution = res,
    algorithm = 1,
    random.seed = 340,
    verbose = FALSE
  )
  field <- switch(as.character(res), "0.3" = "cluster_harmony_r0_3", "0.5" = "cluster_harmony_r0_5", "0.8" = "cluster_harmony_r0_8")
  qc_object[[field]] <- as.character(Idents(qc_object))
}
Idents(qc_object) <- "cluster_harmony_r0_5"
qc_object <- RunUMAP(
  qc_object,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 15,
  min.dist = 0.3,
  seed.use = 340,
  reduction.name = "umap.harmony",
  verbose = FALSE
)
harmony_mix <- mixing_metrics(Embeddings(qc_object, "umap.harmony"), qc_object[[]], "qc_pass_harmony_sample_accession_umap", k = 15L)
write_tsv(rbind(unfiltered_mix, qc_mix, harmony_mix), file.path(out_dir, "gse189926_r_mixing_metrics.tsv"))

integration_audit <- data.frame(
  route_name = c("unfiltered_unintegrated", "qc_pass_unintegrated", "qc_pass_harmony_sample_accession"),
  language = "R",
  package_route = c("Seurat", "Seurat", "Seurat+harmony"),
  cells_used = c(nrow(unfiltered_mix) * 0 + unfiltered_mix$cells_used, qc_mix$cells_used, harmony_mix$cells_used),
  normalization = "NormalizeData LogNormalize scale.factor=10000",
  hvg_rule = c(
    sprintf("FindVariableFeatures vst nfeatures=3000; initial=%d", unfiltered_audit$hvg_initial),
    sprintf("FindVariableFeatures vst nfeatures=3000; initial=%d", qc_route_audit$hvg_initial),
    sprintf("reused QC-pass PCA; initial HVG=%d", qc_route_audit$hvg_initial)
  ),
  technical_gene_exclusion_rule = "Excluded from PCA-driving HVG only: ^MT-, ^RPS/^RPL, hemoglobin regex ^HB[ABDEGMQZ][0-9]*$, MALAT1/MALAT1-, and immunoglobulin constant-region prefixes IGHA/IGHG/IGHM/IGHD/IGHE/IGKC/IGLC. Genes remain in counts, data, and raw marker table.",
  technical_hvg_excluded = c(unfiltered_audit$technical_hvg_excluded, qc_route_audit$technical_hvg_excluded, qc_route_audit$technical_hvg_excluded),
  pca_features = c(unfiltered_audit$pca_features, qc_route_audit$pca_features, qc_route_audit$pca_features),
  pca_components_computed = 40,
  dimensions_used = "1:30",
  neighbor_parameters = "FindNeighbors k.param=15",
  umap_parameters = "RunUMAP n.neighbors=15 min.dist=0.3 seed.use=340",
  correction_field = c("", "", "sample_accession"),
  correction_method = c("none", "none", "Harmony"),
  clustering_resolutions = c("not run", "not run", "0.3;0.5;0.8"),
  notes = c(
    "All cells retained; labels for QC-failed cells will remain QC_failed_not_labeled in the final figure.",
    "QC-pass cells; no sample-aware correction.",
    "QC-pass cells; patient_id was not used as a correction field."
  ),
  stringsAsFactors = FALSE
)
write_tsv(integration_audit, file.path(out_dir, "gse189926_r_integration_strategy_audit.tsv"))

resolution_rows <- rbindlist(lapply(c("cluster_harmony_r0_3", "cluster_harmony_r0_5", "cluster_harmony_r0_8"), function(field) {
  counts_by_cluster <- table(qc_object[[field, drop = TRUE]])
  data.frame(
    cluster_field = field,
    resolution = switch(field, cluster_harmony_r0_3 = 0.3, cluster_harmony_r0_5 = 0.5, cluster_harmony_r0_8 = 0.8),
    cells_used = sum(counts_by_cluster),
    number_of_clusters = length(counts_by_cluster),
    minimum_cluster_cells = min(counts_by_cluster),
    median_cluster_cells = median(as.numeric(counts_by_cluster)),
    maximum_cluster_cells = max(counts_by_cluster),
    clusters_below_50_cells = sum(counts_by_cluster < 50),
    cluster_cell_counts = paste(names(counts_by_cluster), as.integer(counts_by_cluster), sep = ":", collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
write_tsv(resolution_rows, file.path(out_dir, "gse189926_r_cluster_resolution_summary.tsv"))

log_step("Computing full-cell sparse marker summaries for Harmony resolution 0.5")
markers <- marker_tables(qc_object, "cluster_harmony_r0_5")
write_tsv(markers$raw, file.path(out_dir, "gse189926_r_marker_top50_raw.tsv"))
write_tsv(markers$no_technical, file.path(out_dir, "gse189926_r_marker_top50_no_technical.tsv"))
rm(markers)

prelabel_path <- file.path(object_dir, "gse189926_r_qc_pass_prelabel.rds")
log_step(sprintf("Saving QC-pass prelabel object: %s", prelabel_path))
qc_object <- DietSeurat(
  qc_object,
  assays = "RNA",
  layers = c("counts", "data"),
  dimreducs = c("pca.qc", "harmony", "umap.qc", "umap.harmony"),
  graphs = NULL,
  misc = TRUE
)
saveRDS(qc_object, prelabel_path, compress = FALSE)

stage1_status <- data.frame(
  dataset_id = "GSE189926",
  step = "R_sparse_import_qc_seurat_harmony_clustering_marker_stage1",
  status = "completed",
  language = "R",
  package_route = "Seurat 5.4.0; harmony 1.2.4; Matrix 1.7.4; data.table 1.18.2.1; RANN 2.6.2",
  input_files = "22 raw gzip text matrices; gse189926_sample_qc.tsv; gse189926_sample_parse_audit.tsv",
  output_files = "stage-1 audits, mixing metrics, cluster resolution summary, raw and no-technical marker tables, local RDS objects",
  cell_count_input = nrow(unfiltered_mix) * 0 + unfiltered_mix$cells_used,
  cell_count_used = ncol(qc_object),
  notes = "Awaiting evidence-based cluster label review and final figure generation in R stage 2.",
  error_message = "",
  stringsAsFactors = FALSE
)
write_tsv(stage1_status, file.path(out_dir, "stage1_run_status.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "stage1_session_info.txt"))
log_step("Stage 1 completed")
