options(stringsAsFactors = FALSE)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (tolower(Sys.getenv("CODEX_AUTO_INSTALL_R_PACKAGES", "true")) %in% c("1", "true", "yes")) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package missing: ", pkg)
  }
}

for (pkg in c("Matrix", "Seurat", "data.table", "ggplot2", "patchwork")) {
  install_if_missing(pkg)
}

library(Matrix)
library(Seurat)
library(data.table)
library(ggplot2)
library(patchwork)

find_repo_root <- function(start) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, ".git"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Cannot find repository root")
    }
    current <- parent
  }
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  script_path <- normalizePath("bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_stromal_subset_clustering/scripts/pc_stromal_subset_clustering.R", winslash = "/", mustWork = TRUE)
}

request_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- find_repo_root(request_dir)
project_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche")
working_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_working_qc_sparse_objects_20260708")
review_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_stromal_first_pass_review_20260708")
cluster_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_stromal_subset_clustering_20260708")
upload_dir <- file.path(project_dir, "desktop_exchange", "uploads", "20260708_stromal_subset_clustering")
dir.create(cluster_object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  data.table::fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

file_size_or_na <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  as.numeric(file.info(path)$size)
}

safe_feature_symbols <- function(feature_df) {
  if (ncol(feature_df) >= 2) {
    return(as.character(feature_df[[2]]))
  }
  as.character(feature_df[[1]])
}

aggregate_sparse_by_symbol <- function(counts, symbols) {
  if (length(symbols) != nrow(counts)) {
    stop("Feature symbol length does not match count matrix rows")
  }
  symbols <- as.character(symbols)
  keep <- !is.na(symbols) & symbols != ""
  counts <- counts[keep, , drop = FALSE]
  symbols <- symbols[keep]
  group_levels <- unique(symbols)
  group_index <- match(symbols, group_levels)
  mapper <- Matrix::sparseMatrix(
    i = group_index,
    j = seq_along(group_index),
    x = 1,
    dims = c(length(group_levels), length(group_index))
  )
  aggregated <- mapper %*% as(counts, "dgCMatrix")
  rownames(aggregated) <- group_levels
  colnames(aggregated) <- colnames(counts)
  as(aggregated, "dgCMatrix")
}

align_to_features <- function(mat, all_features) {
  mat <- as(mat, "dgCMatrix")
  entries <- summary(mat)
  if (nrow(entries) == 0) {
    aligned <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(length(all_features), ncol(mat)))
  } else {
    feature_index <- match(rownames(mat), all_features)
    aligned <- Matrix::sparseMatrix(
      i = feature_index[entries$i],
      j = entries$j,
      x = entries$x,
      dims = c(length(all_features), ncol(mat))
    )
  }
  rownames(aligned) <- all_features
  colnames(aligned) <- colnames(mat)
  as(aligned, "dgCMatrix")
}

standardize_meta_value <- function(x, n, field_name) {
  if (is.null(x)) {
    return(rep("", n))
  }
  if (length(x) == 1) {
    return(rep(as.character(x), n))
  }
  if (length(x) == n) {
    return(as.character(x))
  }
  stop("Unexpected field length for ", field_name, ": ", length(x), "; expected 1 or ", n)
}

top_n_by_cluster <- function(markers, n) {
  markers <- as.data.table(markers)
  if (nrow(markers) == 0) {
    return(markers)
  }
  score_col <- if ("avg_log2FC" %in% names(markers)) "avg_log2FC" else if ("avg_logFC" %in% names(markers)) "avg_logFC" else NULL
  if (is.null(score_col)) {
    markers[, marker_rank := seq_len(.N), by = cluster]
    return(markers[marker_rank <= n])
  }
  setorderv(markers, c("cluster", score_col, "pct.1"), c(1, -1, -1))
  markers[, marker_rank := seq_len(.N), by = cluster]
  markers[marker_rank <= n]
}

make_count_table <- function(meta, group_cols, value_col = "main_cluster") {
  dt <- as.data.table(meta)
  cols <- c("dataset_id", value_col, group_cols)
  dt[, .N, by = cols][order(dataset_id, get(value_col))]
}

save_umap <- function(obj, group_by, path, title_text, width = 8.5, height = 6.5) {
  p <- DimPlot(obj, reduction = "umap", group.by = group_by, pt.size = 0.18, raster = FALSE) +
    ggtitle(title_text) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    )
  ggsave(path, p, width = width, height = height, dpi = 220, bg = "white")
}

save_marker_dotplot <- function(obj, path, title_text, width = 12, height = 5.8) {
  marker_genes <- c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PRELP",
    "FAP", "ACTA2", "POSTN", "THBS2", "MMP2", "MMP11",
    "CXCL14", "CFD", "C3", "IL6", "CXCL12", "SFRP1",
    "MCAM", "RGS5", "PDGFRB", "CSPG4",
    "PECAM1", "VWF", "KDR",
    "EPCAM", "KRT8", "KRT18", "KRT19",
    "PTPRC", "CD3D", "LYZ",
    "MKI67", "TOP2A"
  )
  marker_genes <- intersect(marker_genes, rownames(obj))
  if (length(marker_genes) == 0) {
    return(invisible(FALSE))
  }
  p <- DotPlot(obj, features = marker_genes, group.by = "main_cluster") +
    ggtitle(title_text) +
    coord_flip() +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    )
  ggsave(path, p, width = width, height = height, dpi = 220, bg = "white")
  invisible(TRUE)
}

review_classes_to_keep <- c("caf_ecm_review", "pericyte_review", "endothelial_review")
dataset_ids <- c("GSE178318", "GSE245552")

load_dataset_subset <- function(dataset_id) {
  working_paths <- sort(list.files(working_object_dir, pattern = "_working_qc_sparse_matrix[.]rds$", full.names = TRUE))
  dataset_mats <- list()
  dataset_meta <- list()
  subset_summary <- list()

  for (working_path in working_paths) {
    object_label <- sub("[.]rds$", "", basename(working_path))
    review_path <- file.path(review_object_dir, paste0(object_label, "_stromal_first_pass_review.rds"))
    if (!file.exists(review_path)) {
      next
    }
    working_obj <- readRDS(working_path)
    review_obj <- readRDS(review_path)
    if (!identical(as.character(working_obj$dataset_id), dataset_id)) {
      next
    }
    required_review_fields <- c("dataset_id", "object_label", "review_table")
    missing_review_fields <- setdiff(required_review_fields, names(review_obj))
    if (length(missing_review_fields) > 0) {
      stop("Review object missing fields: ", paste(missing_review_fields, collapse = ", "), " in ", review_path)
    }
    if (!all(c("counts", "features", "metadata") %in% names(working_obj))) {
      stop("Working object missing counts, features, or metadata: ", working_path)
    }
    review_table <- as.data.frame(review_obj$review_table)
    required_review_table_fields <- c("review_class", "top_marker_program", "score_margin")
    missing_review_table_fields <- setdiff(required_review_table_fields, names(review_table))
    if (length(missing_review_table_fields) > 0) {
      stop("Review table missing fields: ", paste(missing_review_table_fields, collapse = ", "), " in ", review_path)
    }
    counts <- working_obj$counts
    metadata <- as.data.frame(working_obj$metadata)
    if (ncol(counts) != nrow(review_table)) {
      stop("Cell mismatch in ", object_label, ": count columns ", ncol(counts), "; review rows ", nrow(review_table))
    }
    keep <- review_table$review_class %in% review_classes_to_keep
    subset_summary[[length(subset_summary) + 1]] <- data.frame(
      dataset_id = dataset_id,
      object_label = object_label,
      input_cells = ncol(counts),
      retained_cells = sum(keep),
      retained_fraction = sum(keep) / ncol(counts)
    )
    if (!any(keep)) {
      next
    }
    symbols <- safe_feature_symbols(working_obj$features)
    sub_counts <- counts[, keep, drop = FALSE]
    sub_counts <- aggregate_sparse_by_symbol(sub_counts, symbols)
    raw_meta <- metadata[keep, , drop = FALSE]
    sub_review <- review_table[keep, , drop = FALSE]
    cell_index <- which(keep)
    raw_barcode <- if ("barcode" %in% names(raw_meta)) as.character(raw_meta$barcode) else as.character(cell_index)
    cell_id <- paste(dataset_id, object_label, cell_index, sep = "__")
    colnames(sub_counts) <- cell_id
    sample_label <- if ("sample_label" %in% names(sub_review)) standardize_meta_value(sub_review$sample_label, sum(keep), "sample_label") else
      if ("sample_label" %in% names(raw_meta)) standardize_meta_value(raw_meta$sample_label, sum(keep), "sample_label") else
        if (!is.null(working_obj$sample_label)) rep(as.character(working_obj$sample_label), sum(keep)) else rep("", sum(keep))
    tissue_site <- if ("tissue_site" %in% names(sub_review)) standardize_meta_value(sub_review$tissue_site, sum(keep), "tissue_site") else
      if ("tissue_site" %in% names(raw_meta)) standardize_meta_value(raw_meta$tissue_site, sum(keep), "tissue_site") else
        if (!is.null(working_obj$tissue_site)) rep(as.character(working_obj$tissue_site), sum(keep)) else rep("", sum(keep))
    meta_out <- data.frame(
      cell_id = cell_id,
      dataset_id = dataset_id,
      object_label = object_label,
      sample_label = sample_label,
      tissue_site = tissue_site,
      review_class = as.character(sub_review$review_class),
      top_marker_program = as.character(sub_review$top_marker_program),
      score_margin = as.numeric(sub_review$score_margin),
      raw_barcode = raw_barcode,
      stringsAsFactors = FALSE
    )
    rownames(meta_out) <- cell_id
    dataset_mats[[length(dataset_mats) + 1]] <- sub_counts
    dataset_meta[[length(dataset_meta) + 1]] <- meta_out
    rm(working_obj, review_obj, counts, sub_counts)
    gc()
  }

  if (length(dataset_mats) == 0) {
    stop("No retained stromal-supported cells for ", dataset_id)
  }
  all_features <- sort(unique(unlist(lapply(dataset_mats, rownames), use.names = FALSE)))
  aligned <- lapply(dataset_mats, align_to_features, all_features = all_features)
  combined_counts <- do.call(cbind, aligned)
  combined_meta <- do.call(rbind, dataset_meta)
  if (!identical(colnames(combined_counts), rownames(combined_meta))) {
    stop("Combined count matrix and metadata rownames are not aligned for ", dataset_id)
  }
  list(
    counts = combined_counts,
    metadata = combined_meta,
    subset_summary = do.call(rbind, subset_summary)
  )
}

run_dataset <- function(dataset_id) {
  loaded <- load_dataset_subset(dataset_id)
  obj <- CreateSeuratObject(
    counts = loaded$counts,
    meta.data = loaded$metadata,
    project = paste0(dataset_id, "_stromal_subset"),
    min.cells = 0,
    min.features = 0
  )
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
  dims_use <- 1:min(20, ncol(Embeddings(obj, "pca")))
  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  for (resolution_value in c(0.2, 0.4, 0.6)) {
    obj <- FindClusters(obj, resolution = resolution_value, verbose = FALSE)
  }
  cluster_cols <- grep("_snn_res[.]0[.]4$", colnames(obj@meta.data), value = TRUE)
  if (length(cluster_cols) == 0) {
    cluster_cols <- grep("res[.]0[.]4$", colnames(obj@meta.data), value = TRUE)
  }
  if (length(cluster_cols) == 0) {
    cluster_col <- "seurat_clusters"
  } else {
    cluster_col <- cluster_cols[1]
  }
  obj$main_cluster <- as.character(obj@meta.data[[cluster_col]])
  Idents(obj) <- "main_cluster"
  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)

  full_marker_path <- file.path(cluster_object_dir, paste0(tolower(dataset_id), "_stromal_subset_all_markers_res0.4.rds"))
  markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  )
  saveRDS(markers, full_marker_path)

  object_path <- file.path(cluster_object_dir, paste0(tolower(dataset_id), "_stromal_subset_seurat_res0.4.rds"))
  saveRDS(obj, object_path)

  meta <- obj@meta.data
  meta$main_cluster <- as.character(meta$main_cluster)
  cluster_counts <- as.data.table(meta)[, .N, by = .(dataset_id, main_cluster)]
  cluster_counts[, fraction := N / sum(N), by = dataset_id]
  composition_review <- as.data.table(meta)[, .N, by = .(dataset_id, main_cluster, review_class)]
  composition_sample <- as.data.table(meta)[, .N, by = .(dataset_id, main_cluster, sample_label)]
  composition_tissue <- as.data.table(meta)[, .N, by = .(dataset_id, main_cluster, tissue_site)]
  composition_program <- as.data.table(meta)[, .N, by = .(dataset_id, main_cluster, top_marker_program)]

  marker_genes <- c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PRELP",
    "FAP", "ACTA2", "POSTN", "THBS2", "MMP2", "MMP11",
    "CXCL14", "CFD", "C3", "IL6", "CXCL12", "SFRP1",
    "MCAM", "RGS5", "PDGFRB", "CSPG4",
    "PECAM1", "VWF", "KDR",
    "EPCAM", "KRT8", "KRT18", "KRT19",
    "PTPRC", "CD3D", "LYZ",
    "MKI67", "TOP2A"
  )
  marker_genes <- intersect(marker_genes, rownames(obj))
  avg_expr <- AverageExpression(obj, features = marker_genes, group.by = "main_cluster", assays = DefaultAssay(obj), slot = "data", verbose = FALSE)
  avg_mat <- avg_expr[[DefaultAssay(obj)]]
  marker_score <- as.data.table(as.data.frame(as.table(as.matrix(avg_mat))))
  names(marker_score) <- c("gene_symbol", "main_cluster", "mean_log_norm")
  marker_score[, dataset_id := dataset_id]
  setcolorder(marker_score, c("dataset_id", "main_cluster", "gene_symbol", "mean_log_norm"))

  dataset_prefix <- tolower(dataset_id)
  save_umap(
    obj,
    "main_cluster",
    file.path(upload_dir, paste0(dataset_prefix, "_stromal_subset_umap_by_cluster.png")),
    paste0(dataset_id, " stromal subset by cluster")
  )
  save_umap(
    obj,
    "review_class",
    file.path(upload_dir, paste0(dataset_prefix, "_stromal_subset_umap_by_review_class.png")),
    paste0(dataset_id, " stromal subset by review class")
  )
  if (dataset_id == "GSE245552") {
    save_umap(
      obj,
      "tissue_site",
      file.path(upload_dir, paste0(dataset_prefix, "_stromal_subset_umap_by_tissue.png")),
      paste0(dataset_id, " stromal subset by tissue"),
      width = 9.5,
      height = 6.5
    )
  }
  save_marker_dotplot(
    obj,
    file.path(upload_dir, paste0(dataset_prefix, "_stromal_subset_marker_dotplot.png")),
    paste0(dataset_id, " stromal subset marker dot plot")
  )

  data.frame_summary <- data.frame(
    dataset_id = dataset_id,
    total_working_cells_reviewed = sum(loaded$subset_summary$input_cells),
    stromal_subset_cells = ncol(obj),
    stromal_subset_fraction = ncol(obj) / sum(loaded$subset_summary$input_cells),
    features = nrow(obj),
    variable_features = length(VariableFeatures(obj)),
    pca_dims_used = length(dims_use),
    main_cluster_column = cluster_col,
    main_cluster_count = length(unique(obj$main_cluster)),
    stringsAsFactors = FALSE
  )
  manifest <- data.frame(
    dataset_id = dataset_id,
    seurat_object_path = object_path,
    seurat_object_size_bytes = file_size_or_na(object_path),
    all_marker_rds_path = full_marker_path,
    all_marker_rds_size_bytes = file_size_or_na(full_marker_path),
    cells = ncol(obj),
    features = nrow(obj),
    stringsAsFactors = FALSE
  )
  list(
    summary = data.frame_summary,
    subset_summary = loaded$subset_summary,
    manifest = manifest,
    cluster_counts = cluster_counts,
    composition_review = composition_review,
    composition_sample = composition_sample,
    composition_tissue = composition_tissue,
    composition_program = composition_program,
    markers_top50 = top_n_by_cluster(markers, 50),
    markers_top10 = top_n_by_cluster(markers, 10),
    marker_score = marker_score
  )
}

results <- lapply(dataset_ids, run_dataset)

dataset_summary <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)
subset_summary <- rbindlist(lapply(results, `[[`, "subset_summary"), fill = TRUE)
run_manifest <- rbindlist(lapply(results, `[[`, "manifest"), fill = TRUE)
cluster_counts <- rbindlist(lapply(results, `[[`, "cluster_counts"), fill = TRUE)
composition_review <- rbindlist(lapply(results, `[[`, "composition_review"), fill = TRUE)
composition_sample <- rbindlist(lapply(results, `[[`, "composition_sample"), fill = TRUE)
composition_tissue <- rbindlist(lapply(results, `[[`, "composition_tissue"), fill = TRUE)
composition_program <- rbindlist(lapply(results, `[[`, "composition_program"), fill = TRUE)
markers_top50 <- rbindlist(lapply(results, `[[`, "markers_top50"), fill = TRUE)
markers_top10 <- rbindlist(lapply(results, `[[`, "markers_top10"), fill = TRUE)
marker_score <- rbindlist(lapply(results, `[[`, "marker_score"), fill = TRUE)

write_tsv(run_manifest, file.path(upload_dir, "stromal_subset_run_manifest.tsv"))
write_tsv(dataset_summary, file.path(upload_dir, "stromal_subset_dataset_summary.tsv"))
write_tsv(cluster_counts, file.path(upload_dir, "stromal_subset_cluster_counts.tsv"))
write_tsv(composition_review, file.path(upload_dir, "stromal_subset_cluster_composition_by_review_class.tsv"))
write_tsv(composition_sample, file.path(upload_dir, "stromal_subset_cluster_composition_by_sample.tsv"))
write_tsv(composition_tissue, file.path(upload_dir, "stromal_subset_cluster_composition_by_tissue.tsv"))
write_tsv(composition_program, file.path(upload_dir, "stromal_subset_cluster_composition_by_top_marker_program.tsv"))
write_tsv(markers_top50, file.path(upload_dir, "stromal_subset_cluster_markers_top50.tsv"))
write_tsv(markers_top10, file.path(upload_dir, "stromal_subset_cluster_markers_top10.tsv"))
write_tsv(marker_score, file.path(upload_dir, "stromal_subset_marker_score_by_cluster.tsv"))

status_text <- c(
  "# CRLM CAF/ECM stromal subset clustering",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  paste0("Working object directory: `", working_object_dir, "`"),
  paste0("Review object directory: `", review_object_dir, "`"),
  paste0("PC-local clustering object directory: `", cluster_object_dir, "`"),
  "",
  paste0("Datasets analyzed: ", paste(dataset_ids, collapse = ", ")),
  paste0("Total stromal-supported cells analyzed: ", sum(dataset_summary$stromal_subset_cells)),
  paste0("Dataset summary rows: ", nrow(dataset_summary)),
  paste0("Cluster count rows: ", nrow(cluster_counts)),
  paste0("Review-class composition rows: ", nrow(composition_review)),
  paste0("Sample composition rows: ", nrow(composition_sample)),
  paste0("Tissue composition rows: ", nrow(composition_tissue)),
  paste0("Top50 marker rows: ", nrow(markers_top50)),
  paste0("Top10 marker rows: ", nrow(markers_top10)),
  paste0("Marker score rows: ", nrow(marker_score)),
  "",
  "This request produced review-grade clustering and all-gene marker evidence.",
  "It does not finalize cell labels.",
  "PC-local Seurat objects and full marker RDS files should not be committed to GitHub."
)
writeLines(status_text, file.path(upload_dir, "STATUS.md"))
