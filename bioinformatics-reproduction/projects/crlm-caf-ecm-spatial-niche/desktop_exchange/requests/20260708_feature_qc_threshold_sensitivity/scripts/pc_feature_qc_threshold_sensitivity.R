options(stringsAsFactors = FALSE)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Required R package missing: Matrix")
}

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
  script_path <- normalizePath("bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_feature_qc_threshold_sensitivity/scripts/pc_feature_qc_threshold_sensitivity.R", winslash = "/", mustWork = TRUE)
}

request_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- find_repo_root(request_dir)
raw_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche", "desktop_exchange", "requests", "20260708_core_geo_download", "raw_core_geo")
upload_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche", "desktop_exchange", "uploads", "20260708_feature_qc_threshold_sensitivity")
input_metadata <- file.path(request_dir, "inputs", "sample_metadata_from_geo_series.tsv")
dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

metric_summary <- function(x) {
  x <- as.numeric(x)
  data.frame(
    n = length(x),
    min = suppressWarnings(min(x, na.rm = TRUE)),
    q1 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)),
    median = stats::median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    q3 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)),
    max = suppressWarnings(max(x, na.rm = TRUE))
  )
}

read_gz_table_base <- function(path, header = FALSE) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, header = header, check.names = FALSE, sep = "\t")
}

read_lines_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con)
}

mt_percent_from_values <- function(mat, values) {
  mt <- grepl("^MT-", values)
  total <- Matrix::colSums(mat)
  mt_counts <- if (any(mt)) Matrix::colSums(mat[mt, , drop = FALSE]) else rep(0, ncol(mat))
  ifelse(total > 0, mt_counts / total * 100, NA_real_)
}

feature_column_audit <- function(features, dataset_id, gsm, sample_label, tissue_site, feature_member) {
  out <- list()
  for (i in seq_len(ncol(features))) {
    values <- as.character(features[[i]])
    out[[length(out) + 1]] <- data.frame(
      dataset_id = dataset_id,
      gsm = gsm,
      sample_label = sample_label,
      tissue_site = tissue_site,
      feature_file = feature_member,
      feature_rows = nrow(features),
      feature_columns = ncol(features),
      feature_column_index = i,
      non_missing = sum(!is.na(values) & values != ""),
      unique_values = length(unique(values)),
      mt_upper_prefix_count = sum(grepl("^MT-", values)),
      mt_lower_prefix_count = sum(grepl("^mt-", values)),
      mt_title_prefix_count = sum(grepl("^Mt-", values)),
      ensg_prefix_count = sum(grepl("^ENSG", values)),
      example_values = paste(utils::head(values, 8), collapse = " | ")
    )
  }
  do.call(rbind, out)
}

summarise_mito_by_feature_column <- function(mat, features, dataset_id, gsm, sample_label, tissue_site) {
  out <- list()
  for (i in seq_len(ncol(features))) {
    values <- as.character(features[[i]])
    percent_mt <- mt_percent_from_values(mat, values)
    s <- metric_summary(percent_mt)
    out[[length(out) + 1]] <- data.frame(
      dataset_id = dataset_id,
      gsm = gsm,
      sample_label = sample_label,
      tissue_site = tissue_site,
      feature_column_index = i,
      mt_upper_prefix_count = sum(grepl("^MT-", values)),
      s
    )
  }
  do.call(rbind, out)
}

retention_by_group <- function(qc, group_col, dataset_id, min_features_grid, max_percent_mt_grid = NULL, min_counts_grid = NULL) {
  if (is.null(min_counts_grid)) {
    min_counts_grid <- 0
  }
  if (is.null(max_percent_mt_grid)) {
    max_percent_mt_grid <- NA_real_
  }
  out <- list()
  groups <- sort(unique(as.character(qc[[group_col]])))
  for (group_value in groups) {
    idx <- which(as.character(qc[[group_col]]) == group_value)
    n_start <- length(idx)
    for (min_features in min_features_grid) {
      for (min_counts in min_counts_grid) {
        for (max_percent_mt in max_percent_mt_grid) {
          keep <- qc$nFeature_RNA[idx] >= min_features & qc$nCount_RNA[idx] >= min_counts
          if (!is.na(max_percent_mt)) {
            keep <- keep & qc$percent.mt[idx] <= max_percent_mt
          }
          retained <- sum(keep, na.rm = TRUE)
          out[[length(out) + 1]] <- data.frame(
            dataset_id = dataset_id,
            group_field = group_col,
            group_value = group_value,
            start_cells = n_start,
            min_features = min_features,
            min_counts = min_counts,
            max_percent_mt = ifelse(is.na(max_percent_mt), "", as.character(max_percent_mt)),
            retained_cells = retained,
            retained_fraction = if (n_start > 0) retained / n_start else NA_real_
          )
        }
      }
    }
  }
  do.call(rbind, out)
}

sample_meta <- read.delim(input_metadata, sep = "\t", header = TRUE, check.names = FALSE)
status_lines <- c()

gse178318_dir <- file.path(raw_dir, "GSE178318")
gse178318_matrix <- Matrix::readMM(gzfile(file.path(gse178318_dir, "GSE178318_matrix.mtx.gz"), open = "rt"))
gse178318_features <- read_gz_table_base(file.path(gse178318_dir, "GSE178318_genes.tsv.gz"), header = FALSE)
gse178318_barcodes <- read_lines_gz(file.path(gse178318_dir, "GSE178318_barcodes.tsv.gz"))
gse178318_symbols <- if (ncol(gse178318_features) >= 2) as.character(gse178318_features[[2]]) else as.character(gse178318_features[[1]])
gse178318_qc <- data.frame(
  barcode_suffix = sub("^[^_]+_", "", gse178318_barcodes),
  nCount_RNA = as.numeric(Matrix::colSums(gse178318_matrix)),
  nFeature_RNA = as.numeric(Matrix::colSums(gse178318_matrix > 0)),
  percent.mt = mt_percent_from_values(gse178318_matrix, gse178318_symbols)
)
gse178318_retention <- retention_by_group(
  gse178318_qc,
  "barcode_suffix",
  "GSE178318",
  min_features_grid = c(200, 300, 500),
  max_percent_mt_grid = c(10, 15, 20)
)
write_tsv(gse178318_retention, file.path(upload_dir, "gse178318_threshold_retention_by_suffix.tsv"))
status_lines <- c(status_lines, sprintf("GSE178318 threshold rows: %s", nrow(gse178318_retention)))
rm(gse178318_matrix, gse178318_qc)
gc()

gse245552_tar <- file.path(raw_dir, "GSE245552", "GSE245552_RAW.tar")
tar_members <- utils::untar(gse245552_tar, list = TRUE)
matrix_members <- sort(tar_members[grepl("_matrix[.]mtx[.]gz$", tar_members)])
gse245552_temp <- tempfile("gse245552_feature_qc_")
dir.create(gse245552_temp)

feature_audit_rows <- list()
mito_qc_rows <- list()
threshold_rows <- list()

for (matrix_member in matrix_members) {
  prefix <- sub("_matrix[.]mtx[.]gz$", "", matrix_member)
  barcode_member <- paste0(prefix, "_barcodes.tsv.gz")
  feature_member <- paste0(prefix, "_features.tsv.gz")
  utils::untar(gse245552_tar, files = c(matrix_member, barcode_member, feature_member), exdir = gse245552_temp)

  matrix_path <- file.path(gse245552_temp, matrix_member)
  barcode_path <- file.path(gse245552_temp, barcode_member)
  feature_path <- file.path(gse245552_temp, feature_member)
  mat <- Matrix::readMM(gzfile(matrix_path, open = "rt"))
  features <- read_gz_table_base(feature_path, header = FALSE)

  sample_label <- sub("^GSM[0-9]+_", "", prefix)
  gsm <- sub("_.*$", "", prefix)
  meta_row <- sample_meta[sample_meta$sample_id == gsm, , drop = FALSE]
  tissue_site <- if (nrow(meta_row) > 0) meta_row$tissue_site[1] else ""

  feature_audit_rows[[length(feature_audit_rows) + 1]] <- feature_column_audit(
    features,
    "GSE245552",
    gsm,
    sample_label,
    tissue_site,
    feature_member
  )
  mito_qc_rows[[length(mito_qc_rows) + 1]] <- summarise_mito_by_feature_column(
    mat,
    features,
    "GSE245552",
    gsm,
    sample_label,
    tissue_site
  )

  sample_qc <- data.frame(
    sample_label = rep(sample_label, ncol(mat)),
    nCount_RNA = as.numeric(Matrix::colSums(mat)),
    nFeature_RNA = as.numeric(Matrix::colSums(mat > 0)),
    percent.mt = rep(NA_real_, ncol(mat))
  )
  threshold_rows[[length(threshold_rows) + 1]] <- retention_by_group(
    sample_qc,
    "sample_label",
    "GSE245552",
    min_features_grid = c(200, 300, 500),
    min_counts_grid = c(0, 500, 1000),
    max_percent_mt_grid = NA_real_
  )

  rm(mat)
  unlink(c(matrix_path, barcode_path, feature_path))
  gc()
}

feature_audit <- do.call(rbind, feature_audit_rows)
mito_qc <- do.call(rbind, mito_qc_rows)
threshold_retention <- do.call(rbind, threshold_rows)

write_tsv(feature_audit, file.path(upload_dir, "gse245552_feature_column_audit.tsv"))
write_tsv(mito_qc, file.path(upload_dir, "gse245552_mito_qc_by_feature_column.tsv"))
write_tsv(threshold_retention, file.path(upload_dir, "gse245552_threshold_retention_by_sample.tsv"))
status_lines <- c(
  status_lines,
  sprintf("GSE245552 feature audit rows: %s", nrow(feature_audit)),
  sprintf("GSE245552 mitochondrial QC rows: %s", nrow(mito_qc)),
  sprintf("GSE245552 threshold rows: %s", nrow(threshold_retention))
)

status_text <- c(
  "# CRLM CAF/ECM feature QC and threshold sensitivity",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  paste0("Raw input directory: `", raw_dir, "`"),
  "",
  status_lines,
  "",
  "No filtering, normalization, clustering, annotation, or final inclusion decision was performed.",
  "GSE245552 mitochondrial fields are reported by exact feature column and are not used for final filtering in this request."
)
writeLines(status_text, file.path(upload_dir, "STATUS.md"))
