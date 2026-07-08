options(stringsAsFactors = FALSE)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Required R package missing: Matrix")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

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
  script_path <- normalizePath("bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_initial_qc_preflight/scripts/pc_initial_qc_preflight.R", winslash = "/", mustWork = TRUE)
}

request_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- find_repo_root(request_dir)
raw_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche", "desktop_exchange", "requests", "20260708_core_geo_download", "raw_core_geo")
upload_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche", "desktop_exchange", "uploads", "20260708_initial_qc_preflight")
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

summarise_by_group <- function(df, group_col, metric_cols, dataset_id) {
  out <- list()
  groups <- sort(unique(as.character(df[[group_col]])))
  for (group_value in groups) {
    idx <- which(as.character(df[[group_col]]) == group_value)
    for (metric in metric_cols) {
      if (!metric %in% names(df)) next
      s <- metric_summary(df[[metric]][idx])
      out[[length(out) + 1]] <- data.frame(
        dataset_id = dataset_id,
        group_field = group_col,
        group_value = group_value,
        metric = metric,
        s
      )
    }
  }
  if (length(out) == 0) {
    return(data.frame())
  }
  do.call(rbind, out)
}

read_gz_table_base <- function(path, header = TRUE) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, header = header, check.names = FALSE, sep = "\t")
}

read_lines_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con)
}

safe_gene_symbols <- function(feature_df) {
  if (ncol(feature_df) >= 2) {
    return(feature_df[[2]])
  }
  feature_df[[1]]
}

mt_fraction <- function(mat, gene_symbols) {
  mt <- grepl("^MT-", gene_symbols)
  total <- Matrix::colSums(mat)
  mt_counts <- if (any(mt)) Matrix::colSums(mat[mt, , drop = FALSE]) else rep(0, ncol(mat))
  ifelse(total > 0, mt_counts / total * 100, NA_real_)
}

summarise_sparse_matrix <- function(mat, gene_symbols, group_values, dataset_id, group_field) {
  qc <- data.frame(
    group = group_values,
    nCount_RNA = as.numeric(Matrix::colSums(mat)),
    nFeature_RNA = as.numeric(Matrix::colSums(mat > 0)),
    percent.mt = mt_fraction(mat, gene_symbols)
  )
  out <- list()
  for (group_value in sort(unique(as.character(qc$group)))) {
    idx <- which(as.character(qc$group) == group_value)
    for (metric in c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
      s <- metric_summary(qc[[metric]][idx])
      out[[length(out) + 1]] <- data.frame(
        dataset_id = dataset_id,
        group_field = group_field,
        group_value = group_value,
        metric = metric,
        s
      )
    }
  }
  do.call(rbind, out)
}

sample_meta <- read.delim(input_metadata, sep = "\t", header = TRUE, check.names = FALSE)

pc_status <- c()

gse225857_tar <- file.path(raw_dir, "GSE225857", "GSE225857_RAW.tar")
gse225857_temp <- tempfile("gse225857_")
dir.create(gse225857_temp)
untar(gse225857_tar, files = c("GSM7058755_non_immune_meta.txt.gz"), exdir = gse225857_temp)
nonimmune_meta <- read_gz_table_base(file.path(gse225857_temp, "GSM7058755_non_immune_meta.txt.gz"), header = TRUE)

columns_out <- data.frame(
  column_name = names(nonimmune_meta),
  class = vapply(nonimmune_meta, function(x) class(x)[1], character(1)),
  non_missing = vapply(nonimmune_meta, function(x) sum(!is.na(x) & as.character(x) != ""), integer(1))
)
write_tsv(columns_out, file.path(upload_dir, "gse225857_nonimmune_meta_columns.tsv"))

value_count_fields <- intersect(
  c("patients", "organs", "samples", "patients_organ", "cluster", "doublet", "predicted.doublet", "seurat_clusters", "batch"),
  names(nonimmune_meta)
)
value_rows <- list()
for (field in value_count_fields) {
  tab <- sort(table(as.character(nonimmune_meta[[field]])), decreasing = TRUE)
  value_rows[[length(value_rows) + 1]] <- data.frame(
    field = field,
    value = names(tab),
    cell_count = as.integer(tab)
  )
}
write_tsv(do.call(rbind, value_rows), file.path(upload_dir, "gse225857_nonimmune_meta_value_counts.tsv"))

qc_metric_cols <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "doublet.score"), names(nonimmune_meta))
qc_group_fields <- intersect(c("organs", "patients_organ", "samples", "batch", "doublet", "predicted.doublet"), names(nonimmune_meta))
qc_rows <- list()
for (field in qc_group_fields) {
  qc_rows[[length(qc_rows) + 1]] <- summarise_by_group(nonimmune_meta, field, qc_metric_cols, "GSE225857_nonimmune")
}
write_tsv(do.call(rbind, qc_rows), file.path(upload_dir, "gse225857_nonimmune_qc_by_field.tsv"))
pc_status <- c(pc_status, sprintf("GSE225857 nonimmune metadata cells: %s", nrow(nonimmune_meta)))

gse178318_dir <- file.path(raw_dir, "GSE178318")
gse178318_matrix <- Matrix::readMM(gzfile(file.path(gse178318_dir, "GSE178318_matrix.mtx.gz"), open = "rt"))
gse178318_features <- read_gz_table_base(file.path(gse178318_dir, "GSE178318_genes.tsv.gz"), header = FALSE)
gse178318_barcodes <- read_lines_gz(file.path(gse178318_dir, "GSE178318_barcodes.tsv.gz"))
barcode_suffix <- sub("^[^_]+_", "", gse178318_barcodes)
suffix_tab <- sort(table(barcode_suffix), decreasing = TRUE)
write_tsv(data.frame(barcode_suffix = names(suffix_tab), cell_count = as.integer(suffix_tab)), file.path(upload_dir, "gse178318_barcode_suffix_counts.tsv"))
gse178318_qc <- summarise_sparse_matrix(gse178318_matrix, safe_gene_symbols(gse178318_features), barcode_suffix, "GSE178318", "barcode_suffix")
write_tsv(gse178318_qc, file.path(upload_dir, "gse178318_qc_by_suffix.tsv"))
pc_status <- c(pc_status, sprintf("GSE178318 matrix cells: %s", ncol(gse178318_matrix)))
rm(gse178318_matrix)
gc()

gse245552_tar <- file.path(raw_dir, "GSE245552", "GSE245552_RAW.tar")
tar_members <- utils::untar(gse245552_tar, list = TRUE)
matrix_members <- sort(tar_members[grepl("_matrix[.]mtx[.]gz$", tar_members)])
sample_count_rows <- list()
sample_qc_rows <- list()
gse245552_temp <- tempfile("gse245552_")
dir.create(gse245552_temp)
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
  barcodes <- read_lines_gz(barcode_path)
  sample_label <- sub("^GSM[0-9]+_", "", prefix)
  gsm <- sub("_.*$", "", prefix)
  meta_row <- sample_meta[sample_meta$sample_id == gsm, , drop = FALSE]
  tissue_site <- if (nrow(meta_row) > 0) meta_row$tissue_site[1] else ""
  sample_title <- if (nrow(meta_row) > 0) meta_row$lesion_type[1] else ""
  sample_count_rows[[length(sample_count_rows) + 1]] <- data.frame(
    dataset_id = "GSE245552",
    gsm = gsm,
    sample_label = sample_label,
    tissue_site = tissue_site,
    sample_title = sample_title,
    genes = nrow(mat),
    cells = ncol(mat),
    barcodes = length(barcodes)
  )
  qc <- summarise_sparse_matrix(mat, safe_gene_symbols(features), rep(sample_label, ncol(mat)), "GSE245552", "sample_label")
  qc$gsm <- gsm
  qc$tissue_site <- tissue_site
  qc$sample_title <- sample_title
  sample_qc_rows[[length(sample_qc_rows) + 1]] <- qc
  rm(mat)
  unlink(c(matrix_path, barcode_path, feature_path))
  gc()
}
write_tsv(do.call(rbind, sample_count_rows), file.path(upload_dir, "gse245552_sample_cell_counts.tsv"))
write_tsv(do.call(rbind, sample_qc_rows), file.path(upload_dir, "gse245552_sample_qc_summary.tsv"))
pc_status <- c(pc_status, sprintf("GSE245552 sample matrices: %s", length(matrix_members)))

status_text <- c(
  "# CRLM CAF/ECM initial QC preflight",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  paste0("Raw input directory: `", raw_dir, "`"),
  "",
  pc_status,
  "",
  "No filtering, normalization, clustering, annotation, or final inclusion decision was performed.",
  "Only small QC and metadata summary tables were written."
)
writeLines(status_text, file.path(upload_dir, "STATUS.md"))
