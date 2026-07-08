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

install_if_missing("data.table")
library(data.table)

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
  script_path <- normalizePath(
    "bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_stromal_marker_table_fix/scripts/pc_fix_stromal_marker_tables.R",
    winslash = "/",
    mustWork = TRUE
  )
}

request_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- find_repo_root(request_dir)
project_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche")
cluster_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_stromal_subset_clustering_20260708")
upload_dir <- file.path(project_dir, "desktop_exchange", "uploads", "20260708_stromal_marker_table_fix")
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

score_column <- function(markers) {
  if ("avg_log2FC" %in% names(markers)) {
    return("avg_log2FC")
  }
  if ("avg_logFC" %in% names(markers)) {
    return("avg_logFC")
  }
  NULL
}

top_n_by_cluster <- function(markers, n) {
  markers <- as.data.table(markers)
  if (nrow(markers) == 0) {
    markers[, marker_rank := integer()]
    return(markers)
  }
  sc <- score_column(markers)
  markers[, cluster_sort := suppressWarnings(as.integer(cluster))]
  if (is.null(sc)) {
    setorder(markers, dataset_id, cluster_sort, gene)
  } else if ("pct.1" %in% names(markers)) {
    setorderv(markers, c("dataset_id", "cluster_sort", sc, "pct.1"), c(1, 1, -1, -1))
  } else {
    setorderv(markers, c("dataset_id", "cluster_sort", sc), c(1, 1, -1))
  }
  markers[, marker_rank := seq_len(.N), by = .(dataset_id, cluster)]
  markers <- markers[marker_rank <= n]
  markers[, cluster_sort := NULL]
  markers
}

read_dataset_markers <- function(dataset_id) {
  dataset_id_value <- dataset_id
  marker_path <- file.path(
    cluster_object_dir,
    paste0(tolower(dataset_id_value), "_stromal_subset_all_markers_res0.4.rds")
  )
  if (!file.exists(marker_path)) {
    stop("Missing PC-local marker RDS: ", marker_path)
  }
  markers <- as.data.table(readRDS(marker_path))
  required_fields <- c("cluster", "gene")
  missing_fields <- setdiff(required_fields, names(markers))
  if (length(missing_fields) > 0) {
    stop("Marker table missing fields for ", dataset_id_value, ": ", paste(missing_fields, collapse = ", "))
  }
  markers[, dataset_id := dataset_id_value]
  markers[, cluster := as.character(cluster)]
  setcolorder(markers, c("dataset_id", setdiff(names(markers), "dataset_id")))
  markers
}

dataset_ids <- c("GSE178318", "GSE245552")
markers_all <- rbindlist(lapply(dataset_ids, read_dataset_markers), fill = TRUE)

if (!all(c("dataset_id", "cluster", "gene") %in% names(markers_all))) {
  stop("Repaired marker table does not contain dataset_id, cluster, and gene")
}

top10 <- top_n_by_cluster(copy(markers_all), 10)
top50 <- top_n_by_cluster(copy(markers_all), 50)

if (any(is.na(top10$dataset_id)) || any(top10$dataset_id == "")) {
  stop("Top10 marker table has empty dataset_id values")
}
if (any(is.na(top50$dataset_id)) || any(top50$dataset_id == "")) {
  stop("Top50 marker table has empty dataset_id values")
}

summary_dt <- copy(markers_all)
summary_dt[, cluster_sort := suppressWarnings(as.integer(cluster))]
sc <- score_column(summary_dt)
if (is.null(sc)) {
  setorder(summary_dt, dataset_id, cluster_sort, gene)
} else if ("pct.1" %in% names(summary_dt)) {
  setorderv(summary_dt, c("dataset_id", "cluster_sort", sc, "pct.1"), c(1, 1, -1, -1))
} else {
  setorderv(summary_dt, c("dataset_id", "cluster_sort", sc), c(1, 1, -1))
}
summary_dt[, top_gene := gene[1], by = .(dataset_id, cluster)]
marker_count_summary <- summary_dt[, .(
  marker_rows = .N,
  top_gene = top_gene[1],
  max_score = if (is.null(sc)) NA_real_ else max(get(sc), na.rm = TRUE)
), by = .(dataset_id, cluster)]
marker_count_summary[, cluster_sort := suppressWarnings(as.integer(cluster))]
setorder(marker_count_summary, dataset_id, cluster_sort)
marker_count_summary[, cluster_sort := NULL]

manifest <- data.frame(
  dataset_id = dataset_ids,
  all_marker_rds_path = file.path(
    cluster_object_dir,
    paste0(tolower(dataset_ids), "_stromal_subset_all_markers_res0.4.rds")
  ),
  all_marker_rds_size_bytes = vapply(
    file.path(cluster_object_dir, paste0(tolower(dataset_ids), "_stromal_subset_all_markers_res0.4.rds")),
    file_size_or_na,
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

write_tsv(marker_count_summary, file.path(upload_dir, "stromal_subset_cluster_marker_count_summary.tsv"))
write_tsv(top10, file.path(upload_dir, "stromal_subset_cluster_markers_top10_by_dataset.tsv"))
write_tsv(top50, file.path(upload_dir, "stromal_subset_cluster_markers_top50_by_dataset.tsv"))
write_tsv(manifest, file.path(upload_dir, "stromal_subset_marker_table_fix_manifest.tsv"))

status_text <- c(
  "# CRLM CAF/ECM Stromal Marker Table Fix",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  paste0("PC-local marker RDS directory: `", cluster_object_dir, "`"),
  paste0("Datasets repaired: ", paste(dataset_ids, collapse = ", ")),
  paste0("Full marker rows read: ", nrow(markers_all)),
  paste0("Top10 marker rows written: ", nrow(top10)),
  paste0("Top50 marker rows written: ", nrow(top50)),
  paste0("Marker count summary rows written: ", nrow(marker_count_summary)),
  "",
  "This repair only adds explicit dataset_id support to marker review tables.",
  "It does not rerun clustering, marker testing, filtering, or final cell labeling.",
  "PC-local RDS files should not be committed to GitHub."
)
writeLines(status_text, file.path(upload_dir, "STATUS.md"))
