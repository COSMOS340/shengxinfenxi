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

for (pkg in c("Seurat", "data.table", "ggplot2")) {
  install_if_missing(pkg)
}

library(Seurat)
library(data.table)
library(ggplot2)

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
    "bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_stromal_umap_coordinate_export/scripts/pc_export_stromal_umap_coordinates.R",
    winslash = "/",
    mustWork = TRUE
  )
}

request_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- find_repo_root(request_dir)
project_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche")
cluster_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_stromal_subset_clustering_20260708")
upload_dir <- file.path(project_dir, "desktop_exchange", "uploads", "20260708_stromal_umap_coordinate_export")
dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  data.table::fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

object_paths <- data.table(
  dataset_id = c("GSE178318", "GSE245552"),
  object_path = file.path(
    cluster_object_dir,
    c("gse178318_stromal_subset_seurat_res0.4.rds", "gse245552_stromal_subset_seurat_res0.4.rds")
  )
)

required_meta_fields <- c(
  "dataset_id",
  "main_cluster",
  "review_class",
  "sample_label",
  "tissue_site",
  "top_marker_program",
  "score_margin"
)
optional_meta_fields <- c("object_label", "raw_barcode", "nCount_RNA", "nFeature_RNA", "percent.mt")

export_one <- function(dataset_id, object_path) {
  if (!file.exists(object_path)) {
    stop("Missing Seurat object: ", object_path)
  }
  obj <- readRDS(object_path)
  if (!inherits(obj, "Seurat")) {
    stop("RDS is not a Seurat object: ", object_path)
  }
  if (!"umap" %in% names(obj@reductions)) {
    stop("Missing umap reduction in: ", object_path)
  }
  umap <- Seurat::Embeddings(obj, reduction = "umap")
  if (ncol(umap) < 2) {
    stop("UMAP reduction has fewer than two dimensions in: ", object_path)
  }
  meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
  missing_fields <- setdiff(required_meta_fields, names(meta))
  if (length(missing_fields) > 0) {
    stop("Missing metadata fields in ", object_path, ": ", paste(missing_fields, collapse = ", "))
  }
  meta <- meta[match(rownames(umap), cell_id)]
  if (any(is.na(meta$cell_id))) {
    stop("UMAP cell names do not align with metadata rownames in: ", object_path)
  }
  export_fields <- c(required_meta_fields, intersect(optional_meta_fields, names(meta)))
  out <- cbind(
    data.table(cell_id = rownames(umap), umap_1 = as.numeric(umap[, 1]), umap_2 = as.numeric(umap[, 2])),
    meta[, ..export_fields]
  )
  out[, dataset_id := as.character(dataset_id)]
  setcolorder(out, c("dataset_id", "cell_id", "umap_1", "umap_2", setdiff(names(out), c("dataset_id", "cell_id", "umap_1", "umap_2"))))
  center_dt <- out[, .(
    umap_1 = stats::median(umap_1),
    umap_2 = stats::median(umap_2),
    cells = .N
  ), by = .(dataset_id, main_cluster)]
  audit_dt <- data.table(
    dataset_id = dataset_id,
    object_path = object_path,
    cells = nrow(out),
    umap_columns = paste(colnames(umap)[1:2], collapse = ","),
    required_fields = paste(required_meta_fields, collapse = ","),
    optional_fields_present = paste(intersect(optional_meta_fields, names(meta)), collapse = ","),
    main_cluster_count = length(unique(out$main_cluster)),
    review_class_count = length(unique(out$review_class))
  )
  list(coordinates = out, centers = center_dt, audit = audit_dt)
}

results <- lapply(seq_len(nrow(object_paths)), function(i) {
  export_one(object_paths$dataset_id[i], object_paths$object_path[i])
})

coordinates <- rbindlist(lapply(results, `[[`, "coordinates"), fill = TRUE)
centers <- rbindlist(lapply(results, `[[`, "centers"), fill = TRUE)
audit <- rbindlist(lapply(results, `[[`, "audit"), fill = TRUE)

write_tsv(coordinates, file.path(upload_dir, "stromal_subset_umap_coordinates.tsv"))
write_tsv(centers, file.path(upload_dir, "stromal_subset_umap_cluster_centers.tsv"))
write_tsv(audit, file.path(upload_dir, "stromal_subset_umap_export_audit.tsv"))

plot_umap <- function(dt, color_col, path, title_text) {
  p <- ggplot(dt, aes(x = umap_1, y = umap_2, color = .data[[color_col]])) +
    geom_point(size = 0.08, alpha = 0.75) +
    facet_wrap(~dataset_id, scales = "free") +
    coord_equal() +
    labs(title = title_text, x = "UMAP 1", y = "UMAP 2", color = color_col) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold"),
      legend.text = element_text(size = 6.5),
      legend.title = element_text(size = 7, face = "bold")
    )
  ggsave(path, p, width = 9.5, height = 5.5, dpi = 240, bg = "white")
}

plot_umap(
  coordinates,
  "main_cluster",
  file.path(upload_dir, "stromal_subset_umap_by_cluster_preview.png"),
  "CRLM stromal subset UMAP by cluster"
)
plot_umap(
  coordinates,
  "review_class",
  file.path(upload_dir, "stromal_subset_umap_by_review_class_preview.png"),
  "CRLM stromal subset UMAP by review class"
)
plot_umap(
  coordinates,
  "dataset_id",
  file.path(upload_dir, "stromal_subset_umap_by_dataset_preview.png"),
  "CRLM stromal subset UMAP by dataset"
)

status_text <- c(
  "# CRLM CAF/ECM stromal UMAP coordinate export",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  paste0("PC-local clustering object directory: `", cluster_object_dir, "`"),
  paste0("Upload directory: `", upload_dir, "`"),
  "",
  paste0("Datasets exported: ", paste(audit$dataset_id, collapse = ", ")),
  paste0("Cells exported: ", nrow(coordinates)),
  paste0("Cluster center rows: ", nrow(centers)),
  "",
  "This request only exports UMAP coordinates and metadata from prior review objects.",
  "It does not rerun clustering, alter identities, or finalize biological labels."
)
writeLines(status_text, file.path(upload_dir, "STATUS.md"))

message("Wrote UMAP coordinate export to: ", upload_dir)
