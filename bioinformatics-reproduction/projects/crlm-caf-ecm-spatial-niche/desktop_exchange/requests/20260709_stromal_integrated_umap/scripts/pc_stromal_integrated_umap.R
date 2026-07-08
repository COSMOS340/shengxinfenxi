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

for (pkg in c("Seurat", "data.table", "ggplot2", "patchwork")) {
  install_if_missing(pkg)
}

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
  script_path <- normalizePath(
    "bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260709_stromal_integrated_umap/scripts/pc_stromal_integrated_umap.R",
    winslash = "/",
    mustWork = TRUE
  )
}

request_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- find_repo_root(request_dir)
project_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche")
cluster_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_stromal_subset_clustering_20260708")
integrated_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_stromal_integrated_umap_20260709")
upload_dir <- file.path(project_dir, "desktop_exchange", "uploads", "20260709_stromal_integrated_umap")
dir.create(integrated_object_dir, recursive = TRUE, showWarnings = FALSE)
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

object_paths <- data.table(
  dataset_id = c("GSE178318", "GSE245552"),
  object_path = file.path(
    cluster_object_dir,
    c("gse178318_stromal_subset_seurat_res0.4.rds", "gse245552_stromal_subset_seurat_res0.4.rds")
  )
)

load_one_object <- function(dataset_id, object_path) {
  if (!file.exists(object_path)) {
    stop("Missing Seurat object: ", object_path)
  }
  obj <- readRDS(object_path)
  if (!inherits(obj, "Seurat")) {
    stop("RDS is not a Seurat object: ", object_path)
  }
  meta <- obj@meta.data
  missing_fields <- setdiff(required_meta_fields, colnames(meta))
  if (length(missing_fields) > 0) {
    stop("Missing metadata fields in ", object_path, ": ", paste(missing_fields, collapse = ", "))
  }
  dataset_values <- unique(as.character(meta$dataset_id))
  if (!identical(dataset_values, dataset_id)) {
    stop("Unexpected dataset_id values in ", object_path, ": ", paste(dataset_values, collapse = ", "))
  }
  DefaultAssay(obj) <- "RNA"
  obj$source_main_cluster <- as.character(obj$main_cluster)
  obj$source_cluster_key <- paste(as.character(obj$dataset_id), as.character(obj$source_main_cluster), sep = "__")
  obj$source_cell_id <- colnames(obj)
  obj$main_cluster <- NULL
  obj
}

obj_list <- lapply(seq_len(nrow(object_paths)), function(i) {
  load_one_object(object_paths$dataset_id[i], object_paths$object_path[i])
})
names(obj_list) <- object_paths$dataset_id

object_audit <- rbindlist(lapply(names(obj_list), function(dataset_id) {
  obj <- obj_list[[dataset_id]]
  data.table(
    dataset_id = dataset_id,
    input_cells = ncol(obj),
    input_features = nrow(obj),
    metadata_fields_present = paste(intersect(c(required_meta_fields, optional_meta_fields, "source_main_cluster", "source_cluster_key"), colnames(obj@meta.data)), collapse = ",")
  )
}))

set.seed(20260709)
obj_list <- lapply(obj_list, function(obj) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  obj
})

integration_features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = integration_features, dims = 1:30, verbose = FALSE)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30, verbose = FALSE)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
dims_use <- 1:min(30, ncol(Embeddings(integrated, "pca")))
integrated <- RunUMAP(integrated, dims = dims_use, reduction = "pca", verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = dims_use, reduction = "pca", verbose = FALSE)
for (resolution_value in c(0.2, 0.4, 0.6)) {
  integrated <- FindClusters(integrated, resolution = resolution_value, verbose = FALSE)
}

cluster_cols <- grep("_snn_res[.]0[.]4$", colnames(integrated@meta.data), value = TRUE)
if (length(cluster_cols) == 0) {
  cluster_cols <- grep("res[.]0[.]4$", colnames(integrated@meta.data), value = TRUE)
}
if (length(cluster_cols) == 0) {
  integrated_cluster_col <- "seurat_clusters"
} else {
  integrated_cluster_col <- cluster_cols[1]
}
integrated$integrated_cluster <- as.character(integrated@meta.data[[integrated_cluster_col]])

integrated_object_path <- file.path(integrated_object_dir, "stromal_integrated_umap_seurat_res0.4.rds")
saveRDS(integrated, integrated_object_path)

umap <- Embeddings(integrated, reduction = "umap")
meta <- as.data.table(integrated@meta.data, keep.rownames = "cell_id")
meta <- meta[match(rownames(umap), cell_id)]
if (any(is.na(meta$cell_id))) {
  stop("Integrated UMAP cell names do not align with metadata")
}

export_fields <- c(
  "dataset_id",
  "integrated_cluster",
  "source_main_cluster",
  "source_cluster_key",
  "review_class",
  "sample_label",
  "tissue_site",
  "top_marker_program",
  "score_margin",
  "source_cell_id",
  intersect(optional_meta_fields, colnames(meta))
)
coord_dt <- cbind(
  data.table(cell_id = rownames(umap), integrated_umap_1 = as.numeric(umap[, 1]), integrated_umap_2 = as.numeric(umap[, 2])),
  meta[, ..export_fields]
)
setcolorder(coord_dt, c("dataset_id", "cell_id", "integrated_umap_1", "integrated_umap_2", setdiff(colnames(coord_dt), c("dataset_id", "cell_id", "integrated_umap_1", "integrated_umap_2"))))

cluster_centers <- coord_dt[, .(
  integrated_umap_1 = stats::median(integrated_umap_1),
  integrated_umap_2 = stats::median(integrated_umap_2),
  cells = .N
), by = .(integrated_cluster)]

cluster_counts <- coord_dt[, .N, by = .(integrated_cluster)]
cluster_counts[, fraction := N / sum(N)]
composition_dataset <- coord_dt[, .N, by = .(integrated_cluster, dataset_id)]
composition_review <- coord_dt[, .N, by = .(integrated_cluster, review_class)]
composition_program <- coord_dt[, .N, by = .(integrated_cluster, top_marker_program)]
composition_tissue <- coord_dt[, .N, by = .(integrated_cluster, tissue_site)]

audit <- rbind(
  object_audit,
  data.table(
    dataset_id = "integrated",
    input_cells = ncol(integrated),
    input_features = nrow(integrated),
    metadata_fields_present = paste(colnames(integrated@meta.data), collapse = ",")
  ),
  fill = TRUE
)
audit[, integration_features := length(integration_features)]
audit[, pca_dims_used := length(dims_use)]
audit[, integrated_cluster_column := integrated_cluster_col]
audit[, integrated_cluster_count := length(unique(coord_dt$integrated_cluster))]
audit[, integrated_object_path := integrated_object_path]
audit[, integrated_object_size_bytes := file_size_or_na(integrated_object_path)]

write_tsv(audit, file.path(upload_dir, "stromal_integrated_umap_audit.tsv"))
write_tsv(coord_dt, file.path(upload_dir, "stromal_integrated_umap_coordinates.tsv"))
write_tsv(cluster_centers, file.path(upload_dir, "stromal_integrated_umap_cluster_centers.tsv"))
write_tsv(cluster_counts, file.path(upload_dir, "stromal_integrated_cluster_counts.tsv"))
write_tsv(composition_dataset, file.path(upload_dir, "stromal_integrated_cluster_composition_by_dataset.tsv"))
write_tsv(composition_review, file.path(upload_dir, "stromal_integrated_cluster_composition_by_review_class.tsv"))
write_tsv(composition_program, file.path(upload_dir, "stromal_integrated_cluster_composition_by_top_marker_program.tsv"))
write_tsv(composition_tissue, file.path(upload_dir, "stromal_integrated_cluster_composition_by_tissue.tsv"))

source_group <- as.character(coord_dt$tissue_site)
source_group[is.na(source_group) | source_group == ""] <- as.character(coord_dt$dataset_id[is.na(source_group) | source_group == ""])
coord_dt[, source_group := source_group]

plot_umap <- function(dt, color_col, path, title_text, width = 8.0, height = 6.2) {
  p <- ggplot(dt, aes(x = integrated_umap_1, y = integrated_umap_2, color = .data[[color_col]])) +
    geom_point(size = 0.08, alpha = 0.72) +
    coord_equal() +
    labs(title = title_text, x = "Integrated UMAP 1", y = "Integrated UMAP 2", color = color_col) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.text = element_text(size = 6.5),
      legend.title = element_text(size = 7, face = "bold")
    )
  ggsave(path, p, width = width, height = height, dpi = 240, bg = "white")
}

plot_umap(coord_dt, "dataset_id", file.path(upload_dir, "stromal_integrated_umap_by_dataset.png"), "Integrated stromal UMAP by dataset")
plot_umap(coord_dt, "review_class", file.path(upload_dir, "stromal_integrated_umap_by_review_class.png"), "Integrated stromal UMAP by review class")
plot_umap(coord_dt, "integrated_cluster", file.path(upload_dir, "stromal_integrated_umap_by_integrated_cluster.png"), "Integrated stromal UMAP by integrated cluster", width = 9.2)
plot_umap(coord_dt, "top_marker_program", file.path(upload_dir, "stromal_integrated_umap_by_top_marker_program.png"), "Integrated stromal UMAP by top marker program", width = 9.2)

status_text <- c(
  "# CRLM CAF/ECM stromal integrated UMAP",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  paste0("PC-local input object directory: `", cluster_object_dir, "`"),
  paste0("PC-local integrated object directory: `", integrated_object_dir, "`"),
  paste0("Upload directory: `", upload_dir, "`"),
  "",
  paste0("Input datasets: ", paste(object_paths$dataset_id, collapse = ", ")),
  paste0("Cells integrated: ", nrow(coord_dt)),
  paste0("Integration features: ", length(integration_features)),
  paste0("PCA dimensions used: ", length(dims_use)),
  paste0("Integrated cluster column: ", integrated_cluster_col),
  paste0("Integrated clusters at resolution 0.4: ", length(unique(coord_dt$integrated_cluster))),
  "",
  "This request creates one shared integrated UMAP for atlas-style visualization.",
  "It does not finalize cell labels or remove carryover groups.",
  "The integrated Seurat RDS remains PC-local and should not be committed."
)
writeLines(status_text, file.path(upload_dir, "STATUS.md"))

message("Wrote integrated UMAP export to: ", upload_dir)
