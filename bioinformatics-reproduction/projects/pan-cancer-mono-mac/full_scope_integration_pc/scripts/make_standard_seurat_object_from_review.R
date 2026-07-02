#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

input_rds <- Sys.getenv("STANDARDIZE_INPUT_RDS", unset = "")
annotation_path <- Sys.getenv("STANDARDIZE_CLUSTER_ANNOTATION", unset = "")
output_rds <- Sys.getenv("STANDARDIZE_OUTPUT_RDS", unset = "")
dataset_id <- Sys.getenv("STANDARDIZE_DATASET_ID", unset = "")
paper_label <- Sys.getenv("STANDARDIZE_PAPER_LABEL", unset = "")
cluster_column <- Sys.getenv("STANDARDIZE_CLUSTER_COLUMN", unset = "seurat_clusters")
source_cell_column <- Sys.getenv("STANDARDIZE_SOURCE_CELL_COLUMN", unset = "colnames")
patient_column <- Sys.getenv("STANDARDIZE_PATIENT_COLUMN", unset = "")
origin_column <- Sys.getenv("STANDARDIZE_ORIGIN_COLUMN", unset = "")
batch_column <- Sys.getenv("STANDARDIZE_BATCH_COLUMN", unset = "")

required_env <- c(
  STANDARDIZE_INPUT_RDS = input_rds,
  STANDARDIZE_CLUSTER_ANNOTATION = annotation_path,
  STANDARDIZE_OUTPUT_RDS = output_rds,
  STANDARDIZE_DATASET_ID = dataset_id,
  STANDARDIZE_PAPER_LABEL = paper_label,
  STANDARDIZE_PATIENT_COLUMN = patient_column,
  STANDARDIZE_ORIGIN_COLUMN = origin_column,
  STANDARDIZE_BATCH_COLUMN = batch_column
)
missing_env <- names(required_env)[!nzchar(required_env)]
if (length(missing_env) > 0L) {
  stop("Missing required environment variables: ", paste(missing_env, collapse = ", "))
}

if (!file.exists(input_rds)) {
  stop("Missing input RDS: ", input_rds)
}
if (!file.exists(annotation_path)) {
  stop("Missing cluster annotation TSV: ", annotation_path)
}

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("STANDARDIZE_INPUT_RDS does not contain a Seurat object.")
}
if (!"RNA" %in% names(obj@assays)) {
  stop("Input Seurat object lacks RNA assay.")
}
DefaultAssay(obj) <- "RNA"

annotation <- fread(annotation_path)
required_annotation_columns <- c("cluster", "retain_for_full_integration", "manual_annotation", "manual_broad_lineage", "review_notes")
missing_annotation_columns <- setdiff(required_annotation_columns, names(annotation))
if (length(missing_annotation_columns) > 0L) {
  stop("Cluster annotation table missing columns: ", paste(missing_annotation_columns, collapse = ", "))
}
unexpected_retain_values <- setdiff(unique(annotation$retain_for_full_integration), c("TRUE", "FALSE"))
if (length(unexpected_retain_values) > 0L) {
  stop("retain_for_full_integration must contain only TRUE or FALSE. Observed: ", paste(unexpected_retain_values, collapse = ", "))
}
allowed_lineages <- c("Macrophages", "Monocytes", "Classical DCs", "Plasmacytoid DCs", "Mast cells")
retained_annotation <- annotation[retain_for_full_integration == "TRUE"]
if (nrow(retained_annotation) == 0L) {
  stop("Cluster annotation table retains zero clusters.")
}
if (anyNA(retained_annotation$manual_broad_lineage) || any(!nzchar(as.character(retained_annotation$manual_broad_lineage)))) {
  stop("Retained rows must have non-empty manual_broad_lineage values.")
}
outside_lineages <- setdiff(unique(retained_annotation$manual_broad_lineage), allowed_lineages)
if (length(outside_lineages) > 0L) {
  stop("manual_broad_lineage contains values outside the allowed set: ", paste(outside_lineages, collapse = ", "))
}
if (anyDuplicated(as.character(annotation$cluster)) > 0L) {
  stop("Cluster annotation table contains duplicated cluster values.")
}

required_metadata_columns <- c(cluster_column, patient_column, origin_column, batch_column)
if (!identical(source_cell_column, "colnames")) {
  required_metadata_columns <- c(required_metadata_columns, source_cell_column)
}
missing_metadata_columns <- setdiff(required_metadata_columns, names(obj@meta.data))
if (length(missing_metadata_columns) > 0L) {
  stop("Input Seurat metadata missing columns: ", paste(missing_metadata_columns, collapse = ", "))
}

meta <- as.data.table(obj@meta.data, keep.rownames = "cell_name")
meta[, cluster_for_review := as.character(get(cluster_column))]
annotation[, cluster_for_review := as.character(cluster)]
meta <- merge(
  meta,
  annotation[, .(cluster_for_review, retain_for_full_integration, manual_annotation, manual_broad_lineage, review_notes)],
  by = "cluster_for_review",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(meta$retain_for_full_integration)) {
  missing_clusters <- sort(unique(meta[is.na(retain_for_full_integration), cluster_for_review]))
  stop("Input object has clusters not present in annotation table: ", paste(missing_clusters, collapse = ", "))
}

retained_cells <- meta[retain_for_full_integration == "TRUE", cell_name]
if (length(retained_cells) == 0L) {
  stop("No cells retained for full integration.")
}
obj <- subset(obj, cells = retained_cells)
meta <- meta[cell_name %in% retained_cells]
meta <- meta[match(colnames(obj), cell_name)]
if (!identical(meta$cell_name, colnames(obj))) {
  stop("Metadata order did not match Seurat cell order after retention.")
}

source_cell_id <- if (identical(source_cell_column, "colnames")) {
  colnames(obj)
} else {
  as.character(meta[[source_cell_column]])
}
if (anyNA(source_cell_id) || any(!nzchar(source_cell_id))) {
  stop("source_cell_id contains empty values.")
}

obj$paper_dataset_id <- dataset_id
obj$paper_label <- paper_label
obj$patient_or_subject <- as.character(meta[[patient_column]])
obj$origin_or_tissue <- as.character(meta[[origin_column]])
obj$source_cell_id <- source_cell_id
obj$manual_annotation <- as.character(meta$manual_annotation)
obj$manual_broad_lineage <- as.character(meta$manual_broad_lineage)
obj$batch_id <- as.character(meta[[batch_column]])
obj$review_notes <- as.character(meta$review_notes)

standard_columns <- c(
  "paper_dataset_id", "paper_label", "patient_or_subject", "origin_or_tissue",
  "source_cell_id", "manual_annotation", "manual_broad_lineage", "batch_id"
)
for (column in standard_columns) {
  if (anyNA(obj@meta.data[[column]]) || any(!nzchar(as.character(obj@meta.data[[column]])))) {
    stop("Standard metadata column contains empty values: ", column)
  }
}
if (anyDuplicated(obj$source_cell_id) > 0L) {
  stop("Retained object has duplicated source_cell_id values.")
}

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, output_rds)

summary <- data.table(
  metric = c("dataset_id", "paper_label", "retained_cells", "genes", "clusters_retained", "output_rds"),
  value = c(dataset_id, paper_label, as.character(ncol(obj)), as.character(nrow(obj)), as.character(uniqueN(obj[[cluster_column]][, 1])), output_rds)
)
summary_path <- sub("\\.rds$", "_standard_object_summary.tsv", output_rds)
fwrite(summary, summary_path, sep = "\t")
message("Standard Seurat object written: ", output_rds)
