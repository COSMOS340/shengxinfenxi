#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

input_rds <- Sys.getenv(
  "TNBC_PORTABLE_INPUT_RDS",
  unset = "outputs/04_data_processed/single_cell/full_scope_standard_objects/set2_TNBC_standard_myeloid.rds"
)
output_rds <- Sys.getenv(
  "TNBC_PORTABLE_OUTPUT_RDS",
  unset = "outputs/04_data_processed/single_cell/full_scope_standard_objects/set2_TNBC_standard_myeloid_pre_supp_fig_s1_portable.rds"
)
table_dir <- Sys.getenv(
  "TNBC_PORTABLE_TABLE_DIR",
  unset = "../desktop_exchange/uploads/20260702_tnbc_pre_supp_fig_s1_portable/tables"
)
log_path <- Sys.getenv(
  "TNBC_PORTABLE_LOG_PATH",
  unset = "../desktop_exchange/uploads/20260702_tnbc_pre_supp_fig_s1_portable/logs/export_tnbc_pre_supp_fig_s1_portable.log"
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

write_log <- function(step, status, note = "") {
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    step = step,
    status = status,
    note = note,
    stringsAsFactors = FALSE
  )
  write.table(row, log_path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = !file.exists(log_path), append = file.exists(log_path))
}

fig_s1_labels <- c(
  "Pre_P001", "Pre_P002", "Pre_P004", "Pre_P005", "Pre_P007", "Pre_P008",
  "Pre_P010", "Pre_P011", "Pre_P012", "Pre_P013", "Pre_P014", "Pre_P016",
  "Pre_P017", "Pre_P018", "Pre_P019", "Pre_P020", "Pre_P022", "Pre_P023",
  "Pre_P024", "Pre_P025", "Pre_P028"
)
fig_s1_patient_ids <- sub("^Pre_", "", fig_s1_labels)
expression_source_value <- "standard_seurat_object_pc_handoff_scope_refined_to_supp_fig_s1_pre"
required_metadata_columns <- c(
  "paper_dataset_id", "paper_label", "patient_or_subject", "origin_or_tissue",
  "source_cell_id", "manual_annotation", "manual_broad_lineage", "batch_id"
)
required_tnbc_columns <- c("patient_from_title", "timepoint_from_title", "batch_id", "source_cell_id")
allowed_lineages <- c("Macrophages", "Monocytes", "Classical DCs", "Plasmacytoid DCs", "Mast cells")

to_dgCMatrix <- function(mat, layer_name) {
  write_log("materialize_layer_start", "running", layer_name)
  sparse <- as(mat, "dgCMatrix")
  if (!inherits(sparse, "dgCMatrix")) {
    sparse <- as(sparse, "dgCMatrix")
  }
  write_log("materialize_layer_done", "complete", paste0(layer_name, "; rows=", nrow(sparse), "; cols=", ncol(sparse), "; nnz=", length(sparse@x)))
  sparse
}

validate_standard_object <- function(obj) {
  missing_columns <- setdiff(required_metadata_columns, names(obj@meta.data))
  if (length(missing_columns) > 0L) {
    stop("Missing required metadata columns: ", paste(missing_columns, collapse = ", "))
  }
  for (column in required_metadata_columns) {
    values <- as.character(obj@meta.data[[column]])
    if (anyNA(values) || any(!nzchar(values))) {
      stop("Metadata column contains empty values: ", column)
    }
  }
  outside_lineages <- setdiff(sort(unique(as.character(obj$manual_broad_lineage))), allowed_lineages)
  if (length(outside_lineages) > 0L) {
    stop("Unexpected manual_broad_lineage values: ", paste(outside_lineages, collapse = ", "))
  }
  duplicate_index <- anyDuplicated(as.character(obj$source_cell_id))
  if (duplicate_index > 0L) {
    stop("Duplicated source_cell_id at position: ", duplicate_index)
  }
}

write_log("script_start", "running", paste0("input=", input_rds, "; output=", output_rds))
if (!file.exists(input_rds)) {
  stop("Missing input RDS: ", input_rds)
}

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input RDS does not contain a Seurat object.")
}
DefaultAssay(obj) <- "RNA"
missing_tnbc_columns <- setdiff(required_tnbc_columns, names(obj@meta.data))
if (length(missing_tnbc_columns) > 0L) {
  stop("Input TNBC object missing columns: ", paste(missing_tnbc_columns, collapse = ", "))
}

meta <- as.data.table(obj@meta.data, keep.rownames = "cell_name")
keep_cells <- meta[
  patient_from_title %in% fig_s1_patient_ids & timepoint_from_title == "Pre",
  cell_name
]
if (length(keep_cells) == 0L) {
  stop("Supplementary Fig. S1 TNBC filter retained zero cells.")
}
filtered <- subset(obj, cells = keep_cells)
filtered$paper_fig_s1_patient_label <- paste0("Pre_", as.character(filtered$patient_from_title))
filtered$paper_scope_source <- "supplementary_fig_s1_patient_legend_pre_p_labels"
filtered$expression_source <- expression_source_value
validate_standard_object(filtered)

retained_patient_ids <- sort(unique(as.character(filtered$patient_from_title)))
missing_patient_ids <- setdiff(fig_s1_patient_ids, retained_patient_ids)
extra_patient_ids <- setdiff(retained_patient_ids, fig_s1_patient_ids)
if (length(missing_patient_ids) > 0L || length(extra_patient_ids) > 0L) {
  stop(
    "Retained patient IDs do not match Fig. S1 labels. Missing: ",
    paste(missing_patient_ids, collapse = ";"),
    "; extra: ",
    paste(extra_patient_ids, collapse = ";")
  )
}

counts_sparse <- to_dgCMatrix(GetAssayData(filtered, assay = "RNA", layer = "counts"), "counts")
data_sparse <- to_dgCMatrix(GetAssayData(filtered, assay = "RNA", layer = "data"), "data")
if (!identical(rownames(counts_sparse), rownames(data_sparse)) || !identical(colnames(counts_sparse), colnames(data_sparse))) {
  stop("Materialized counts and data layers do not have identical dimnames.")
}

portable_meta <- filtered@meta.data[colnames(counts_sparse), , drop = FALSE]
portable <- CreateSeuratObject(counts = counts_sparse, meta.data = portable_meta, project = "TNBC")
portable <- SetAssayData(portable, assay = "RNA", layer = "data", new.data = data_sparse)
DefaultAssay(portable) <- "RNA"
validate_standard_object(portable)
saveRDS(portable, output_rds)
write_log("save_rds", "complete", output_rds)

check_obj <- readRDS(output_rds)
DefaultAssay(check_obj) <- "RNA"
layer_check <- rbindlist(lapply(c("counts", "data"), function(layer_name) {
  layer_mat <- GetAssayData(check_obj, assay = "RNA", layer = layer_name)
  class_text <- paste(class(layer_mat), collapse = ";")
  small_slice_status <- tryCatch({
    invisible(as.matrix(layer_mat[1:min(3, nrow(layer_mat)), 1:min(3, ncol(layer_mat))]))
    "ok"
  }, error = function(e) paste("error:", conditionMessage(e)))
  data.table(
    layer = layer_name,
    matrix_class = class_text,
    rows = nrow(layer_mat),
    cols = ncol(layer_mat),
    contains_bpcells_class = grepl("BPCells|Iterable|RenameDims", class_text),
    small_slice_status = small_slice_status
  )
}))
fwrite(layer_check, file.path(table_dir, "set2_TNBC_pre_supp_fig_s1_portable_layer_check.tsv"), sep = "\t", quote = FALSE)
if (any(layer_check$contains_bpcells_class) || any(layer_check$small_slice_status != "ok")) {
  stop("Saved portable RDS failed layer materialization check.")
}

summary <- data.table(
  metric = c(
    "dataset_id",
    "paper_label",
    "retained_cells",
    "genes",
    "retained_patients",
    "retained_batches",
    "lineages",
    "expression_source",
    "output_rds"
  ),
  value = c(
    "set2_TNBC",
    "TNBC",
    as.character(ncol(check_obj)),
    as.character(nrow(check_obj)),
    as.character(uniqueN(as.character(check_obj$patient_or_subject))),
    as.character(uniqueN(as.character(check_obj$batch_id))),
    paste(sort(unique(as.character(check_obj$manual_broad_lineage))), collapse = ";"),
    paste(sort(unique(as.character(check_obj$expression_source))), collapse = ";"),
    output_rds
  )
)
fwrite(summary, file.path(table_dir, "set2_TNBC_pre_supp_fig_s1_portable_summary.tsv"), sep = "\t", quote = FALSE)

if (ncol(check_obj) != 30209L) {
  stop("Expected 30209 retained cells; observed ", ncol(check_obj))
}
if (uniqueN(as.character(check_obj$patient_or_subject)) != 21L) {
  stop("Expected 21 retained patients.")
}
if (uniqueN(as.character(check_obj$batch_id)) != 36L) {
  stop("Expected 36 retained batches.")
}

write_log("script_done", "complete", paste0("cells=", ncol(check_obj), "; genes=", nrow(check_obj), "; patients=21; batches=36"))
message("Portable TNBC Supplementary Fig. S1 pre-only standard object complete.")
