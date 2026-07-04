#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg) == 0L) stop("Cannot resolve script path from commandArgs().")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
request_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

output_dir <- Sys.getenv(
  "FIG4_OUTPUT_DIR",
  unset = file.path(getwd(), "desktop_exchange", "uploads", "20260704_fig4a_gse306459_formal_snatac", "pc_formal_fig4a_results")
)

required_table <- fread(file.path(request_dir, "inputs", "fig4a_formal_required_outputs.tsv"))
required_table <- required_table[required == "yes"]
required_table[, local_path := file.path(output_dir, relative_path)]
required_table[, exists := file.exists(local_path)]
required_table[, size_bytes := fifelse(exists, file.size(local_path), NA_real_)]
required_table[, nonempty := exists & size_bytes > 0]

dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
fwrite(
  required_table[, .(relative_path, exists, size_bytes, nonempty, reason)],
  file.path(output_dir, "tables", "fig4a_formal_required_file_check.tsv"),
  sep = "\t",
  quote = FALSE
)

if (any(!required_table$nonempty)) {
  missing <- required_table[!nonempty, relative_path]
  stop("Missing or empty required outputs: ", paste(missing, collapse = "; "))
}

umap_path <- file.path(output_dir, "tables", "fig4a_snatac_umap_labels.tsv.gz")
label_path <- file.path(output_dir, "tables", "fig4a_label_counts.tsv")
score_path <- file.path(output_dir, "tables", "fig4a_prediction_score_summary_by_label.tsv")
qc_path <- file.path(output_dir, "tables", "fig4a_qc_metrics_by_cell.tsv.gz")
threshold_path <- file.path(output_dir, "tables", "fig4a_qc_thresholds.tsv")

umap_dt <- fread(umap_path)
label_dt <- fread(label_path)
score_dt <- fread(score_path)
qc_dt <- fread(qc_path, nrows = 5)
threshold_dt <- fread(threshold_path)

required_umap_cols <- c("cell_id", "UMAP_1", "UMAP_2", "Library", "tissue_group", "predicted_label", "prediction.score.max")
missing_umap_cols <- setdiff(required_umap_cols, names(umap_dt))
if (length(missing_umap_cols) > 0L) {
  stop("UMAP table missing columns: ", paste(missing_umap_cols, collapse = ", "))
}

summary_dt <- data.table(
  check = c(
    "required_files_nonempty",
    "umap_rows",
    "transferred_labels",
    "libraries",
    "qc_metric_columns",
    "qc_threshold_rows",
    "score_rows"
  ),
  value = c(
    nrow(required_table),
    nrow(umap_dt),
    uniqueN(umap_dt$predicted_label),
    uniqueN(umap_dt$Library),
    paste(names(qc_dt), collapse = ";"),
    nrow(threshold_dt),
    nrow(score_dt)
  )
)
fwrite(summary_dt, file.path(output_dir, "tables", "fig4a_formal_postrun_check.tsv"), sep = "\t", quote = FALSE)

message("PASS: formal Fig4A upload check completed for ", output_dir)
