suppressPackageStartupMessages(library(data.table))

out <- Sys.getenv("FIG4_OUTPUT_DIR", unset = "fig4_snatac/outputs_fig4g_hi_all")
tab <- file.path(out, "tables")
fig <- file.path(out, "figures")
logs <- file.path(out, "logs")
files <- file.path(out, "files")

required <- c(
  file.path(logs, "fig4_gse306459_snatac.log"),
  file.path(tab, "fig4_gse306459_run_summary.tsv"),
  file.path(tab, "fig4_motif_lookup.tsv"),
  file.path(tab, "fig4_requested_motif_ids.tsv"),
  file.path(tab, "fig4e_chromvar_motif_activity_by_label.tsv"),
  file.path(tab, "fig4g_track_region_summary.tsv"),
  file.path(fig, "fig4g_gene_coverage_tracks.png"),
  file.path(fig, "fig4g_gene_coverage_tracks.pdf"),
  file.path(fig, "fig4h_motif_accessibility_umaps.png"),
  file.path(fig, "fig4h_motif_accessibility_umaps.pdf"),
  file.path(fig, "fig4i_motif_logos.png"),
  file.path(fig, "fig4i_motif_logos.pdf"),
  file.path(files, "fig4_gse306459_snatac_processed.rds.not_saved.txt")
)

required_info <- file.info(required)
missing <- required[!file.exists(required) | required_info$size <= 0]
if (length(missing)) {
  stop("Missing or empty required files: ", paste(missing, collapse = "; "))
}

summary_dt <- fread(file.path(tab, "fig4_gse306459_run_summary.tsv"))
motifs <- fread(file.path(tab, "fig4_requested_motif_ids.tsv"))
motif_activity <- fread(file.path(tab, "fig4e_chromvar_motif_activity_by_label.tsv"))
track_summary <- fread(file.path(tab, "fig4g_track_region_summary.tsv"))

required_summary_metrics <- c(
  "run_scope",
  "cells",
  "peaks",
  "libraries",
  "samples",
  "transferred_labels",
  "transfer_features",
  "differential_peak_rows",
  "requested_motifs_found"
)
missing_metrics <- setdiff(required_summary_metrics, summary_dt$metric)
if (length(missing_metrics)) {
  stop("Missing run summary metrics: ", paste(missing_metrics, collapse = ", "))
}

run_scope <- summary_dt[metric == "run_scope", value][1]
if (!identical(run_scope, "all")) stop("run_scope is not all: ", run_scope)

requested_motifs_found <- as.integer(summary_dt[metric == "requested_motifs_found", value][1])
if (!identical(requested_motifs_found, 6L)) {
  stop("requested_motifs_found is not 6: ", requested_motifs_found)
}
if (nrow(motifs) != 6L) stop("fig4_requested_motif_ids.tsv row count is not 6: ", nrow(motifs))
if (nrow(track_summary) != 4L) stop("fig4g_track_region_summary.tsv row count is not 4: ", nrow(track_summary))
if (nrow(motif_activity) <= 0L) stop("fig4e_chromvar_motif_activity_by_label.tsv is empty")

checks <- data.table(
  check = c(
    "required_files_nonempty",
    "run_scope",
    "requested_motifs_found",
    "requested_motif_rows",
    "motif_activity_rows",
    "track_gene_rows",
    "processed_rds_policy"
  ),
  result = "PASS",
  value = c(
    length(required),
    run_scope,
    requested_motifs_found,
    nrow(motifs),
    nrow(motif_activity),
    nrow(track_summary),
    "not_saved"
  )
)
fwrite(checks, file.path(tab, "fig4_postrun_checks_all.tsv"), sep = "\t", quote = FALSE)

cat("fig4_g_h_i_all_checks_passed\n")
print(checks)
