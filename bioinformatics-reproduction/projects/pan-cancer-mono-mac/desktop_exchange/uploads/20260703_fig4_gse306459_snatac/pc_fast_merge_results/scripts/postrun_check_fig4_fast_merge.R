suppressPackageStartupMessages(library(data.table))

out <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-cancer-mono-mac/fig4_snatac/outputs_fast_merge"
tab <- file.path(out, "tables")
fig <- file.path(out, "figures")
logs <- file.path(out, "logs")
files <- file.path(out, "files")

required <- c(
  file.path(logs, "fig4_gse306459_snatac.log"),
  file.path(tab, "fig4_gse306459_geo_download_manifest.tsv"),
  file.path(tab, "fig4_table_s1_barcode_counts_by_library.tsv"),
  file.path(tab, "fig4a_snatac_umap_labels.tsv.gz"),
  file.path(tab, "fig4b_label_tissue_proportions.tsv"),
  file.path(tab, "fig4b_wilcoxon_by_label.tsv"),
  file.path(tab, "fig4_differential_peaks_by_label.tsv.gz"),
  file.path(tab, "fig4c_chipseeker_peak_annotation_full.tsv.gz"),
  file.path(tab, "fig4c_chipseeker_annotation_counts.tsv"),
  file.path(fig, "fig4a_snatac_umap_label_transfer.png"),
  file.path(fig, "fig4a_snatac_umap_label_transfer.pdf"),
  file.path(fig, "fig4b_tissue_proportion_boxplots.png"),
  file.path(fig, "fig4b_tissue_proportion_boxplots.pdf"),
  file.path(fig, "fig4c_chipseeker_annotation_pie.png"),
  file.path(fig, "fig4c_chipseeker_annotation_pie.pdf"),
  file.path(files, "fig4_gse306459_snatac_processed.rds.not_saved.txt")
)
required_info <- file.info(required)
missing <- required[!file.exists(required) | required_info$size <= 0]
if (length(missing)) {
  stop("Missing or empty required files: ", paste(missing, collapse = "; "))
}

log_dt <- fread(file.path(logs, "fig4_gse306459_snatac.log"), sep = "\t", header = FALSE, fill = TRUE)
setnames(log_dt, c("timestamp", "step", "status", "note")[seq_len(ncol(log_dt))])
umap <- fread(file.path(tab, "fig4a_snatac_umap_labels.tsv.gz"))
da <- fread(file.path(tab, "fig4_differential_peaks_by_label.tsv.gz"))
anno <- fread(file.path(tab, "fig4c_chipseeker_annotation_counts.tsv"))
manifest <- fread(file.path(tab, "fig4_gse306459_geo_download_manifest.tsv"))
prop <- fread(file.path(tab, "fig4b_label_tissue_proportions.tsv"))
wilcox <- fread(file.path(tab, "fig4b_wilcoxon_by_label.tsv"))

expected_cols <- list(
  umap = c("cell_id", "UMAP_1", "UMAP_2", "raw_barcode", "Library", "Tissue", "Sample", "tissue_group", "seurat_clusters", "predicted_label", "prediction.score.max"),
  da = c("peak", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "predicted_label"),
  anno = c("annotation", "N", "fraction")
)
for (nm in names(expected_cols)) {
  observed <- names(get(nm))
  missing_cols <- setdiff(expected_cols[[nm]], observed)
  if (length(missing_cols)) {
    stop(nm, " missing columns: ", paste(missing_cols, collapse = ", "), "; observed: ", paste(observed, collapse = ", "))
  }
}

if (nrow(umap) != 122218L) stop("UMAP row count mismatch: ", nrow(umap))
if (nrow(da) <= 0L) stop("No differential peak rows")
if (nrow(anno) <= 0L) stop("No ChIPSeeker annotation rows")
frac_sum <- sum(anno$fraction)
if (abs(frac_sum - 1) > 1e-6) stop("Annotation fractions do not sum to 1: ", frac_sum)

merge_note <- log_dt[step == "merge" & status == "complete", tail(note, 1)]
peaks <- sub(".*peaks=([0-9]+).*", "\\1", merge_note)
transfer_note <- log_dt[step == "label_transfer" & status == "running", tail(note, 1)]
transfer_features <- sub("features=([0-9]+);.*", "\\1", transfer_note)

summary <- data.table(
  metric = c(
    "run_scope",
    "cells",
    "peaks",
    "libraries",
    "samples",
    "transferred_labels",
    "transfer_features",
    "differential_peak_rows",
    "unique_differential_peaks",
    "chipseeker_annotation_rows",
    "requested_motifs_found",
    "processed_rds"
  ),
  value = c(
    "fig4c",
    nrow(umap),
    peaks,
    uniqueN(umap$Library),
    uniqueN(umap$Sample),
    uniqueN(umap$predicted_label),
    transfer_features,
    nrow(da),
    uniqueN(da$peak),
    nrow(anno),
    "NA",
    "not_saved"
  )
)
fwrite(summary, file.path(tab, "fig4_gse306459_run_summary.tsv"), sep = "\t", quote = FALSE)

label_counts <- umap[, .N, by = predicted_label][order(-N)]
fwrite(label_counts, file.path(tab, "fig4a_label_counts.tsv"), sep = "\t", quote = FALSE)
da_counts <- da[, .N, by = predicted_label][order(predicted_label)]
fwrite(da_counts, file.path(tab, "fig4c_differential_peak_rows_by_label.tsv"), sep = "\t", quote = FALSE)

checks <- data.table(
  check = c(
    "required_files_nonempty",
    "umap_rows",
    "umap_columns",
    "da_rows_positive",
    "annotation_rows_positive",
    "annotation_fraction_sum",
    "libraries",
    "samples",
    "labels",
    "geo_h5_files",
    "wilcoxon_rows",
    "proportion_rows",
    "processed_rds_policy"
  ),
  result = rep("PASS", 13),
  value = c(
    length(required),
    nrow(umap),
    paste(names(umap), collapse = ";"),
    nrow(da),
    nrow(anno),
    sprintf("%.12f", frac_sum),
    uniqueN(umap$Library),
    uniqueN(umap$Sample),
    uniqueN(umap$predicted_label),
    nrow(manifest[grepl("raw_peak_bc_matrix.h5$", file_name)]),
    nrow(wilcox),
    nrow(prop),
    "not_saved"
  )
)
fwrite(checks, file.path(tab, "fig4_postrun_checks.tsv"), sep = "\t", quote = FALSE)

cat("checks_passed\n")
print(summary)
print(checks[, .(check, result, value)])
