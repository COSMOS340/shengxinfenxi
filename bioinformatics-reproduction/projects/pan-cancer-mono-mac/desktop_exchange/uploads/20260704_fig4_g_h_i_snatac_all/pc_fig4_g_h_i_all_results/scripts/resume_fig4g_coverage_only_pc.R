#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

repo <- getwd()
request_dir <- Sys.getenv("FIG4_REQUEST_DIR", unset = file.path(repo, "desktop_exchange", "requests", "20260704_fig4_g_h_i_snatac_all"))
geo_filelist <- file.path(request_dir, "inputs", "gse306459_geo_filelist.tsv")
geo_dir <- Sys.getenv("FIG4_GEO_DIR", unset = file.path(repo, "fig4_snatac", "raw", "GSE306459"))
supp_dir <- Sys.getenv("FIG4_SUPP_DIR", unset = file.path(repo, "desktop_exchange", "uploads", "20260703_fig4_gse306459_snatac", "mac_provided_inputs", "supplementary"))
out_dir <- Sys.getenv("FIG4_OUTPUT_DIR", unset = file.path(repo, "fig4_snatac", "outputs_fig4g_hi_all"))
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")
file_dir <- file.path(out_dir, "files")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(log_dir, "fig4_gse306459_snatac.log")

write_log <- function(step, status, note = "") {
  line <- paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), step, status, note, sep = "\t")
  cat(line, "\n", file = log_path, append = TRUE)
  message(step, " | ", status, " | ", note)
}

stop_missing_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
}

assert_columns <- function(data, expected, label) {
  missing <- setdiff(expected, names(data))
  if (length(missing)) {
    stop(label, " missing columns: ", paste(missing, collapse = ", "),
         ". Observed columns: ", paste(names(data), collapse = ", "))
  }
}

get_assay_data <- function(object, assay, layer_name) {
  tryCatch(
    GetAssayData(object, assay = assay, layer = layer_name),
    error = function(e) GetAssayData(object, assay = assay, slot = layer_name)
  )
}

save_plot <- function(plot, stem, width, height, dpi = 300) {
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, bg = "white", limitsize = FALSE)
  write_log("save_plot", "complete", png_path)
}

write_log("fig4g_resume", "running", "coverage_only")
stop_missing_file(geo_filelist, "GSE306459 GEO file list")
table_s1 <- file.path(supp_dir, "cir-24-1255_table_s1_suppst1.xlsx")
stop_missing_file(table_s1, "Supplementary Table S1")

geo_files <- fread(geo_filelist)
assert_columns(geo_files, c("record_type", "file_name", "size_bytes", "file_type", "url"), "GSE306459 file list")
geo_files <- geo_files[record_type == "File"]
geo_files[, Library := sub("^GSM[0-9]+_([^_]+)_(fragments.tsv.gz|raw_peak_bc_matrix.h5)$", "\\1", file_name)]
if (any(geo_files$Library == geo_files$file_name)) {
  stop("Could not parse Library from file names: ", paste(geo_files[Library == file_name, file_name], collapse = ", "))
}
file_pairs <- dcast(geo_files, Library ~ file_type, value.var = "file_name")
assert_columns(file_pairs, c("Library", "H5", "TSV"), "GSE306459 paired files")
setorder(file_pairs, Library)

barcode_meta <- as.data.table(read_xlsx(table_s1, sheet = "barcode", skip = 1))
assert_columns(barcode_meta, c("Cell barcode", "Library", "Tissue", "Sample"), "Table S1 barcode sheet")
setnames(barcode_meta, "Cell barcode", "Cell_barcode")
barcode_meta[, tissue_group := fifelse(grepl("N$", Tissue), "Adjacent normal", fifelse(grepl("T$", Tissue), "Tumor", NA_character_))]

read_snatac_library <- function(library_id) {
  row <- file_pairs[Library == library_id]
  h5_path <- file.path(geo_dir, row$H5)
  fragment_path <- file.path(geo_dir, row$TSV)
  stop_missing_file(h5_path, paste0(library_id, " H5"))
  stop_missing_file(fragment_path, paste0(library_id, " fragments"))
  write_log("fig4g_read_library", "running", library_id)
  counts <- Read10X_h5(h5_path)
  if (is.list(counts)) counts <- counts[[1]]
  cells_keep <- intersect(colnames(counts), barcode_meta[Library == library_id, Cell_barcode])
  if (length(cells_keep) == 0L) stop("No Table S1 barcodes found in H5 for ", library_id)
  counts <- counts[, cells_keep, drop = FALSE]
  raw_barcodes <- colnames(counts)
  cell_ids <- paste(library_id, raw_barcodes, sep = "__")
  colnames(counts) <- cell_ids

  lib_meta <- barcode_meta[Library == library_id]
  lib_meta <- lib_meta[match(raw_barcodes, Cell_barcode)]
  if (any(is.na(lib_meta$Cell_barcode))) stop("Metadata match failed for ", library_id)
  meta <- data.frame(
    raw_barcode = raw_barcodes,
    Library = lib_meta$Library,
    Tissue = lib_meta$Tissue,
    Sample = lib_meta$Sample,
    tissue_group = lib_meta$tissue_group,
    row.names = cell_ids,
    check.names = FALSE
  )
  fragment_cells <- raw_barcodes
  names(fragment_cells) <- cell_ids
  fragment_obj <- CreateFragmentObject(
    path = fragment_path,
    cells = fragment_cells,
    validate.fragments = FALSE,
    verbose = FALSE
  )
  write_log("fig4g_read_library", "complete", paste0(library_id, "; cells=", ncol(counts), "; peaks=", nrow(counts)))
  list(counts = counts, meta = meta, fragment = fragment_obj)
}

snatac_list <- lapply(file_pairs$Library, read_snatac_library)
names(snatac_list) <- file_pairs$Library
write_log("fig4g_merge_counts", "running", paste0("libraries=", length(snatac_list)))
counts_list <- lapply(snatac_list, `[[`, "counts")
fragment_list <- lapply(snatac_list, `[[`, "fragment")
combined_counts <- Reduce(SeuratObject:::RowMergeSparseMatrices, counts_list)
combined_meta <- do.call(rbind, lapply(snatac_list, `[[`, "meta"))
snatac_list <- NULL
counts_list <- NULL
gc(verbose = FALSE)

combined_assay <- CreateChromatinAssay(
  counts = combined_counts,
  sep = c(":", "-"),
  min.cells = 1,
  min.features = 1,
  fragments = fragment_list,
  validate.fragments = FALSE
)
combined <- CreateSeuratObject(counts = combined_assay, assay = "peaks", project = "GSE306459", meta.data = combined_meta)
combined_counts <- NULL
combined_assay <- NULL
combined_meta <- NULL
gc(verbose = FALSE)
DefaultAssay(combined) <- "peaks"
write_log("fig4g_merge", "complete", paste0("cells=", ncol(combined), "; peaks=", nrow(combined)))

write_log("fig4g_annotation", "running", "EnsDb.Hsapiens.v86")
annotations <- Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
old_seqlevels <- seqlevels(annotations)
new_seqlevels <- ifelse(grepl("^chr", old_seqlevels), old_seqlevels, paste0("chr", old_seqlevels))
names(new_seqlevels) <- old_seqlevels
annotations <- GenomeInfoDb::renameSeqlevels(annotations, new_seqlevels)
standard_chromosomes <- paste0("chr", c(seq_len(22), "X", "Y", "M"))
annotations <- GenomeInfoDb::keepSeqlevels(annotations, intersect(standard_chromosomes, seqlevels(annotations)), pruning.mode = "coarse")
genome(annotations) <- "hg38"
Annotation(combined) <- annotations
write_log("fig4g_annotation", "complete", paste0("features=", length(annotations)))

umap_labels <- fread(file.path(tab_dir, "fig4a_snatac_umap_labels.tsv.gz"))
assert_columns(umap_labels, c("cell_id", "predicted_label"), "Fig4A UMAP label table")
combined$predicted_label <- umap_labels$predicted_label[match(colnames(combined), umap_labels$cell_id)]
if (any(is.na(combined$predicted_label))) {
  stop("Missing predicted_label for cells: ", sum(is.na(combined$predicted_label)))
}

track_genes <- c("C1QB", "CHID1", "IL1B", "SLC25A37")
track_summary <- data.table(gene = track_genes, status = "requested")
fwrite(track_summary, file.path(tab_dir, "fig4g_track_region_summary.tsv"), sep = "\t", quote = FALSE)
track_plots <- vector("list", length(track_genes))
for (i in seq_along(track_genes)) {
  gene <- track_genes[[i]]
  write_log("fig4g_coverage_plot", "running", gene)
  track_plots[[i]] <- CoveragePlot(
    object = combined,
    region = gene,
    group.by = "predicted_label",
    extend.upstream = 8000,
    extend.downstream = 8000,
    peaks = TRUE,
    annotation = TRUE
  ) + ggtitle(gene)
  write_log("fig4g_coverage_plot", "complete", gene)
  gc(verbose = FALSE)
}
save_plot(wrap_plots(track_plots, ncol = 2), "fig4g_gene_coverage_tracks", width = 10.5, height = 7.5)

log_lines <- readLines(log_path, warn = FALSE)
transfer_line <- tail(grep("label_transfer\\trunning\\tfeatures=", log_lines, value = TRUE), 1)
transfer_features <- sub(".*features=([0-9]+).*", "\\1", transfer_line)
if (!grepl("^[0-9]+$", transfer_features)) transfer_features <- NA_character_
da_rows <- nrow(fread(file.path(tab_dir, "fig4_differential_peaks_by_label.tsv.gz")))
requested_motifs_found <- nrow(fread(file.path(tab_dir, "fig4_requested_motif_ids.tsv")))
run_summary <- data.table(
  metric = c(
    "run_scope",
    "cells",
    "peaks",
    "libraries",
    "samples",
    "transferred_labels",
    "transfer_features",
    "differential_peak_rows",
    "requested_motifs_found"
  ),
  value = c(
    "all",
    ncol(combined),
    nrow(combined[["peaks"]]),
    length(unique(combined$Library)),
    length(unique(combined$Sample)),
    length(unique(combined$predicted_label)),
    transfer_features,
    da_rows,
    requested_motifs_found
  )
)
fwrite(run_summary, file.path(tab_dir, "fig4_gse306459_run_summary.tsv"), sep = "\t", quote = FALSE)
writeLines(
  "processed RDS not saved; set FIG4_SAVE_RDS=true to enable this large output",
  file.path(file_dir, "fig4_gse306459_snatac_processed.rds.not_saved.txt")
)
write_log("complete", "complete", paste0("output_dir=", out_dir))
