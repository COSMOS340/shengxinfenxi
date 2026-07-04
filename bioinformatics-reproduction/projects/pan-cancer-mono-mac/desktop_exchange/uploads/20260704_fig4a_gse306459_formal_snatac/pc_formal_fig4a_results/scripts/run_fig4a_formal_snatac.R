suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(scales)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg) == 0L) stop("Cannot resolve script path from commandArgs().")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)

repo_dir <- normalizePath(Sys.getenv("FIG4_REPO_DIR", unset = getwd()), mustWork = TRUE)
request_dir <- normalizePath(
  Sys.getenv(
    "FIG4_REQUEST_DIR",
    unset = file.path(repo_dir, "desktop_exchange", "requests", "20260704_fig4a_gse306459_formal_snatac")
  ),
  mustWork = TRUE
)
geo_filelist <- file.path(repo_dir, "desktop_exchange", "requests", "20260703_fig4_gse306459_snatac", "inputs", "gse306459_geo_filelist.tsv")
geo_dir <- Sys.getenv("FIG4_GEO_DIR", unset = file.path(repo_dir, "fig4_snatac", "raw", "GSE306459"))
supp_dir <- Sys.getenv("FIG4_SUPP_DIR", unset = file.path(repo_dir, "desktop_exchange", "uploads", "20260703_fig4_gse306459_snatac", "mac_provided_inputs", "supplementary"))
out_dir <- Sys.getenv("FIG4_OUTPUT_DIR", unset = file.path(repo_dir, "desktop_exchange", "uploads", "20260704_fig4a_gse306459_formal_snatac", "pc_formal_fig4a_results"))
reference_rds <- Sys.getenv("FIG4_REFERENCE_RDS", unset = file.path(repo_dir, "desktop_exchange", "uploads", "20260703_fig4_gse306459_snatac", "mac_provided_inputs", "files", "fig4_coad_read_reference_with_fig4_display_label.rds"))
reference_label_column <- Sys.getenv("FIG4_REFERENCE_LABEL_COLUMN", unset = "fig4_display_label")
npcs <- as.integer(Sys.getenv("FIG4_NPCS", unset = "30"))
cluster_resolution <- as.numeric(Sys.getenv("FIG4_RESOLUTION", unset = "0.4"))

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
file_dir <- file.path(out_dir, "files")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(log_dir, "fig4a_formal_snatac.log")
if (file.exists(log_path)) file.remove(log_path)
write_log <- function(step, status, note = "") {
  line <- paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), step, status, note, sep = "\t")
  cat(line, "\n", file = log_path, append = TRUE)
  message(step, " | ", status, " | ", note)
}

options(error = function() {
  msg <- trimws(geterrmessage())
  cat(paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "error", "failed", msg, sep = "\t"), "\n", file = log_path, append = TRUE)
  traceback(2, max.lines = 20)
  quit(save = "no", status = 1)
})

stop_missing_file <- function(path, label) {
  if (!file.exists(path) || file.size(path) == 0) stop(label, " missing or empty: ", path)
}

assert_columns <- function(data, expected, label) {
  missing <- setdiff(expected, names(data))
  if (length(missing) > 0L) {
    stop(label, " missing columns: ", paste(missing, collapse = ", "), ". Observed columns: ", paste(names(data), collapse = ", "))
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

index_fragment_if_needed <- function(fragment_path) {
  tbi_path <- paste0(fragment_path, ".tbi")
  csi_path <- paste0(fragment_path, ".csi")
  if (file.exists(tbi_path) || file.exists(csi_path)) {
    return("existing")
  }
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Rsamtools is required to index fragment files.")
  }
  write_log("index_fragment", "running", basename(fragment_path))
  Rsamtools::indexTabix(fragment_path, format = "bed")
  if (!file.exists(tbi_path) && !file.exists(csi_path)) {
    stop("Fragment index was not created for ", fragment_path)
  }
  "created"
}

read_snatac_library <- function(library_id, file_pairs, barcode_meta) {
  row <- file_pairs[Library == library_id]
  h5_path <- file.path(geo_dir, row$H5)
  fragment_path <- file.path(geo_dir, row$TSV)
  stop_missing_file(h5_path, paste0(library_id, " H5"))
  stop_missing_file(fragment_path, paste0(library_id, " fragments"))
  observed_fragment_size <- file.size(fragment_path)
  expected_fragment_size <- as.numeric(row$TSV_size_bytes)
  if (!is.na(expected_fragment_size) && observed_fragment_size != expected_fragment_size) {
    stop("Fragment size mismatch for ", library_id, ": observed ", observed_fragment_size, ", expected ", expected_fragment_size)
  }

  index_status <- index_fragment_if_needed(fragment_path)
  write_log("read_library", "running", paste0(library_id, "; fragment_index=", index_status))
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
    cell_id = cell_ids,
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

  write_log("read_library", "complete", paste0(library_id, "; cells=", ncol(counts), "; peaks=", nrow(counts)))
  list(counts = counts, meta = meta, fragment = fragment_obj, fragment_path = fragment_path, raw_barcodes = raw_barcodes, cell_ids = cell_ids)
}

count_fragments_for_library <- function(lib_data, library_id) {
  write_log("count_fragments", "running", library_id)
  counts <- CountFragments(
    fragments = lib_data$fragment_path,
    cells = lib_data$raw_barcodes,
    max_lines = NULL,
    verbose = TRUE
  )
  counts <- as.data.table(counts)
  write_log("count_fragments", "columns", paste(names(counts), collapse = ";"))
  if (!"CB" %in% names(counts)) stop("CountFragments output missing CB column for ", library_id, ". Observed columns: ", paste(names(counts), collapse = ", "))
  total_col <- intersect(c("frequency_count", "read_count", "reads_count", "count", "fragments", "passed_filters"), names(counts))
  if (length(total_col) == 0L) {
    numeric_cols <- names(counts)[vapply(counts, is.numeric, logical(1))]
    if (length(numeric_cols) == 0L) stop("CountFragments output has no numeric count columns for ", library_id)
    total_col <- numeric_cols[[length(numeric_cols)]]
  } else {
    total_col <- total_col[[1]]
  }
  counts[, cell_id := paste(library_id, CB, sep = "__")]
  counts[, total_fragments := as.numeric(get(total_col))]
  out <- counts[, .(cell_id, total_fragments)]
  write_log("count_fragments", "complete", paste0(library_id, "; cells=", nrow(out), "; total_col=", total_col))
  out
}

write_log("start", "running", "Formal Fig4A GSE306459 snATAC workflow")
stop_missing_file(geo_filelist, "GSE306459 GEO file list")
stop_missing_file(reference_rds, "Reference Seurat RDS")
table_s1 <- file.path(supp_dir, "cir-24-1255_table_s1_suppst1.xlsx")
stop_missing_file(table_s1, "Supplementary Table S1")

required_namespaces <- c("hdf5r", "ensembldb", "EnsDb.Hsapiens.v86", "Rsamtools")
missing_namespaces <- required_namespaces[!vapply(required_namespaces, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_namespaces) > 0L) {
  stop("Missing R/Bioconductor packages: ", paste(missing_namespaces, collapse = ", "))
}

geo_files <- fread(geo_filelist)
assert_columns(geo_files, c("record_type", "file_name", "size_bytes", "file_type", "url"), "GSE306459 file list")
geo_files <- geo_files[record_type == "File"]
geo_files[, Library := sub("^GSM[0-9]+_([^_]+)_(fragments.tsv.gz|raw_peak_bc_matrix.h5)$", "\\1", file_name)]
if (any(geo_files$Library == geo_files$file_name)) {
  stop("Could not parse Library from file names: ", paste(geo_files[Library == file_name, file_name], collapse = ", "))
}
file_pairs <- dcast(geo_files, Library ~ file_type, value.var = "file_name")
size_pairs <- dcast(geo_files, Library ~ file_type, value.var = "size_bytes")
setnames(size_pairs, c("H5", "TSV"), c("H5_size_bytes", "TSV_size_bytes"), skip_absent = TRUE)
file_pairs <- merge(file_pairs, size_pairs, by = "Library", all.x = TRUE)
assert_columns(file_pairs, c("Library", "H5", "TSV", "H5_size_bytes", "TSV_size_bytes"), "GSE306459 paired files")
setorder(file_pairs, Library)

barcode_meta <- as.data.table(read_xlsx(table_s1, sheet = "barcode", skip = 1))
assert_columns(barcode_meta, c("Cell barcode", "Library", "Tissue", "Sample"), "Table S1 barcode sheet")
setnames(barcode_meta, "Cell barcode", "Cell_barcode")
barcode_meta[, tissue_group := fifelse(grepl("N$", Tissue), "Adjacent normal", fifelse(grepl("T$", Tissue), "Tumor", NA_character_))]
if (any(is.na(barcode_meta$tissue_group))) {
  stop("Table S1 Tissue values not ending in N/T: ", paste(unique(barcode_meta[is.na(tissue_group), Tissue]), collapse = ", "))
}

fragment_missing <- file_pairs[!file.exists(file.path(geo_dir, TSV)) | file.size(file.path(geo_dir, TSV)) != as.numeric(TSV_size_bytes)]
if (nrow(fragment_missing) > 0L) {
  stop("Missing or incomplete fragment files: ", paste(fragment_missing$TSV, collapse = "; "))
}

snatac_list <- lapply(file_pairs$Library, read_snatac_library, file_pairs = file_pairs, barcode_meta = barcode_meta)
names(snatac_list) <- file_pairs$Library

write_log("merge_counts", "running", paste0("libraries=", length(snatac_list)))
counts_list <- lapply(snatac_list, `[[`, "counts")
combined_counts <- Reduce(SeuratObject:::RowMergeSparseMatrices, counts_list)
combined_meta <- do.call(rbind, lapply(snatac_list, `[[`, "meta"))
rownames(combined_meta) <- combined_meta$cell_id
combined_meta$cell_id <- NULL
fragment_list <- lapply(snatac_list, `[[`, "fragment")
write_log("merge_counts", "complete", paste0("cells=", ncol(combined_counts), "; peaks=", nrow(combined_counts)))
if (!identical(colnames(combined_counts), rownames(combined_meta))) {
  stop("Merged count cell names do not match metadata row names.")
}

fragment_counts <- rbindlist(Map(count_fragments_for_library, snatac_list, names(snatac_list)), use.names = TRUE)
combined_meta$total_fragments <- fragment_counts$total_fragments[match(rownames(combined_meta), fragment_counts$cell_id)]
if (any(is.na(combined_meta$total_fragments))) {
  stop("Missing total fragment counts for cells: ", sum(is.na(combined_meta$total_fragments)))
}

snatac_list <- NULL
counts_list <- NULL
gc(verbose = FALSE)

write_log("create_chromatin_assay", "running", paste0("cells=", ncol(combined_counts), "; peaks=", nrow(combined_counts)))
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
fragment_list <- NULL
gc(verbose = FALSE)
DefaultAssay(combined) <- "peaks"
write_log("create_chromatin_assay", "complete", paste0("cells=", ncol(combined), "; peaks=", nrow(combined[["peaks"]])))

write_log("annotation", "running", "EnsDb.Hsapiens.v86")
annotations <- Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
old_seqlevels <- seqlevels(annotations)
new_seqlevels <- ifelse(grepl("^chr", old_seqlevels), old_seqlevels, paste0("chr", old_seqlevels))
names(new_seqlevels) <- old_seqlevels
annotations <- GenomeInfoDb::renameSeqlevels(annotations, new_seqlevels)
standard_chromosomes <- paste0("chr", c(seq_len(22), "X", "Y", "M"))
annotations <- GenomeInfoDb::keepSeqlevels(annotations, intersect(standard_chromosomes, seqlevels(annotations)), pruning.mode = "coarse")
genome(annotations) <- "hg38"
Annotation(combined) <- annotations
write_log("annotation", "complete", paste0("features=", length(annotations)))

write_log("qc_metrics", "running", "FRiP")
combined <- FRiP(combined, assay = "peaks", total.fragments = "total_fragments", col.name = "pct_reads_in_peaks", verbose = TRUE)
combined$pct_reads_in_peaks <- combined$pct_reads_in_peaks * 100

write_log("qc_metrics", "running", "NucleosomeSignal")
combined <- NucleosomeSignal(combined, assay = "peaks", verbose = TRUE)

write_log("qc_metrics", "running", "TSSEnrichment")
combined <- TSSEnrichment(combined, assay = "peaks", fast = TRUE, process_n = 2000, verbose = TRUE)

combined$blacklist_ratio <- NA_real_
blacklist_note <- "not_available"

qc_dt <- as.data.table(combined@meta.data, keep.rownames = "cell_id")
qc_cols <- c("cell_id", "raw_barcode", "Library", "Tissue", "Sample", "tissue_group", "nCount_peaks", "nFeature_peaks", "total_fragments", "pct_reads_in_peaks", "nucleosome_signal", "TSS.enrichment", "blacklist_ratio")
assert_columns(qc_dt, qc_cols, "QC metadata")
fwrite(qc_dt[, ..qc_cols], file.path(tab_dir, "fig4a_qc_metrics_by_cell.tsv.gz"), sep = "\t", quote = FALSE)
write_log("qc_metrics", "complete", paste0("cells=", nrow(qc_dt)))

threshold_rows <- list(
  data.table(metric = "nCount_peaks", operator = ">=", value = as.numeric(quantile(qc_dt$nCount_peaks, 0.01, na.rm = TRUE)), rule = "q01 lower bound"),
  data.table(metric = "nCount_peaks", operator = "<=", value = as.numeric(quantile(qc_dt$nCount_peaks, 0.995, na.rm = TRUE)), rule = "q995 upper bound"),
  data.table(metric = "nFeature_peaks", operator = ">=", value = as.numeric(quantile(qc_dt$nFeature_peaks, 0.01, na.rm = TRUE)), rule = "q01 lower bound"),
  data.table(metric = "nFeature_peaks", operator = "<=", value = as.numeric(quantile(qc_dt$nFeature_peaks, 0.995, na.rm = TRUE)), rule = "q995 upper bound"),
  data.table(metric = "TSS.enrichment", operator = ">=", value = max(2, as.numeric(quantile(qc_dt$TSS.enrichment, 0.10, na.rm = TRUE))), rule = "max(2, q10)"),
  data.table(metric = "nucleosome_signal", operator = "<=", value = as.numeric(quantile(qc_dt$nucleosome_signal, 0.95, na.rm = TRUE)), rule = "q95 upper bound"),
  data.table(metric = "pct_reads_in_peaks", operator = ">=", value = max(5, as.numeric(quantile(qc_dt$pct_reads_in_peaks, 0.05, na.rm = TRUE))), rule = "max(5, q05)"),
  data.table(metric = "blacklist_ratio", operator = "not_applied", value = NA_real_, rule = blacklist_note)
)
thresholds <- rbindlist(threshold_rows, use.names = TRUE)
fwrite(thresholds, file.path(tab_dir, "fig4a_qc_thresholds.tsv"), sep = "\t", quote = FALSE)

nc_low <- thresholds[metric == "nCount_peaks" & operator == ">=", value][1]
nc_high <- thresholds[metric == "nCount_peaks" & operator == "<=", value][1]
nf_low <- thresholds[metric == "nFeature_peaks" & operator == ">=", value][1]
nf_high <- thresholds[metric == "nFeature_peaks" & operator == "<=", value][1]
tss_low <- thresholds[metric == "TSS.enrichment" & operator == ">=", value][1]
nuc_high <- thresholds[metric == "nucleosome_signal" & operator == "<=", value][1]
frip_low <- thresholds[metric == "pct_reads_in_peaks" & operator == ">=", value][1]

qc_dt[, pass_qc := nCount_peaks >= nc_low &
  nCount_peaks <= nc_high &
  nFeature_peaks >= nf_low &
  nFeature_peaks <= nf_high &
  TSS.enrichment >= tss_low &
  nucleosome_signal <= nuc_high &
  pct_reads_in_peaks >= frip_low]

before_after <- qc_dt[, .(before_qc = .N, after_qc = sum(pass_qc)), by = Library][order(Library)]
before_after[, retention_fraction := after_qc / before_qc]
fwrite(before_after, file.path(tab_dir, "fig4a_cells_before_after_qc_by_library.tsv"), sep = "\t", quote = FALSE)
write_log("qc_filter", "complete", paste0("before=", nrow(qc_dt), "; after=", sum(qc_dt$pass_qc)))
if (sum(qc_dt$pass_qc) < 1000L) stop("QC filtering retained fewer than 1000 cells.")

qc_plot_dt <- melt(
  qc_dt[sample.int(nrow(qc_dt), min(nrow(qc_dt), 50000L))],
  id.vars = c("cell_id", "pass_qc"),
  measure.vars = c("nCount_peaks", "nFeature_peaks", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks"),
  variable.name = "metric",
  value.name = "value"
)
p_qc <- ggplot(qc_plot_dt, aes(pass_qc, value, fill = pass_qc)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.size = 0.1, linewidth = 0.2) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(title = "Fig. 4A formal snATAC QC metrics", x = "Pass QC", y = NULL, fill = "Pass QC") +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 10), legend.position = "top")
save_plot(p_qc, "fig4a_formal_qc_violin", width = 9.5, height = 3.4)

combined <- subset(combined, cells = qc_dt[pass_qc == TRUE, cell_id])
DefaultAssay(combined) <- "peaks"
write_log("lsi_umap_cluster", "running", paste0("npcs=", npcs))
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q0")
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:npcs)
combined <- FindNeighbors(combined, reduction = "lsi", dims = 2:npcs)
combined <- FindClusters(combined, resolution = cluster_resolution)
write_log("lsi_umap_cluster", "complete", paste0("cells=", ncol(combined), "; clusters=", length(unique(combined$seurat_clusters))))

reference <- readRDS(reference_rds)
ref_meta <- reference@meta.data
if (!reference_label_column %in% names(ref_meta)) {
  stop("Reference label column not found: ", reference_label_column, ". Observed metadata columns: ", paste(names(ref_meta), collapse = ", "))
}
ref_labels <- ref_meta[[reference_label_column]]
if (any(is.na(ref_labels))) stop("Reference label column contains NA values: ", reference_label_column)
reference_assay <- DefaultAssay(reference)
query_assay <- "RNA"
reference_features <- VariableFeatures(reference)
if (length(reference_features) < 500L) {
  reference_features <- rownames(get_assay_data(reference, reference_assay, "data"))
}
annotation_gene_names <- unique(Annotation(combined[["peaks"]])$gene_name)
annotation_gene_names <- annotation_gene_names[!is.na(annotation_gene_names) & nzchar(annotation_gene_names)]
transfer_features <- intersect(reference_features, annotation_gene_names)
if (length(transfer_features) < 500L) stop("Too few transfer features shared between reference and query: ", length(transfer_features))
if (length(transfer_features) > 3000L) transfer_features <- transfer_features[seq_len(3000L)]

write_log("gene_activity", "running", paste0("GeneActivity; features=", length(transfer_features)))
gene_activity <- GeneActivity(combined, features = transfer_features, process_n = 500)
write_log("gene_activity", "complete", paste0("genes=", nrow(gene_activity), "; cells=", ncol(gene_activity)))
combined[["RNA"]] <- CreateAssayObject(counts = gene_activity)
combined <- NormalizeData(combined, assay = "RNA")
query_features <- rownames(get_assay_data(combined, query_assay, "data"))
transfer_features <- intersect(transfer_features, query_features)
if (length(transfer_features) < 500L) stop("Too few transfer features retained after GeneActivity: ", length(transfer_features))

write_log("label_transfer", "running", paste0("features=", length(transfer_features), "; reduction=pcaproject"))
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = combined,
  reference.assay = reference_assay,
  query.assay = query_assay,
  features = transfer_features,
  reduction = "pcaproject",
  reference.reduction = "pca",
  dims = 1:min(30L, npcs)
)
predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = ref_labels,
  weight.reduction = combined[["lsi"]],
  dims = 2:npcs
)
combined <- AddMetaData(combined, predictions)
combined$predicted_label <- combined$predicted.id
write_log("label_transfer", "complete", paste0("labels=", paste(sort(unique(combined$predicted_label)), collapse = ";")))

umap_dt <- as.data.table(Embeddings(combined, "umap"), keep.rownames = "cell_id")
umap_columns <- setdiff(names(umap_dt), "cell_id")
if (length(umap_columns) < 2L) stop("UMAP embedding table has fewer than two coordinate columns.")
setnames(umap_dt, umap_columns[seq_len(2L)], c("UMAP_1", "UMAP_2"))
meta_dt <- as.data.table(combined@meta.data, keep.rownames = "cell_id")
umap_dt <- merge(
  umap_dt,
  meta_dt[, .(cell_id, raw_barcode, Library, Tissue, Sample, tissue_group, seurat_clusters, predicted_label, prediction.score.max, nCount_peaks, nFeature_peaks, total_fragments, pct_reads_in_peaks, nucleosome_signal, TSS.enrichment, blacklist_ratio)],
  by = "cell_id",
  all.x = TRUE,
  all.y = FALSE
)
fwrite(umap_dt, file.path(tab_dir, "fig4a_snatac_umap_labels.tsv.gz"), sep = "\t", quote = FALSE)

label_counts <- umap_dt[, .N, by = predicted_label][order(-N)]
fwrite(label_counts, file.path(tab_dir, "fig4a_label_counts.tsv"), sep = "\t", quote = FALSE)
score_summary <- umap_dt[, .(
  cells = .N,
  min_score = min(prediction.score.max, na.rm = TRUE),
  q25_score = as.numeric(quantile(prediction.score.max, 0.25, na.rm = TRUE)),
  median_score = median(prediction.score.max, na.rm = TRUE),
  mean_score = mean(prediction.score.max, na.rm = TRUE),
  q75_score = as.numeric(quantile(prediction.score.max, 0.75, na.rm = TRUE)),
  max_score = max(prediction.score.max, na.rm = TRUE)
), by = predicted_label][order(-cells)]
fwrite(score_summary, file.path(tab_dir, "fig4a_prediction_score_summary_by_label.tsv"), sep = "\t", quote = FALSE)
composition <- umap_dt[, .N, by = .(Library, predicted_label)][order(Library, -N)]
composition[, library_total := sum(N), by = Library]
composition[, fraction_in_library := N / library_total]
fwrite(composition, file.path(tab_dir, "fig4a_library_label_composition.tsv"), sep = "\t", quote = FALSE)

plot_point_size <- if (nrow(umap_dt) > 100000L) 0.16 else 0.25
p_label <- ggplot(umap_dt, aes(UMAP_1, UMAP_2, color = predicted_label)) +
  geom_point(size = plot_point_size, alpha = 0.75, stroke = 0) +
  coord_equal() +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
  labs(title = "Fig. 4A formal GSE306459 snATAC label transfer", x = "UMAP1", y = "UMAP2", color = NULL) +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10), legend.text = element_text(size = 7))
save_plot(p_label, "fig4a_formal_umap_label_transfer", width = 6.8, height = 5.6)

p_library <- ggplot(umap_dt, aes(UMAP_1, UMAP_2, color = Library)) +
  geom_point(size = plot_point_size, alpha = 0.7, stroke = 0) +
  coord_equal() +
  labs(title = "Fig. 4A formal UMAP by library", x = "UMAP1", y = "UMAP2", color = "Library") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10), legend.text = element_text(size = 6))
save_plot(p_library, "fig4a_formal_umap_by_library", width = 7.4, height = 5.8)

p_tissue <- ggplot(umap_dt, aes(UMAP_1, UMAP_2, color = tissue_group)) +
  geom_point(size = plot_point_size, alpha = 0.75, stroke = 0) +
  coord_equal() +
  scale_color_manual(values = c("Adjacent normal" = "#E97872", "Tumor" = "#23A7A6")) +
  labs(title = "Fig. 4A formal UMAP by tissue group", x = "UMAP1", y = "UMAP2", color = NULL) +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10))
save_plot(p_tissue, "fig4a_formal_umap_by_tissue_group", width = 6.4, height = 5.4)

p_score <- ggplot(umap_dt, aes(UMAP_1, UMAP_2, color = prediction.score.max)) +
  geom_point(size = plot_point_size, alpha = 0.8, stroke = 0) +
  coord_equal() +
  scale_color_viridis_c(option = "magma") +
  labs(title = "Fig. 4A formal UMAP prediction score", x = "UMAP1", y = "UMAP2", color = "Score") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10))
save_plot(p_score, "fig4a_formal_umap_prediction_score", width = 6.4, height = 5.4)

run_summary <- data.table(
  metric = c(
    "run_scope",
    "cells_before_qc",
    "cells_after_qc",
    "peaks_after_assay_creation",
    "libraries_before_qc",
    "libraries_after_qc",
    "samples_after_qc",
    "clusters_after_qc",
    "transferred_labels",
    "transfer_features",
    "blacklist_ratio"
  ),
  value = c(
    "fig4a_formal",
    nrow(qc_dt),
    ncol(combined),
    nrow(combined[["peaks"]]),
    uniqueN(qc_dt$Library),
    uniqueN(umap_dt$Library),
    uniqueN(umap_dt$Sample),
    length(unique(combined$seurat_clusters)),
    uniqueN(umap_dt$predicted_label),
    length(transfer_features),
    blacklist_note
  )
)
fwrite(run_summary, file.path(tab_dir, "fig4a_formal_run_summary.tsv"), sep = "\t", quote = FALSE)
writeLines("Full processed Seurat object not uploaded; object is expected to be too large for GitHub handoff. Re-run script locally with the same inputs to recreate it.", file.path(file_dir, "fig4a_formal_processed_object_not_uploaded.txt"))

write_log("complete", "complete", paste0("output_dir=", out_dir))
