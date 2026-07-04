#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg) == 0L) stop("Cannot resolve script path from commandArgs().")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
request_dir_env <- Sys.getenv("FIG4_REQUEST_DIR", unset = "")
request_dir <- if (nzchar(request_dir_env)) {
  normalizePath(request_dir_env, mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}

geo_filelist <- file.path(request_dir, "inputs", "gse306459_geo_filelist.tsv")
geo_dir <- Sys.getenv("FIG4_GEO_DIR", unset = file.path(getwd(), "fig4_snatac", "raw", "GSE306459"))
supp_dir <- Sys.getenv("FIG4_SUPP_DIR", unset = file.path(getwd(), "fig4_snatac", "inputs", "supplementary"))
out_dir <- Sys.getenv("FIG4_OUTPUT_DIR", unset = file.path(getwd(), "fig4_snatac", "outputs"))
reference_rds <- Sys.getenv("FIG4_REFERENCE_RDS", unset = "")
reference_label_column <- Sys.getenv("FIG4_REFERENCE_LABEL_COLUMN", unset = "")
run_scope <- Sys.getenv("FIG4_RUN_SCOPE", unset = "fig4c")
npcs <- as.integer(Sys.getenv("FIG4_NPCS", unset = "30"))
cluster_resolution <- as.numeric(Sys.getenv("FIG4_RESOLUTION", unset = "0.4"))
if (!run_scope %in% c("fig4c", "all")) stop("FIG4_RUN_SCOPE must be 'fig4c' or 'all'.")

required_namespaces <- c(
  "hdf5r",
  "ensembldb",
  "EnsDb.Hsapiens.v86",
  "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db"
)
if (run_scope == "all") {
  required_namespaces <- c(
    required_namespaces,
    "BSgenome.Hsapiens.UCSC.hg38",
    "JASPAR2020",
    "TFBSTools",
    "ggseqlogo"
  )
}
missing_namespaces <- required_namespaces[!vapply(required_namespaces, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_namespaces) > 0L) {
  stop("Missing R/Bioconductor packages: ", paste(missing_namespaces, collapse = ", "))
}

table_s1 <- file.path(supp_dir, "cir-24-1255_table_s1_suppst1.xlsx")
table_s4 <- file.path(supp_dir, "cir-24-1255_table_s4_suppst4.xlsx")

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
file_dir <- file.path(out_dir, "files")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(log_dir, "fig4_gse306459_snatac.log")

write_log <- function(step, status, note = "") {
  line <- paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), step, status, note, sep = "\t")
  cat(line, "\n", file = log_path, append = TRUE)
  message(step, " | ", status, " | ", note)
}

options(error = function() {
  msg <- geterrmessage()
  line <- paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "error", "failed", trimws(msg), sep = "\t")
  cat(line, "\n", file = log_path, append = TRUE)
  traceback(2, max.lines = 20)
  quit(save = "no", status = 1)
})

stop_missing_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
}

assert_columns <- function(data, expected, label) {
  missing <- setdiff(expected, names(data))
  if (length(missing) > 0L) {
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

compute_gene_activity_from_peaks <- function(object, annotation_gr) {
  peak_counts <- get_assay_data(object, "peaks", "counts")
  peak_gr <- StringToGRanges(rownames(peak_counts), sep = c(":", "-"))
  genome(peak_gr) <- "hg38"

  gene_names <- as.character(mcols(annotation_gr)$gene_name)
  valid <- !is.na(gene_names) & nzchar(gene_names)
  annotation_gr <- annotation_gr[valid]
  gene_names <- gene_names[valid]

  gene_grl <- reduce(split(annotation_gr, gene_names), ignore.strand = TRUE)
  gene_ranges <- unlist(gene_grl, use.names = FALSE)
  range_gene <- rep(names(gene_grl), lengths(gene_grl))
  overlaps <- findOverlaps(peak_gr, gene_ranges, ignore.strand = TRUE)
  if (length(overlaps) == 0L) {
    stop("No peak-gene overlaps found for fragments-free gene activity.")
  }

  overlap_dt <- unique(data.table(
    peak_index = queryHits(overlaps),
    gene_name = range_gene[subjectHits(overlaps)]
  ))
  gene_levels <- sort(unique(overlap_dt$gene_name))
  gene_peak <- sparseMatrix(
    i = match(overlap_dt$gene_name, gene_levels),
    j = overlap_dt$peak_index,
    x = 1,
    dims = c(length(gene_levels), nrow(peak_counts)),
    dimnames = list(gene_levels, rownames(peak_counts))
  )
  gene_activity <- gene_peak %*% peak_counts
  gene_activity <- as(gene_activity, "dgCMatrix")
  gene_activity[rowSums(gene_activity) > 0, , drop = FALSE]
}

save_plot <- function(plot, stem, width, height, dpi = 300) {
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, bg = "white", limitsize = FALSE)
  write_log("save_plot", "complete", png_path)
  invisible(png_path)
}

finalize_run <- function(requested_motifs_found = NA_integer_) {
  processed_rds <- file.path(file_dir, "fig4_gse306459_snatac_processed.rds")
  save_processed_rds <- identical(tolower(Sys.getenv("FIG4_SAVE_RDS", unset = "false")), "true")
  if (save_processed_rds) {
    saveRDS(combined, processed_rds)
    rds_sha <- tools::sha256sum(processed_rds)
    writeLines(paste(unname(rds_sha), basename(processed_rds)), file.path(file_dir, "fig4_gse306459_snatac_processed.rds.sha256"))
  } else {
    writeLines(
      "processed RDS not saved; set FIG4_SAVE_RDS=true to enable this large output",
      file.path(file_dir, "fig4_gse306459_snatac_processed.rds.not_saved.txt")
    )
  }

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
      run_scope,
      ncol(combined),
      nrow(combined[["peaks"]]),
      length(unique(combined$Library)),
      length(unique(combined$Sample)),
      length(unique(combined$predicted_label)),
      length(transfer_features),
      if (exists("da_dt")) nrow(da_dt) else NA_integer_,
      requested_motifs_found
    )
  )
  fwrite(run_summary, file.path(tab_dir, "fig4_gse306459_run_summary.tsv"), sep = "\t", quote = FALSE)
  write_log("complete", "complete", paste0("output_dir=", out_dir))
}

write_log("start", "running", "Fig. 4 GSE306459 snATAC workflow")
stop_missing_file(geo_filelist, "GSE306459 GEO file list")
stop_missing_file(table_s1, "Supplementary Table S1")
stop_missing_file(table_s4, "Supplementary Table S4")
if (!nzchar(reference_rds)) stop("FIG4_REFERENCE_RDS is required.")
if (!nzchar(reference_label_column)) stop("FIG4_REFERENCE_LABEL_COLUMN is required.")
stop_missing_file(reference_rds, "Reference Seurat RDS")

geo_files <- fread(geo_filelist)
assert_columns(geo_files, c("record_type", "file_name", "size_bytes", "file_type", "url"), "GSE306459 file list")
archive_row <- geo_files[record_type == "Archive" & file_name == "GSE306459_RAW.tar"]
if (nrow(archive_row) != 1L) stop("Expected exactly one GSE306459_RAW.tar archive row in GEO file list.")
raw_archive <- file.path(geo_dir, archive_row$file_name)
geo_files <- geo_files[record_type == "File"]
geo_files[, url := paste0("archive_member:", archive_row$file_name, ":", file_name)]
geo_files[, Library := sub("^GSM[0-9]+_([^_]+)_(fragments.tsv.gz|raw_peak_bc_matrix.h5)$", "\\1", file_name)]
if (any(geo_files$Library == geo_files$file_name)) {
  stop("Could not parse Library from file names: ", paste(geo_files[Library == file_name, file_name], collapse = ", "))
}

needed_file_type <- if (run_scope == "all") c("H5", "TSV") else "H5"
needed_files <- geo_files[file_type %in% needed_file_type]
download_manifest <- rbind(
  archive_row[, .(record_type, file_name, url, size_bytes, local_path = raw_archive, required_for_scope = "download_archive")],
  needed_files[, .(record_type, file_name, url, size_bytes, local_path = file.path(geo_dir, file_name), required_for_scope = run_scope)],
  use.names = TRUE,
  fill = TRUE
)

missing_needed <- needed_files[!file.exists(file.path(geo_dir, file_name)) | file.size(file.path(geo_dir, file_name)) == 0]
if (nrow(missing_needed) > 0L) {
  if (file.exists(raw_archive) && file.size(raw_archive) > 0) {
    observed_archive_size <- file.size(raw_archive)
    if (!is.na(archive_row$size_bytes) && observed_archive_size != as.numeric(archive_row$size_bytes)) {
      stop("Existing archive size mismatch: observed ", observed_archive_size, ", expected ", archive_row$size_bytes)
    }
    write_log("download_archive", "skipped_existing", raw_archive)
  } else {
    write_log("download_archive", "running", archive_row$file_name)
    download.file(archive_row$url, raw_archive, mode = "wb", quiet = FALSE)
    observed_archive_size <- file.size(raw_archive)
    if (!is.na(archive_row$size_bytes) && observed_archive_size != as.numeric(archive_row$size_bytes)) {
      stop("Downloaded archive size mismatch: observed ", observed_archive_size, ", expected ", archive_row$size_bytes)
    }
    write_log("download_archive", "complete", paste0(archive_row$file_name, "; bytes=", observed_archive_size))
  }

  archive_listing <- utils::untar(raw_archive, list = TRUE)
  extract_names <- missing_needed$file_name
  archive_name_map <- setNames(archive_listing[match(extract_names, basename(archive_listing))], extract_names)
  if (any(is.na(archive_name_map))) {
    stop("Archive does not contain required files: ", paste(names(archive_name_map)[is.na(archive_name_map)], collapse = ", "))
  }
  write_log("extract_archive", "running", paste0("files=", length(archive_name_map)))
  utils::untar(raw_archive, files = unname(archive_name_map), exdir = geo_dir)
  for (internal_name in unname(archive_name_map)) {
    extracted_path <- file.path(geo_dir, internal_name)
    final_path <- file.path(geo_dir, basename(internal_name))
    if (!identical(normalizePath(extracted_path, mustWork = FALSE), normalizePath(final_path, mustWork = FALSE))) {
      file.rename(extracted_path, final_path)
    }
  }
  write_log("extract_archive", "complete", paste0("files=", length(archive_name_map)))
}

for (i in seq_len(nrow(needed_files))) {
  target_path <- file.path(geo_dir, needed_files$file_name[i])
  if (!file.exists(target_path) || file.size(target_path) == 0) stop("Required extracted file missing: ", target_path)
  observed_size <- file.size(target_path)
  if (!is.na(needed_files$size_bytes[i]) && observed_size != as.numeric(needed_files$size_bytes[i])) {
    stop("Extracted size mismatch for ", needed_files$file_name[i], ": observed ", observed_size, ", expected ", needed_files$size_bytes[i])
  }
}
fwrite(download_manifest, file.path(tab_dir, "fig4_gse306459_geo_download_manifest.tsv"), sep = "\t", quote = FALSE)

file_pairs <- dcast(geo_files, Library ~ file_type, value.var = "file_name")
assert_columns(file_pairs, c("Library", "H5"), "GSE306459 paired files")
if (run_scope == "all") assert_columns(file_pairs, c("Library", "H5", "TSV"), "GSE306459 paired files")
setorder(file_pairs, Library)

barcode_meta <- as.data.table(read_xlsx(table_s1, sheet = "barcode", skip = 1))
assert_columns(barcode_meta, c("Cell barcode", "Library", "Tissue", "Sample"), "Table S1 barcode sheet")
setnames(barcode_meta, "Cell barcode", "Cell_barcode")
barcode_meta[, tissue_group := fifelse(grepl("N$", Tissue), "Adjacent normal",
                                 fifelse(grepl("T$", Tissue), "Tumor", NA_character_))]
if (any(is.na(barcode_meta$tissue_group))) {
  stop("Table S1 Tissue values not ending in N/T: ", paste(unique(barcode_meta[is.na(tissue_group), Tissue]), collapse = ", "))
}

library_count <- barcode_meta[, .N, by = .(Library, tissue_group)][order(Library, tissue_group)]
fwrite(library_count, file.path(tab_dir, "fig4_table_s1_barcode_counts_by_library.tsv"), sep = "\t", quote = FALSE)

missing_libraries <- setdiff(unique(barcode_meta$Library), file_pairs$Library)
if (length(missing_libraries) > 0L) {
  stop("Libraries in Table S1 but not in GEO file list: ", paste(missing_libraries, collapse = ", "))
}

read_snatac_library <- function(library_id) {
  row <- file_pairs[Library == library_id]
  h5_path <- file.path(geo_dir, row$H5)
  stop_missing_file(h5_path, paste0(library_id, " H5"))
  fragment_path <- NULL
  if (run_scope == "all") {
    fragment_path <- file.path(geo_dir, row$TSV)
    stop_missing_file(fragment_path, paste0(library_id, " fragments"))
  }
  write_log("read_library", "running", library_id)

  counts <- Read10X_h5(h5_path)
  if (is.list(counts)) {
    counts <- counts[[1]]
  }
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
  write_log("read_library", "complete", paste0(library_id, "; cells=", ncol(counts), "; peaks=", nrow(counts)))
  list(counts = counts, meta = meta)
}

snatac_list <- lapply(file_pairs$Library, read_snatac_library)
names(snatac_list) <- file_pairs$Library
write_log("merge_counts", "running", paste0("libraries=", length(snatac_list)))
counts_list <- lapply(snatac_list, `[[`, "counts")
combined_counts <- Reduce(SeuratObject:::RowMergeSparseMatrices, counts_list)
write_log("merge_counts", "complete", paste0("cells=", ncol(combined_counts), "; peaks=", nrow(combined_counts)))
combined_meta <- do.call(rbind, lapply(snatac_list, `[[`, "meta"))
snatac_list <- NULL
counts_list <- NULL
gc(verbose = FALSE)
write_log("create_chromatin_assay", "running", paste0("cells=", ncol(combined_counts), "; peaks=", nrow(combined_counts)))
combined_assay <- CreateChromatinAssay(
  counts = combined_counts,
  sep = c(":", "-"),
  min.cells = 1,
  min.features = 1
)
combined <- CreateSeuratObject(
  counts = combined_assay,
  assay = "peaks",
  project = "GSE306459",
  meta.data = combined_meta
)
combined_counts <- NULL
combined_assay <- NULL
combined_meta <- NULL
gc(verbose = FALSE)
DefaultAssay(combined) <- "peaks"
write_log("merge", "complete", paste0("cells=", ncol(combined), "; peaks=", nrow(combined)))

write_log("annotation", "running", "EnsDb.Hsapiens.v86")
annotations <- Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
old_seqlevels <- seqlevels(annotations)
new_seqlevels <- ifelse(grepl("^chr", old_seqlevels), old_seqlevels, paste0("chr", old_seqlevels))
names(new_seqlevels) <- old_seqlevels
annotations <- GenomeInfoDb::renameSeqlevels(annotations, new_seqlevels)
standard_chromosomes <- paste0("chr", c(seq_len(22), "X", "Y", "M"))
annotations <- GenomeInfoDb::keepSeqlevels(
  annotations,
  intersect(standard_chromosomes, seqlevels(annotations)),
  pruning.mode = "coarse"
)
genome(annotations) <- "hg38"
Annotation(combined) <- annotations
write_log("annotation", "complete", paste0("features=", length(annotations)))

write_log("lsi_umap_cluster", "running", paste0("npcs=", npcs))
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q0")
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:npcs)
combined <- FindNeighbors(combined, reduction = "lsi", dims = 2:npcs)
combined <- FindClusters(combined, resolution = cluster_resolution)
write_log("lsi_umap_cluster", "complete", paste0("clusters=", length(unique(combined$seurat_clusters))))

write_log("gene_activity", "running", "fragments_free_peak_overlap")
if (run_scope == "fig4c") {
  gene_activity <- compute_gene_activity_from_peaks(combined, annotations)
} else {
  gene_activity <- GeneActivity(combined)
}
write_log("gene_activity", "complete", paste0("genes=", nrow(gene_activity), "; cells=", ncol(gene_activity)))
combined[["RNA"]] <- CreateAssayObject(counts = gene_activity)
combined <- NormalizeData(combined, assay = "RNA")

reference <- readRDS(reference_rds)
ref_meta <- reference@meta.data
if (!reference_label_column %in% names(ref_meta)) {
  stop("Reference label column not found: ", reference_label_column,
       ". Observed metadata columns: ", paste(names(ref_meta), collapse = ", "))
}
ref_labels <- ref_meta[[reference_label_column]]
if (any(is.na(ref_labels))) {
  stop("Reference label column contains NA values: ", reference_label_column)
}

reference_assay <- DefaultAssay(reference)
query_assay <- "RNA"
query_features <- rownames(get_assay_data(combined, query_assay, "data"))
reference_features <- VariableFeatures(reference)
if (length(reference_features) < 500L) {
  reference_features <- rownames(get_assay_data(reference, reference_assay, "data"))
}
transfer_features <- intersect(reference_features, query_features)
if (length(transfer_features) < 500L) {
  stop("Too few transfer features shared between reference and query: ", length(transfer_features))
}
if (length(transfer_features) > 3000L) {
  transfer_features <- transfer_features[seq_len(3000L)]
}
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
  meta_dt[, .(cell_id, raw_barcode, Library, Tissue, Sample, tissue_group, seurat_clusters, predicted_label, prediction.score.max)],
  by = "cell_id",
  all.x = TRUE,
  all.y = FALSE
)
fwrite(umap_dt, file.path(tab_dir, "fig4a_snatac_umap_labels.tsv.gz"), sep = "\t", quote = FALSE)

p_umap <- ggplot(umap_dt, aes(UMAP_1, UMAP_2, color = predicted_label)) +
  geom_point(size = 0.18, alpha = 0.75, stroke = 0) +
  coord_equal() +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
  labs(title = "Fig. 4A GSE306459 snATAC label transfer", x = "UMAP1", y = "UMAP2", color = NULL) +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10), legend.text = element_text(size = 7))
save_plot(p_umap, "fig4a_snatac_umap_label_transfer", width = 6.5, height = 5.5)

prop_dt <- umap_dt[, .N, by = .(Sample, Tissue, tissue_group, predicted_label)]
prop_dt[, sample_total := sum(N), by = .(Sample, Tissue, tissue_group)]
prop_dt[, proportion := N / sample_total]
fwrite(prop_dt, file.path(tab_dir, "fig4b_label_tissue_proportions.tsv"), sep = "\t", quote = FALSE)

wilcox_dt <- prop_dt[, {
  group_counts <- table(tissue_group)
  if (!all(c("Adjacent normal", "Tumor") %in% names(group_counts)) || any(group_counts[c("Adjacent normal", "Tumor")] < 2)) {
    .(p_value = NA_real_, n_adjacent_normal = unname(group_counts["Adjacent normal"]), n_tumor = unname(group_counts["Tumor"]))
  } else {
    test <- wilcox.test(proportion ~ tissue_group, data = .SD, exact = FALSE)
    .(p_value = test$p.value, n_adjacent_normal = unname(group_counts["Adjacent normal"]), n_tumor = unname(group_counts["Tumor"]))
  }
}, by = predicted_label]
wilcox_dt[, p_adj := p.adjust(p_value, method = "BH")]
fwrite(wilcox_dt, file.path(tab_dir, "fig4b_wilcoxon_by_label.tsv"), sep = "\t", quote = FALSE)

plot_labels <- unique(wilcox_dt[!is.na(p_value), predicted_label])
p_box <- ggplot(prop_dt[predicted_label %in% plot_labels], aes(tissue_group, proportion, fill = tissue_group)) +
  geom_boxplot(width = 0.65, outlier.size = 0.6, linewidth = 0.25) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.65) +
  facet_wrap(~ predicted_label, scales = "free_y", nrow = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Adjacent normal" = "#E97872", "Tumor" = "#23A7A6")) +
  labs(title = "Fig. 4B adjacent normal versus tumor proportions", x = NULL, y = "Proportion", fill = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "top")
save_plot(p_box, "fig4b_tissue_proportion_boxplots", width = 8.5, height = 5.2)

DefaultAssay(combined) <- "peaks"
Idents(combined) <- combined$predicted_label
peak_counts_for_da <- get_assay_data(combined, "peaks", "counts")
da_list <- lapply(sort(unique(combined$predicted_label)), function(label) {
  write_log("differential_peaks", "running", paste0(label, "; method=sparse_fold_change"))
  cells_1 <- colnames(combined)[combined$predicted_label == label]
  cells_2 <- colnames(combined)[combined$predicted_label != label]
  if (length(cells_1) == 0L || length(cells_2) == 0L) return(NULL)

  mat_1 <- peak_counts_for_da[, cells_1, drop = FALSE]
  mat_2 <- peak_counts_for_da[, cells_2, drop = FALSE]
  pct_1 <- Matrix::rowMeans(mat_1 > 0)
  pct_2 <- Matrix::rowMeans(mat_2 > 0)
  mean_1 <- Matrix::rowMeans(mat_1)
  mean_2 <- Matrix::rowMeans(mat_2)
  keep <- which(pct_1 >= 0.05 & (pct_1 - pct_2) > 0 & log2((mean_1 + 1e-6) / (mean_2 + 1e-6)) > 0.25)
  if (length(keep) == 0L) return(NULL)

  dt <- data.table(
    peak = rownames(peak_counts_for_da)[keep],
    p_val = NA_real_,
    avg_log2FC = log2((mean_1[keep] + 1e-6) / (mean_2[keep] + 1e-6)),
    pct.1 = pct_1[keep],
    pct.2 = pct_2[keep],
    p_val_adj = NA_real_,
    predicted_label = label
  )
  setorder(dt, -avg_log2FC, -pct.1)
  head(dt, 5000L)
})
da_dt <- rbindlist(da_list, use.names = TRUE, fill = TRUE)
peak_counts_for_da <- NULL
gc(verbose = FALSE)
fwrite(da_dt, file.path(tab_dir, "fig4_differential_peaks_by_label.tsv.gz"), sep = "\t", quote = FALSE)
if (nrow(da_dt) == 0L) stop("No significant positive differential peaks were found for Fig. 4C annotation.")

if (nrow(da_dt) > 0L) {
  peak_gr <- StringToGRanges(unique(da_dt$peak), sep = c(":", "-"))
  peak_anno <- ChIPseeker::annotatePeak(
    peak_gr,
    tssRegion = c(-3000, 3000),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db"
  )
  anno_dt <- as.data.table(as.data.frame(peak_anno))
  anno_dt[, broad_annotation := fifelse(
    startsWith(annotation, "Promoter"),
    "Promoter",
    fifelse(
      startsWith(annotation, "Intron"),
      "Intron",
      fifelse(
        startsWith(annotation, "Exon"),
        "Exon",
        fifelse(
          startsWith(annotation, "Downstream"),
          "Downstream",
          annotation
        )
      )
    )
  )]
  anno_counts <- anno_dt[, .N, by = .(annotation = broad_annotation)][order(-N)]
  anno_counts[, fraction := N / sum(N)]
  fwrite(anno_dt, file.path(tab_dir, "fig4c_chipseeker_peak_annotation_full.tsv.gz"), sep = "\t", quote = FALSE)
  fwrite(anno_counts, file.path(tab_dir, "fig4c_chipseeker_annotation_counts.tsv"), sep = "\t", quote = FALSE)

  anno_counts[, annotation := factor(annotation, levels = annotation)]
  p_pie <- ggplot(anno_counts, aes(x = "", y = fraction, fill = annotation)) +
    geom_col(width = 1, color = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    guides(fill = guide_legend(ncol = 1, title = "Annotation")) +
    labs(title = "Fig. 4C ChIPSeeker annotation", x = NULL, y = NULL) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  save_plot(p_pie, "fig4c_chipseeker_annotation_pie", width = 6.2, height = 4.8)
}

if (run_scope == "fig4c") {
  finalize_run(requested_motifs_found = NA_integer_)
  quit(save = "no", status = 0)
}

pfm <- TFBSTools::getMatrixSet(
  JASPAR2020::JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)
combined <- AddMotifs(combined, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
combined <- RunChromVAR(combined, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
motif_obj <- Motifs(combined)
motif_name_vector <- motif_obj@motif.names
motif_lookup <- data.table(motif_id = names(motif_name_vector), motif_name = as.character(motif_name_vector))
target_motifs <- c("TFEC", "TFE3", "BHLHE41", "BACH1", "Nfe2l2", "FOSL2")
selected_motifs <- motif_lookup[motif_name %in% target_motifs]
fwrite(motif_lookup, file.path(tab_dir, "fig4_motif_lookup.tsv"), sep = "\t", quote = FALSE)
fwrite(selected_motifs, file.path(tab_dir, "fig4_requested_motif_ids.tsv"), sep = "\t", quote = FALSE)
if (nrow(selected_motifs) == 0L) stop("No requested motif names found in JASPAR motif lookup.")
missing_requested_motifs <- setdiff(target_motifs, selected_motifs$motif_name)
if (length(missing_requested_motifs) > 0L) {
  stop(
    "Requested motif names not found exactly: ",
    paste(missing_requested_motifs, collapse = ", "),
    ". Upload tables/fig4_motif_lookup.tsv and do not substitute motif names automatically."
  )
}

chromvar_data <- get_assay_data(combined, "chromvar", "data")
motif_activity <- rbindlist(lapply(seq_len(nrow(selected_motifs)), function(i) {
  motif_id <- selected_motifs$motif_id[i]
  motif_name <- selected_motifs$motif_name[i]
  values <- as.numeric(chromvar_data[motif_id, colnames(combined)])
  data.table(cell_id = colnames(combined), motif_id = motif_id, motif_name = motif_name, chromvar_z = values)
}), use.names = TRUE)
motif_activity <- merge(motif_activity, meta_dt[, .(cell_id, predicted_label)], by = "cell_id", all.x = TRUE)
motif_by_label <- motif_activity[, .(mean_chromvar_z = mean(chromvar_z, na.rm = TRUE)), by = .(motif_id, motif_name, predicted_label)]
fwrite(motif_by_label, file.path(tab_dir, "fig4e_chromvar_motif_activity_by_label.tsv"), sep = "\t", quote = FALSE)

DefaultAssay(combined) <- "chromvar"
p_motif_umaps <- FeaturePlot(
  combined,
  features = selected_motifs$motif_id,
  reduction = "umap",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  ncol = 3
)
save_plot(p_motif_umaps, "fig4h_motif_accessibility_umaps", width = 8.2, height = 5.4)

logo_list <- lapply(selected_motifs$motif_id, function(id) {
  as.matrix(TFBSTools::profileMatrix(pfm[[id]]))
})
names(logo_list) <- selected_motifs$motif_name
p_logo <- ggseqlogo::ggseqlogo(logo_list, ncol = 3) +
  labs(title = "Fig. 4I motif logos") +
  theme(plot.title = element_text(face = "bold", size = 10))
save_plot(p_logo, "fig4i_motif_logos", width = 8.2, height = 4.8)

DefaultAssay(combined) <- "peaks"
track_genes <- c("C1QB", "CHID1", "IL1B", "SLC25A37")
track_summary <- data.table(gene = track_genes, status = "requested")
fwrite(track_summary, file.path(tab_dir, "fig4g_track_region_summary.tsv"), sep = "\t", quote = FALSE)
track_plots <- lapply(track_genes, function(gene) {
  CoveragePlot(
    object = combined,
    region = gene,
    group.by = "predicted_label",
    extend.upstream = 8000,
    extend.downstream = 8000,
    peaks = TRUE,
    annotation = TRUE
  ) + ggtitle(gene)
})
p_tracks <- wrap_plots(track_plots, ncol = 2)
save_plot(p_tracks, "fig4g_gene_coverage_tracks", width = 10.5, height = 7.5)

finalize_run(requested_motifs_found = nrow(selected_motifs))
