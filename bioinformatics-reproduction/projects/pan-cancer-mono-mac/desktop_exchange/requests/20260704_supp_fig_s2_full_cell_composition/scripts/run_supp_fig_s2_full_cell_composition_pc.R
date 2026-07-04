#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(irlba)
  library(uwot)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", file_arg[1] %||% "")
request_root <- if (nzchar(script_path)) normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE) else normalizePath(getwd(), mustWork = FALSE)
if (!dir.exists(file.path(request_root, "scripts"))) request_root <- normalizePath(getwd(), mustWork = FALSE)

data_dir <- file.path(request_root, "data")
out_dir <- file.path(request_root, "outputs")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(log_dir, "run_supp_fig_s2_full_cell_composition_pc.log")

write_log <- function(step, status, note = "") {
  row <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    step = step,
    status = status,
    note = note
  )
  fwrite(row, log_path, sep = "\t", quote = FALSE, append = file.exists(log_path), col.names = !file.exists(log_path))
  message(step, " [", status, "] ", note)
}

download_file <- function(url, destfile) {
  if (file.exists(destfile) && file.info(destfile)$size > 0) {
    write_log("download", "exists", destfile)
    return(destfile)
  }
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  write_log("download", "running", paste(url, "->", destfile))
  utils::download.file(url, destfile, mode = "wb", quiet = FALSE)
  if (!file.exists(destfile) || file.info(destfile)$size == 0) stop("Downloaded file is missing or empty: ", destfile)
  destfile
}

download_geo_files <- function(accession, files) {
  base <- sub("^(GSE[0-9]+).*$", "\\1", accession)
  n <- as.integer(sub("^GSE([0-9]+)$", "\\1", base))
  series_dir <- sprintf("GSE%03dnnn", floor(n / 1000))
  url_base <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", series_dir, "/", base, "/suppl/")
  acc_dir <- file.path(data_dir, base)
  dir.create(acc_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(acc_dir, files)
  for (i in seq_along(files)) download_file(paste0(url_base, files[i]), out[i])
  out
}

safe_untar <- function(archive, exdir) {
  marker <- file.path(exdir, paste0(".extracted_", basename(archive), ".done"))
  if (file.exists(marker)) return(invisible(TRUE))
  write_log("extract", "running", archive)
  utils::untar(archive, exdir = exdir)
  file.create(marker)
  invisible(TRUE)
}

cell_type_palette <- c(
  "T cells" = "#8F8AC2",
  "Myeloid cells" = "#F4A340",
  "Epithelial cells" = "#9FD44F",
  "B cells" = "#5BC4B7",
  "Endothelial cells" = "#F4A9CC",
  "Fibroblasts" = "#F2D94E",
  "Mast cells" = "#E77761",
  "Plasma cells" = "#A65AB6",
  "NK cells" = "#9BC7E5",
  "Prol. epithelial cells" = "#A9DFA3",
  "pDCs" = "#6FA8D6",
  "Other cells" = "#CFCFCF"
)
display_order <- names(cell_type_palette)

run_seurat_umap <- function(counts, meta, npcs = 30) {
  meta_df <- as.data.frame(meta)
  rownames(meta_df) <- meta_df$cell_id
  obj <- CreateSeuratObject(counts = counts, meta.data = meta_df)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs, verbose = FALSE)
  obj <- RunUMAP(obj, dims = seq_len(min(npcs, ncol(Embeddings(obj, "pca")))), n.neighbors = 30, min.dist = 0.35, verbose = FALSE)
  emb <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "cell_id")
  umap_cols <- setdiff(names(emb), "cell_id")
  setnames(emb, umap_cols[1:2], c("UMAP_1", "UMAP_2"))
  emb[, global_cell_type := obj$global_cell_type[cell_id]]
  list(object = obj, umap = emb[])
}

composition_table <- function(meta, dataset) {
  out <- as.data.table(meta)[, .N, by = global_cell_type]
  out[, dataset := dataset]
  out[, percent := 100 * N / sum(N)]
  setcolorder(out, c("dataset", "global_cell_type", "N", "percent"))
  out[]
}

make_umap_plot <- function(dt, title) {
  d <- copy(dt)
  d[, global_cell_type := factor(global_cell_type, levels = display_order)]
  labs <- d[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), n = .N), by = global_cell_type]
  labs <- labs[!is.na(global_cell_type) & global_cell_type != "Other cells" & n >= max(20, 0.004 * nrow(d))]
  ggplot(d, aes(UMAP_1, UMAP_2, color = global_cell_type)) +
    geom_point(size = 0.10, alpha = 0.78, stroke = 0) +
    ggrepel::geom_text_repel(
      data = labs,
      aes(UMAP_1, UMAP_2, label = as.character(global_cell_type)),
      inherit.aes = FALSE,
      size = 2.4,
      fontface = ifelse(as.character(labs$global_cell_type) == "Myeloid cells", "bold", "plain"),
      color = "black",
      box.padding = 0.14,
      point.padding = 0.05,
      min.segment.length = 0.02,
      segment.size = 0.14,
      max.overlaps = Inf,
      seed = 20260704
    ) +
    scale_color_manual(values = cell_type_palette, drop = FALSE) +
    labs(title = title, x = "UMAP1", y = "UMAP2") +
    coord_equal() +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
      legend.position = "none",
      axis.title.x = element_text(size = 7, color = "black"),
      axis.title.y = element_text(size = 7, color = "black", angle = 90),
      plot.margin = margin(2, 2, 2, 2)
    )
}

make_pie_plot <- function(comp) {
  d <- copy(comp)
  d[, global_cell_type := factor(global_cell_type, levels = display_order)]
  setorder(d, global_cell_type)
  d[, label := ifelse(percent >= 5, sprintf("%.2f%%", percent), "")]
  ggplot(d, aes(x = 1, y = percent, fill = global_cell_type)) +
    geom_col(width = 1, color = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.55) +
    scale_fill_manual(values = cell_type_palette, drop = FALSE) +
    theme_void(base_size = 8) +
    theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))
}

write_source_manifest <- function() {
  manifest <- data.table(
    dataset = c(
      rep("LUAD_GSE131907", 2),
      rep("LIHC_CHOL_GSE125449", 8),
      "Healthy_lung_GSE122960",
      rep("Healthy_liver_GSE115469", 2),
      "NSCLC_GSE148071",
      rep("HCC_GSE156625", 4)
    ),
    file = c(
      "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
      "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
      "GSE125449_Set1_matrix.mtx.gz",
      "GSE125449_Set1_genes.tsv.gz",
      "GSE125449_Set1_barcodes.tsv.gz",
      "GSE125449_Set1_samples.txt.gz",
      "GSE125449_Set2_matrix.mtx.gz",
      "GSE125449_Set2_genes.tsv.gz",
      "GSE125449_Set2_barcodes.tsv.gz",
      "GSE125449_Set2_samples.txt.gz",
      "GSE122960_RAW.tar",
      "GSE115469_CellClusterType.txt.gz",
      "GSE115469_Data.csv.gz",
      "GSE148071_RAW.tar",
      "GSE156625_HCCmatrix.mtx.gz",
      "GSE156625_HCCbarcodes.tsv.gz",
      "GSE156625_HCCgenes.tsv.gz",
      "GSE156625_HCCscanpyobj.h5ad.gz"
    )
  )
  fwrite(manifest, file.path(tab_dir, "supp_fig_s2_dataset_manifest.tsv"), sep = "\t", quote = FALSE)
}

run_lihc_chol <- function() {
  write_log("LIHC-CHOL", "running")
  files <- c(
    "GSE125449_Set1_matrix.mtx.gz", "GSE125449_Set1_genes.tsv.gz", "GSE125449_Set1_barcodes.tsv.gz", "GSE125449_Set1_samples.txt.gz",
    "GSE125449_Set2_matrix.mtx.gz", "GSE125449_Set2_genes.tsv.gz", "GSE125449_Set2_barcodes.tsv.gz", "GSE125449_Set2_samples.txt.gz"
  )
  download_geo_files("GSE125449", files)
  base <- file.path(data_dir, "GSE125449")
  read_one <- function(set_id) {
    mtx <- Matrix::readMM(file.path(base, paste0("GSE125449_", set_id, "_matrix.mtx.gz")))
    genes <- fread(file.path(base, paste0("GSE125449_", set_id, "_genes.tsv.gz")), header = FALSE)
    barcodes <- fread(file.path(base, paste0("GSE125449_", set_id, "_barcodes.tsv.gz")), header = FALSE)
    samples <- fread(file.path(base, paste0("GSE125449_", set_id, "_samples.txt.gz")))
    rownames(mtx) <- make.unique(as.character(genes[[1]]))
    colnames(mtx) <- paste(set_id, as.character(barcodes[[1]]), sep = "__")
    samples[, cell_id := paste(set_id, `Cell Barcode`, sep = "__")]
    samples[, set_id := set_id]
    list(counts = mtx, meta = samples)
  }
  s1 <- read_one("Set1")
  s2 <- read_one("Set2")
  genes <- intersect(rownames(s1$counts), rownames(s2$counts))
  counts <- cbind(s1$counts[genes, ], s2$counts[genes, ])
  meta <- rbindlist(list(s1$meta, s2$meta), fill = TRUE)
  meta[, global_cell_type := fcase(
    Type == "T cell", "T cells",
    Type == "TAM", "Myeloid cells",
    Type %in% c("Malignant cell", "HPC-like"), "Epithelial cells",
    Type == "B cell", "B cells",
    Type == "TEC", "Endothelial cells",
    Type == "CAF", "Fibroblasts",
    default = "Other cells"
  )]
  res <- run_seurat_umap(counts, meta[, .(cell_id, global_cell_type)])
  list(dataset = "LIHC-CHOL", meta = meta[, .(cell_id, global_cell_type)], umap = res$umap, object = res$object)
}

write_log("start", "running", "Supplementary Fig S2 PC request scaffold.")
write_source_manifest()

annotation_audit <- data.table(
  dataset = c("GSE125449", "GSE131907", "GSE115469", "GSE148071", "GSE156625"),
  exact_annotation_source_to_check = c(
    "GSE125449_Set1_samples.txt.gz and GSE125449_Set2_samples.txt.gz column Type",
    "GSE131907_Lung_Cancer_cell_annotation.txt.gz columns Cell_type, Cell_type.refined, Cell_subtype, Sample_Origin",
    "GSE115469_CellClusterType.txt.gz column CellType",
    "Inspect extracted GSE148071 *_exp.txt.gz files; GEO processed expression files may not include cell type labels",
    "Inspect h5ad obs columns from GSE156625_HCCscanpyobj.h5ad.gz before choosing an annotation field"
  )
)
fwrite(annotation_audit, file.path(tab_dir, "supp_fig_s2_annotation_columns_audit.tsv"), sep = "\t", quote = FALSE)

if (Sys.getenv("S2_RUN_LIHC_ONLY", "0") == "1") {
  result <- run_lihc_chol()
  comp <- composition_table(result$meta, result$dataset)
  fwrite(comp, file.path(tab_dir, "supp_fig_s2_cell_type_composition_lihc_chol.tsv"), sep = "\t", quote = FALSE)
  fwrite(result$umap, file.path(tab_dir, "supp_fig_s2_umap_coordinates_lihc_chol.tsv.gz"), sep = "\t", quote = FALSE)
  p <- cowplot::plot_grid(make_umap_plot(result$umap, result$dataset), make_pie_plot(comp), rel_widths = c(1.2, 0.8), nrow = 1)
  ggsave(file.path(fig_dir, "supp_fig_s2_lihc_chol_smoke_test.png"), p, width = 6, height = 3.2, dpi = 300, bg = "white")
}

write_log("manual_completion_needed", "pending", "Complete LUAD, healthy lung, healthy liver, NSCLC-GSE148071, and HCC-GSE156625 as specified in README, then upload the required outputs.")
write_log("complete", "complete", "Scaffold finished.")

