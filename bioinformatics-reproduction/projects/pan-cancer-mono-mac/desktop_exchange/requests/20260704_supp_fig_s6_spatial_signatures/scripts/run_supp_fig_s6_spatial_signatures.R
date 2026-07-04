#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(scales)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", file_arg[1] %||% "")
if (nzchar(script_path)) {
  request_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
  request_root <- normalizePath(Sys.getenv("S6_REQUEST_DIR", unset = getwd()), mustWork = FALSE)
}
if (!dir.exists(file.path(request_root, "scripts"))) {
  request_root <- normalizePath(Sys.getenv("S6_REQUEST_DIR", unset = getwd()), mustWork = FALSE)
}

data_dir <- Sys.getenv("S6_DATA_DIR", unset = file.path(request_root, "data"))
out_dir <- Sys.getenv("S6_OUTPUT_DIR", unset = file.path(request_root, "outputs"))
signature_path <- file.path(request_root, "inputs", "fig5_signature_gene_sets_unique_current_scope.tsv")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(log_dir, "run_supp_fig_s6_spatial_signatures.log")

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
  message(step, " [", status, "] ", note)
}

options(error = function() {
  write_log("error", "failed", trimws(geterrmessage()))
  traceback(2, max.lines = 20)
  quit(save = "no", status = 1)
})

get_assay_data <- function(object, assay, layer_name) {
  tryCatch(
    GetAssayData(object, assay = assay, layer = layer_name),
    error = function(e) GetAssayData(object, assay = assay, slot = layer_name)
  )
}

download_file <- function(url, destfile, expected_size = NA_real_) {
  if (file.exists(destfile) && (is.na(expected_size) || file.info(destfile)$size == expected_size)) {
    write_log("download", "exists", destfile)
    return(destfile)
  }
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  write_log("download", "running", paste(url, "->", destfile))
  status <- tryCatch({
    utils::download.file(url, destfile, mode = "wb", quiet = FALSE)
    0L
  }, error = function(e) {
    write_log("download", "error", paste(destfile, conditionMessage(e)))
    1L
  })
  if (!identical(status, 0L)) stop("Download failed: ", url)
  if (!is.na(expected_size) && file.info(destfile)$size != expected_size) {
    stop("Downloaded file has unexpected size: ", destfile, "; expected=", expected_size, "; observed=", file.info(destfile)$size)
  }
  destfile
}

download_10x_spatial <- function(sample_id, base_url) {
  sample_dir <- file.path(data_dir, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  h5_url <- paste0(base_url, "/", sample_id, "_filtered_feature_bc_matrix.h5")
  spatial_url <- paste0(base_url, "/", sample_id, "_spatial.tar.gz")
  image_url <- paste0(base_url, "/", sample_id, "_image.tif")
  download_file(h5_url, file.path(sample_dir, "filtered_feature_bc_matrix.h5"))
  spatial_tar <- download_file(spatial_url, file.path(sample_dir, "spatial.tar.gz"))
  download_file(image_url, file.path(sample_dir, paste0(sample_id, "_image.tif")))
  if (!dir.exists(file.path(sample_dir, "spatial"))) {
    utils::untar(spatial_tar, exdir = sample_dir)
  }
  sample_dir
}

safe_untar <- function(archive, exdir) {
  marker <- file.path(exdir, paste0(".extracted_", basename(archive), ".done"))
  if (file.exists(marker)) return(invisible(TRUE))
  write_log("extract", "running", archive)
  utils::untar(archive, exdir = exdir)
  file.create(marker)
  invisible(TRUE)
}

extract_nested_archives <- function(root_dir, max_rounds = 3L) {
  for (round_i in seq_len(max_rounds)) {
    archives <- list.files(root_dir, pattern = "\\.(tar|tar.gz|tgz)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(archives) == 0L) break
    extracted_any <- FALSE
    for (archive in archives) {
      marker <- file.path(dirname(archive), paste0(".extracted_", basename(archive), ".done"))
      if (!file.exists(marker)) {
        safe_untar(archive, dirname(archive))
        extracted_any <- TRUE
      }
    }
    if (!extracted_any) break
  }
}

download_geo_raw <- function(accession, expected_size = NA_real_) {
  acc_dir <- file.path(data_dir, accession)
  dir.create(acc_dir, recursive = TRUE, showWarnings = FALSE)
  tar_path <- file.path(acc_dir, paste0(accession, "_RAW.tar"))
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", accession, "&format=file")
  download_file(url, tar_path, expected_size = expected_size)
  safe_untar(tar_path, acc_dir)
  extract_nested_archives(acc_dir)
  acc_dir
}

find_visium_dirs <- function(root_dir) {
  spatial_dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
  spatial_dirs <- spatial_dirs[basename(spatial_dirs) == "spatial"]
  parents <- unique(dirname(spatial_dirs))
  parents[vapply(parents, function(p) {
    any(file.exists(file.path(p, c("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5")))) ||
      file.exists(file.path(p, "filtered_feature_bc_matrix", "matrix.mtx.gz")) ||
      file.exists(file.path(p, "raw_feature_bc_matrix", "matrix.mtx.gz")) ||
      length(list.files(p, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)) > 0L
  }, logical(1))]
}

standardize_visium_dir <- function(sample_dir, label) {
  target <- file.path(data_dir, paste0("standardized_s6_", label))
  if (dir.exists(target)) unlink(target, recursive = TRUE)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  spatial_src <- file.path(sample_dir, "spatial")
  if (!dir.exists(spatial_src)) stop("No spatial directory: ", sample_dir)
  file.copy(spatial_src, target, recursive = TRUE)
  h5 <- file.path(sample_dir, "filtered_feature_bc_matrix.h5")
  if (!file.exists(h5)) h5 <- file.path(sample_dir, "raw_feature_bc_matrix.h5")
  if (file.exists(h5)) {
    file.copy(h5, file.path(target, "filtered_feature_bc_matrix.h5"))
    return(target)
  }
  matrix_dirs <- unique(dirname(list.files(sample_dir, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)))
  if (length(matrix_dirs) == 0L) stop("No matrix file found under: ", sample_dir)
  file.copy(matrix_dirs[1], file.path(target, "filtered_feature_bc_matrix"), recursive = TRUE)
  target
}

load_spatial_object <- function(sample_dir, dataset_id) {
  write_log("load_spatial", "running", paste(dataset_id, sample_dir))
  obj <- Load10X_Spatial(data.dir = sample_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
  DefaultAssay(obj) <- "Spatial"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj$spatial_dataset_id <- dataset_id
  obj
}

signature_score <- function(expr, genes, score_name) {
  measured <- intersect(genes, rownames(expr))
  if (length(measured) == 0L) stop(score_name, " has no measured genes in spatial object.")
  score <- Matrix::colMeans(expr[measured, , drop = FALSE])
  list(score = score, measured = measured, missing = setdiff(genes, measured))
}

extract_plot_data <- function(obj, dataset_id, tumor_marker, spp1_genes, mdsc_genes) {
  DefaultAssay(obj) <- "Spatial"
  expr <- get_assay_data(obj, "Spatial", "data")
  coords <- as.data.table(GetTissueCoordinates(obj))
  if (!all(c("x", "y", "cell") %in% names(coords))) {
    stop(dataset_id, " tissue coordinate columns are unsupported. Observed columns: ", paste(names(coords), collapse = ", "))
  }
  image_obj <- obj@images[[1]]
  lowres_scale <- image_obj@scale.factors$lowres %||% 1
  raster_img <- as.raster(image_obj@image)
  img_height <- dim(image_obj@image)[1]
  img_width <- dim(image_obj@image)[2]
  coords[, x_plot := x * lowres_scale]
  coords[, y_plot := img_height - y * lowres_scale]

  spp1 <- signature_score(expr, spp1_genes, paste0(dataset_id, " SPP1+ TAMs top10"))
  mdsc <- signature_score(expr, mdsc_genes, paste0(dataset_id, " MDSC 19 genes"))
  tumor_values <- if (tumor_marker %in% rownames(expr)) as.numeric(expr[tumor_marker, coords$cell]) else rep(NA_real_, nrow(coords))
  dt <- copy(coords)
  dt[, dataset_id := dataset_id]
  dt[, tumor_marker := tumor_marker]
  dt[, tumor_marker_expr := tumor_values]
  dt[, spp1_tam_signature := as.numeric(spp1$score[cell])]
  dt[, mdsc_signature := as.numeric(mdsc$score[cell])]
  list(
    dt = dt,
    raster_img = raster_img,
    img_width = img_width,
    img_height = img_height,
    coverage = rbindlist(list(
      data.table(dataset_id = dataset_id, signature = "tumor_marker", gene = tumor_marker, measured = tumor_marker %in% rownames(expr)),
      data.table(dataset_id = dataset_id, signature = "SPP1+ TAMs top10", gene = spp1_genes, measured = spp1_genes %in% rownames(expr)),
      data.table(dataset_id = dataset_id, signature = "MDSC 19 genes", gene = mdsc_genes, measured = mdsc_genes %in% rownames(expr))
    ), fill = TRUE)
  )
}

feature_palette <- c("#2C338B", "#2CA25F", "#FDE725", "#F98E52", "#B40426")

make_spatial_plot <- function(dt, raster_img, img_width, img_height, value_column, title) {
  plot_dt <- copy(dt)
  plot_dt[, value := as.numeric(get(value_column))]
  finite_values <- plot_dt[is.finite(value), value]
  if (length(finite_values) == 0L) {
    return(ggplot() + annotate("text", x = 0, y = 0, label = paste(title, "missing"), size = 4) + theme_void())
  }
  lower <- as.numeric(quantile(finite_values, probs = 0.005, na.rm = TRUE))
  upper <- as.numeric(quantile(finite_values, probs = 0.995, na.rm = TRUE))
  if (!is.finite(lower)) lower <- 0
  if (!is.finite(upper) || upper <= lower) upper <- max(finite_values, na.rm = TRUE)
  if (!is.finite(upper) || upper <= lower) upper <- lower + 1

  ggplot(plot_dt, aes(x_plot, y_plot)) +
    annotation_raster(raster_img, xmin = 0, xmax = img_width, ymin = 0, ymax = img_height) +
    geom_point(aes(color = value), size = 0.38, alpha = 0.95, stroke = 0) +
    scale_color_gradientn(colors = feature_palette, limits = c(lower, upper), oob = scales::squish) +
    coord_fixed(xlim = c(0, img_width), ylim = c(0, img_height), expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
      legend.position = "right",
      legend.key.height = unit(11, "mm"),
      legend.key.width = unit(2, "mm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 5),
      plot.margin = margin(2, 2, 2, 2)
    )
}

format_p <- function(p_value) {
  if (!is.finite(p_value)) return("p = NA")
  if (p_value < 2.2e-16) return("p < 2.2e-16")
  if (p_value < 1e-4) return(paste0("p = ", formatC(p_value, format = "e", digits = 1)))
  paste0("p = ", signif(p_value, 3))
}

make_scatter <- function(dt, dataset_label) {
  ct <- suppressWarnings(cor.test(dt$spp1_tam_signature, dt$mdsc_signature, method = "spearman", exact = FALSE))
  label <- paste0("R = ", signif(unname(ct$estimate), 2), ", ", format_p(ct$p.value))
  p <- ggplot(dt, aes(spp1_tam_signature, mdsc_signature)) +
    geom_point(color = "#3045A6", size = 0.55, alpha = 0.60, stroke = 0) +
    geom_smooth(method = "lm", se = FALSE, color = "#E31A1C", linewidth = 0.9) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 1.25, label = label, size = 2.6) +
    labs(x = "SPP1+ TAMs signature", y = "MDSC signature", title = dataset_label) +
    theme_classic(base_size = 7) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8))
  list(
    plot = p,
    stats = data.table(dataset_id = unique(dt$dataset_id), dataset_label = dataset_label, spearman_rho = unname(ct$estimate), p_value = ct$p.value, spots = nrow(dt))
  )
}

write_log("start", "running", "Supplementary Fig. S6 spatial signature workflow")
if (!file.exists(signature_path)) stop("Missing signature table: ", signature_path)
sig <- fread(signature_path)
if (!all(c("gene", "components") %in% names(sig))) {
  stop("Signature table missing required columns. Observed columns: ", paste(names(sig), collapse = ", "))
}
spp1_genes <- sig[grepl("SPP1+ TAMs top10", components, fixed = TRUE), unique(gene)]
mdsc_genes <- sig[grepl("MDSC 19 genes", components, fixed = TRUE), unique(gene)]
if (length(spp1_genes) != 10L) stop("Expected 10 SPP1+ TAMs top10 genes, observed ", length(spp1_genes))
if (length(mdsc_genes) != 19L) stop("Expected 19 MDSC genes, observed ", length(mdsc_genes))

dataset_plan <- data.table(
  dataset_id = c("10x_GBM", "GSE226997_CRC", "GSE274103_PDAC", "10x_CRC"),
  dataset_label = c("10x GBM", "GSE226997 Colorectal cancer", "GSE274103 PDAC", "10x Colorectal cancer"),
  tumor_marker = c("EGFR", "EPCAM", "EPCAM", "EPCAM"),
  source = c("10x", "GEO", "GEO", "10x")
)
fwrite(dataset_plan, file.path(tab_dir, "supp_fig_s6_dataset_manifest.tsv"), sep = "\t", quote = FALSE)

objects <- list()
selected_sources <- list()

gbm_id <- "CytAssist_11mm_FFPE_Human_Glioblastoma"
gbm_url <- "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma"
gbm_dir <- download_10x_spatial(gbm_id, gbm_url)
objects[["10x_GBM"]] <- load_spatial_object(gbm_dir, "10x_GBM")
selected_sources[["10x_GBM"]] <- data.table(dataset_id = "10x_GBM", selected_path = gbm_dir, status = "loaded")

crc_id <- "CytAssist_11mm_FFPE_Human_Colorectal_Cancer"
crc_url <- "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Colorectal_Cancer"
crc_dir <- download_10x_spatial(crc_id, crc_url)
objects[["10x_CRC"]] <- load_spatial_object(crc_dir, "10x_CRC")
selected_sources[["10x_CRC"]] <- data.table(dataset_id = "10x_CRC", selected_path = crc_dir, status = "loaded")

if (Sys.getenv("S6_SKIP_GSE274103", "0") != "1") {
  gse274_dir <- download_geo_raw("GSE274103", expected_size = 186603520)
  dirs <- sort(find_visium_dirs(gse274_dir))
  fwrite(data.table(dataset_id = "GSE274103_PDAC", discovered_path = dirs), file.path(tab_dir, "supp_fig_s6_gse274103_discovered_visium_dirs.tsv"), sep = "\t", quote = FALSE)
  selected <- Sys.getenv("S6_GSE274103_SAMPLE", "")
  if (selected == "") selected <- dirs[1]
  if (is.na(selected) || !dir.exists(selected)) stop("No usable GSE274103 Visium directory discovered.")
  std <- standardize_visium_dir(selected, "GSE274103_PDAC")
  objects[["GSE274103_PDAC"]] <- load_spatial_object(std, "GSE274103_PDAC")
  selected_sources[["GSE274103_PDAC"]] <- data.table(dataset_id = "GSE274103_PDAC", selected_path = selected, status = "loaded")
} else {
  selected_sources[["GSE274103_PDAC"]] <- data.table(dataset_id = "GSE274103_PDAC", selected_path = NA_character_, status = "skipped_by_env")
}

if (Sys.getenv("S6_SKIP_GSE226997", "0") != "1") {
  gse226_dir <- download_geo_raw("GSE226997")
  dirs <- sort(find_visium_dirs(gse226_dir))
  fwrite(data.table(dataset_id = "GSE226997_CRC", discovered_path = dirs), file.path(tab_dir, "supp_fig_s6_gse226997_discovered_visium_dirs.tsv"), sep = "\t", quote = FALSE)
  selected <- Sys.getenv("S6_GSE226997_SAMPLE", "")
  if (selected == "") selected <- dirs[1]
  if (is.na(selected) || !dir.exists(selected)) stop("No usable GSE226997 Visium directory discovered.")
  std <- standardize_visium_dir(selected, "GSE226997_CRC")
  objects[["GSE226997_CRC"]] <- load_spatial_object(std, "GSE226997_CRC")
  selected_sources[["GSE226997_CRC"]] <- data.table(dataset_id = "GSE226997_CRC", selected_path = selected, status = "loaded")
} else {
  selected_sources[["GSE226997_CRC"]] <- data.table(dataset_id = "GSE226997_CRC", selected_path = NA_character_, status = "skipped_by_env")
}

selected_sources_dt <- rbindlist(selected_sources, fill = TRUE)
fwrite(selected_sources_dt, file.path(tab_dir, "supp_fig_s6_selected_spatial_sources.tsv"), sep = "\t", quote = FALSE)

plot_data <- list()
coverage <- list()
for (id in names(objects)) {
  marker <- dataset_plan[dataset_id == id, tumor_marker][1]
  plot_data[[id]] <- extract_plot_data(objects[[id]], id, marker, spp1_genes, mdsc_genes)
  coverage[[id]] <- plot_data[[id]]$coverage
}
coverage_dt <- rbindlist(coverage, fill = TRUE)
fwrite(coverage_dt, file.path(tab_dir, "supp_fig_s6_signature_gene_coverage.tsv"), sep = "\t", quote = FALSE)
gene_presence <- coverage_dt[signature == "tumor_marker", .(dataset_id, gene, present = measured)]
fwrite(gene_presence, file.path(tab_dir, "supp_fig_s6_gene_presence.tsv"), sep = "\t", quote = FALSE)

score_dt <- rbindlist(lapply(plot_data, `[[`, "dt"), fill = TRUE)
fwrite(score_dt, file.path(tab_dir, "supp_fig_s6_signature_scores.tsv.gz"), sep = "\t", quote = FALSE)

row_plots <- list()
scatter_plots <- list()
cor_stats <- list()
for (id in dataset_plan$dataset_id) {
  label <- dataset_plan[dataset_id == id, dataset_label][1]
  marker <- dataset_plan[dataset_id == id, tumor_marker][1]
  if (!(id %in% names(plot_data))) {
    row_plots[[id]] <- ggplot() + annotate("text", x = 0, y = 0, label = paste(label, "not loaded"), size = 5) + theme_void()
    scatter_plots[[id]] <- ggplot() + annotate("text", x = 0, y = 0, label = paste(label, "not loaded"), size = 5) + theme_void()
    next
  }
  pd <- plot_data[[id]]
  row_plots[[id]] <- plot_grid(
    make_spatial_plot(pd$dt, pd$raster_img, pd$img_width, pd$img_height, "tumor_marker_expr", marker),
    make_spatial_plot(pd$dt, pd$raster_img, pd$img_width, pd$img_height, "spp1_tam_signature", "SPP1+ TAMs signature"),
    make_spatial_plot(pd$dt, pd$raster_img, pd$img_width, pd$img_height, "mdsc_signature", "MDSC signature"),
    nrow = 1,
    labels = NULL
  )
  row_plots[[id]] <- plot_grid(
    ggdraw() + draw_label(label, angle = 90, fontface = "bold", size = 8),
    row_plots[[id]],
    nrow = 1,
    rel_widths = c(0.08, 1)
  )
  sc <- make_scatter(pd$dt, label)
  scatter_plots[[id]] <- sc$plot
  cor_stats[[id]] <- sc$stats
}
cor_dt <- rbindlist(cor_stats, fill = TRUE)
fwrite(cor_dt, file.path(tab_dir, "supp_fig_s6_spearman_correlations.tsv"), sep = "\t", quote = FALSE)

spatial_panel <- wrap_plots(row_plots, ncol = 1)
scatter_panel <- wrap_plots(scatter_plots, ncol = 4)
combined <- wrap_plots(spatial_panel, scatter_panel, ncol = 1, heights = c(4, 1.35)) +
  plot_annotation(
    title = "Supplementary Fig. S6 spatial THBS1+ MDSC - SPP1+ TAM lineage",
    subtitle = "SPP1+ TAM top10 and MDSC 19-gene signatures; skipped datasets are shown as placeholders.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 9, color = "#555555")
    )
  )

png_path <- file.path(fig_dir, "supp_fig_s6_spatial_signatures.png")
pdf_path <- file.path(fig_dir, "supp_fig_s6_spatial_signatures.pdf")
ggsave(png_path, combined, width = 12.5, height = 15.0, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(pdf_path, combined, width = 12.5, height = 15.0, bg = "white", limitsize = FALSE)
write_log("save_figure", "complete", png_path)

run_summary <- data.table(
  metric = c("datasets_loaded", "datasets_planned", "spp1_signature_genes", "mdsc_signature_genes", "correlation_panels"),
  value = c(length(objects), nrow(dataset_plan), length(spp1_genes), length(mdsc_genes), nrow(cor_dt))
)
fwrite(run_summary, file.path(tab_dir, "supp_fig_s6_run_summary.tsv"), sep = "\t", quote = FALSE)
write_log("done", "complete", paste("datasets_loaded=", length(objects)))
