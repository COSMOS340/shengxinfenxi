#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
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
  request_root <- normalizePath(getwd(), mustWork = FALSE)
}
if (!dir.exists(file.path(request_root, "scripts"))) {
  request_root <- normalizePath(getwd(), mustWork = FALSE)
}

data_dir <- file.path(request_root, "data")
out_dir <- file.path(request_root, "outputs")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(log_dir, "run_supp_fig_s4_s5_spatial_featureplots.log")

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
  target <- file.path(data_dir, paste0("standardized_", label))
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

plot_spatial_features <- function(obj, dataset_label, features) {
  present <- features[features %in% rownames(obj)]
  missing <- setdiff(features, present)
  if (length(missing) > 0L) {
    write_log("gene_presence", "missing", paste(dataset_label, paste(missing, collapse = ",")))
  }
  plots <- lapply(features, function(gene) {
    if (!(gene %in% rownames(obj))) {
      ggplot() +
        annotate("text", x = 0, y = 0, label = paste0(gene, "\nmissing"), size = 5) +
        labs(title = gene) +
        theme_void()
    } else {
      SpatialFeaturePlot(obj, features = gene, pt.size.factor = 1.4, alpha = c(0.15, 1), crop = TRUE) +
        scale_fill_gradientn(colors = c("#2C2C84", "#3DBB88", "#F5E84A", "#F26B2A", "#A50044"), oob = squish) +
        labs(title = gene) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 7)
        )
    }
  })
  row <- wrap_plots(plots, nrow = 1)
  row + plot_annotation(title = dataset_label, theme = theme(plot.title = element_text(angle = 90, face = "bold", size = 10)))
}

dataset_plan <- data.table(
  dataset_id = c("10x_GBM", "GSE226997_CRC", "GSE274103_PDAC", "10x_CRC"),
  dataset_label = c("10x GBM", "GSE226997 Colorectal cancer", "GSE274103 PDAC", "10x Colorectal cancer"),
  genes = c("EGFR,C1QC,FOLR2", "EPCAM,C1QC,FOLR2", "EPCAM,C1QC,PLTP", "EPCAM,C1QC,FOLR2"),
  source = c("10x", "GEO", "GEO", "10x")
)
fwrite(dataset_plan, file.path(tab_dir, "supp_fig_s4_s5_dataset_manifest.tsv"), sep = "\t", quote = FALSE)

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

if (Sys.getenv("S4S5_SKIP_GSE274103", "0") != "1") {
  gse274_dir <- download_geo_raw("GSE274103", expected_size = 186603520)
  dirs <- sort(find_visium_dirs(gse274_dir))
  fwrite(data.table(dataset_id = "GSE274103_PDAC", discovered_path = dirs), file.path(tab_dir, "supp_fig_s4_s5_gse274103_discovered_visium_dirs.tsv"), sep = "\t", quote = FALSE)
  selected <- Sys.getenv("S4S5_GSE274103_SAMPLE", "")
  if (selected == "") selected <- dirs[1]
  if (is.na(selected) || !dir.exists(selected)) stop("No usable GSE274103 Visium directory discovered.")
  std <- standardize_visium_dir(selected, "GSE274103_PDAC")
  objects[["GSE274103_PDAC"]] <- load_spatial_object(std, "GSE274103_PDAC")
  selected_sources[["GSE274103_PDAC"]] <- data.table(dataset_id = "GSE274103_PDAC", selected_path = selected, status = "loaded")
} else {
  selected_sources[["GSE274103_PDAC"]] <- data.table(dataset_id = "GSE274103_PDAC", selected_path = NA_character_, status = "skipped_by_env")
}

if (Sys.getenv("S4S5_SKIP_GSE226997", "0") != "1") {
  gse226_dir <- download_geo_raw("GSE226997")
  dirs <- sort(find_visium_dirs(gse226_dir))
  fwrite(data.table(dataset_id = "GSE226997_CRC", discovered_path = dirs), file.path(tab_dir, "supp_fig_s4_s5_gse226997_discovered_visium_dirs.tsv"), sep = "\t", quote = FALSE)
  selected <- Sys.getenv("S4S5_GSE226997_SAMPLE", "")
  if (selected == "") selected <- dirs[1]
  if (is.na(selected) || !dir.exists(selected)) stop("No usable GSE226997 Visium directory discovered.")
  std <- standardize_visium_dir(selected, "GSE226997_CRC")
  objects[["GSE226997_CRC"]] <- load_spatial_object(std, "GSE226997_CRC")
  selected_sources[["GSE226997_CRC"]] <- data.table(dataset_id = "GSE226997_CRC", selected_path = selected, status = "loaded")
} else {
  selected_sources[["GSE226997_CRC"]] <- data.table(dataset_id = "GSE226997_CRC", selected_path = NA_character_, status = "skipped_by_env")
}

selected_sources_dt <- rbindlist(selected_sources, fill = TRUE)
fwrite(selected_sources_dt, file.path(tab_dir, "supp_fig_s4_s5_selected_spatial_sources.tsv"), sep = "\t", quote = FALSE)

gene_presence <- rbindlist(lapply(names(objects), function(id) {
  genes <- strsplit(dataset_plan[dataset_id == id, genes], ",", fixed = TRUE)[[1]]
  data.table(dataset_id = id, gene = genes, present = genes %in% rownames(objects[[id]]), spots = ncol(objects[[id]]), features = nrow(objects[[id]]))
}))
fwrite(gene_presence, file.path(tab_dir, "supp_fig_s4_s5_gene_presence.tsv"), sep = "\t", quote = FALSE)

row_plots <- list()
for (id in dataset_plan$dataset_id) {
  if (!(id %in% names(objects))) {
    row_plots[[id]] <- ggplot() +
      annotate("text", x = 0, y = 0, label = paste(id, "not loaded"), size = 5) +
      theme_void()
    next
  }
  genes <- strsplit(dataset_plan[dataset_id == id, genes], ",", fixed = TRUE)[[1]]
  label <- dataset_plan[dataset_id == id, dataset_label]
  row_plots[[id]] <- plot_spatial_features(objects[[id]], label, genes)
}

combined <- wrap_plots(row_plots, ncol = 1) +
  plot_annotation(
    title = "Supplementary Fig. S4/S5 spatial feature plots",
    subtitle = "S4 and S5 source PDFs have the same panel content and caption text except the figure number.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 9, color = "#555555")
    )
  )

png_path <- file.path(fig_dir, "supp_fig_s4_s5_spatial_featureplots.png")
pdf_path <- file.path(fig_dir, "supp_fig_s4_s5_spatial_featureplots.pdf")
ggsave(png_path, combined, width = 10.5, height = 12.5, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(pdf_path, combined, width = 10.5, height = 12.5, bg = "white", limitsize = FALSE)
write_log("save_figure", "complete", png_path)

write_log("done", "complete", paste("datasets_loaded=", length(objects)))
