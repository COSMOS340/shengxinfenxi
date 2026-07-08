options(stringsAsFactors = FALSE)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (tolower(Sys.getenv("CODEX_AUTO_INSTALL_R_PACKAGES", "true")) %in% c("1", "true", "yes")) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package missing: ", pkg)
  }
}

for (pkg in c("Matrix", "data.table", "ggplot2")) {
  install_if_missing(pkg)
}

library(Matrix)
library(data.table)
library(ggplot2)

find_repo_root <- function(start) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, ".git"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Cannot find repository root")
    }
    current <- parent
  }
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  script_path <- normalizePath(
    "bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_gse225857_nonimmune_spatial_mapping/scripts/pc_gse225857_nonimmune_spatial_mapping.R",
    winslash = "/",
    mustWork = TRUE
  )
}

request_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- find_repo_root(request_dir)
project_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche")
raw_tar <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "raw_core_geo", "GSE225857", "GSE225857_RAW.tar")
upload_dir <- file.path(project_dir, "desktop_exchange", "uploads", "20260708_gse225857_nonimmune_spatial_mapping")
dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  data.table::fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

read_gz_table <- function(path, header = TRUE) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, header = header, sep = "\t", check.names = FALSE))
}

strip_quotes <- function(x) {
  gsub('^"|"$', "", x)
}

metric_summary <- function(x) {
  x <- as.numeric(x)
  data.table(
    n = sum(!is.na(x)),
    min = suppressWarnings(min(x, na.rm = TRUE)),
    q1 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)),
    median = stats::median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    q3 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)),
    max = suppressWarnings(max(x, na.rm = TRUE))
  )
}

marker_programs <- data.table(
  program = c(
    rep("fibroblast_ecm", 6),
    rep("caf_activation", 6),
    rep("f01_prelp", 5),
    rep("f02_mcam", 6),
    rep("f03_cxcl14", 5),
    rep("f04_c3", 4),
    rep("f05_coch", 3),
    rep("f06_cycling", 3),
    rep("pericyte_smc", 5),
    rep("endothelial", 3),
    rep("epithelial_tumor", 4),
    rep("immune_context", 3)
  ),
  gene_symbol = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PRELP",
    "FAP", "ACTA2", "POSTN", "THBS2", "MMP2", "MMP11",
    "PRELP", "DCN", "LUM", "COL1A1", "COL1A2",
    "MCAM", "RGS5", "PDGFRB", "CSPG4", "ACTA2", "TAGLN",
    "CXCL14", "CFD", "C3", "CXCL12", "SFRP1",
    "C3", "CFD", "CXCL12", "CCL2",
    "COCH", "COL14A1", "PI16",
    "MKI67", "TOP2A", "CENPF",
    "RGS5", "PDGFRB", "MCAM", "CSPG4", "ACTA2",
    "PECAM1", "VWF", "KDR",
    "EPCAM", "KRT8", "KRT18", "KRT19",
    "PTPRC", "CD3D", "LYZ"
  )
)
target_genes <- sort(unique(marker_programs$gene_symbol))

sample_table <- data.table(
  gsm = c("GSM7058756", "GSM7058757", "GSM7058758", "GSM7058759", "GSM7058760", "GSM7058761"),
  spatial_sample = c("C1", "C2", "C3", "C4", "L1", "L2"),
  tissue_site = c("colorectal tissue", "colorectal tissue", "colorectal tissue", "colorectal tissue", "liver tissue", "liver tissue"),
  lesion_type = c(
    "primary tumor sample 1",
    "primary tumor sample 2",
    "primary tumor sample 3",
    "primary tumor sample 4",
    "liver metastatic tumor sample 1",
    "liver metastatic tumor sample 2"
  )
)

required_members <- c(
  "GSM7058755_non_immune_counts.txt.gz",
  "GSM7058755_non_immune_meta.txt.gz",
  paste0(sample_table$gsm, "_", sample_table$spatial_sample, ".matrix.mtx.gz"),
  paste0(sample_table$gsm, "_", sample_table$spatial_sample, ".features.tsv.gz"),
  paste0(sample_table$gsm, "_", sample_table$spatial_sample, ".barcodes.tsv.gz"),
  paste0(sample_table$gsm, "_", sample_table$spatial_sample, "_tissue_positions_list.csv.gz")
)
tar_members <- utils::untar(raw_tar, list = TRUE)
missing_members <- setdiff(required_members, tar_members)
if (length(missing_members) > 0) {
  stop("Missing GSE225857 raw tar members: ", paste(missing_members, collapse = ", "))
}

work_dir <- tempfile("gse225857_nonimmune_spatial_")
dir.create(work_dir)
utils::untar(raw_tar, files = required_members, exdir = work_dir)

meta_path <- file.path(work_dir, "GSM7058755_non_immune_meta.txt.gz")
nonimmune_meta <- read_gz_table(meta_path, header = TRUE)
required_meta_fields <- c("cluster", "organs", "patients", "patients_organ", "nCount_RNA", "predicted.doublet", "doublet")
missing_meta_fields <- setdiff(required_meta_fields, names(nonimmune_meta))
if (length(missing_meta_fields) > 0) {
  stop("Nonimmune metadata missing fields: ", paste(missing_meta_fields, collapse = ", "))
}
cell_id_col <- names(nonimmune_meta)[1]
nonimmune_meta[, cell_id := as.character(get(cell_id_col))]
nonimmune_meta[, cluster := as.character(cluster)]
nonimmune_meta[, organs := as.character(organs)]
nonimmune_meta[, patients := as.character(patients)]
nonimmune_meta[, patients_organ := as.character(patients_organ)]

cluster_counts <- nonimmune_meta[, .N, by = .(cluster, organs)]
cluster_counts[, total_cluster_cells := sum(N), by = cluster]
cluster_counts[, organ_fraction_in_cluster := N / total_cluster_cells]
setorder(cluster_counts, cluster, organs)
write_tsv(cluster_counts, file.path(upload_dir, "gse225857_nonimmune_cluster_counts.tsv"))

fibroblast_clusters <- c(
  "F01_fibroblast_PRELP",
  "F02_fibrblast_MCAM",
  "F03_fibroblast_CXCL14",
  "F04_fibroblast_C3",
  "F05_fibroblast_COCH",
  "F06_cycling_MKI67"
)
fib_counts <- nonimmune_meta[cluster %in% fibroblast_clusters, .N, by = .(cluster, organs, patients_organ)]
setorder(fib_counts, cluster, organs, patients_organ)
write_tsv(fib_counts, file.path(upload_dir, "gse225857_nonimmune_fibroblast_cluster_counts_by_organ.tsv"))

stream_marker_counts <- function(count_path, wanted_genes) {
  con <- gzfile(count_path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1)
  if (length(header) != 1) {
    stop("Cannot read count header")
  }
  header_parts <- strsplit(header, "\t", fixed = TRUE)[[1]]
  cell_ids <- strip_quotes(header_parts[-1])
  rows <- list()
  while (length(line <- readLines(con, n = 1)) > 0) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2) {
      next
    }
    gene <- strip_quotes(parts[1])
    if (gene %in% wanted_genes) {
      values <- suppressWarnings(as.numeric(strip_quotes(parts[-1])))
      if (length(values) != length(cell_ids)) {
        stop("Count row length does not match header for gene: ", gene)
      }
      rows[[gene]] <- values
    }
  }
  if (length(rows) == 0) {
    mat <- matrix(numeric(), nrow = 0, ncol = length(cell_ids))
  } else {
    mat <- do.call(rbind, rows)
  }
  colnames(mat) <- cell_ids
  mat
}

nonimmune_marker_counts <- stream_marker_counts(file.path(work_dir, "GSM7058755_non_immune_counts.txt.gz"), target_genes)
nonimmune_gene_coverage <- data.table(
  gene_symbol = target_genes,
  present = target_genes %in% rownames(nonimmune_marker_counts)
)
nonimmune_gene_coverage <- merge(marker_programs, nonimmune_gene_coverage, by = "gene_symbol", all.x = TRUE)
setorder(nonimmune_gene_coverage, program, gene_symbol)
write_tsv(nonimmune_gene_coverage, file.path(upload_dir, "gse225857_nonimmune_marker_gene_coverage.tsv"))

common_cells <- intersect(colnames(nonimmune_marker_counts), nonimmune_meta$cell_id)
if (length(common_cells) == 0) {
  stop("No overlap between nonimmune count columns and metadata cell IDs")
}
nonimmune_marker_counts <- nonimmune_marker_counts[, common_cells, drop = FALSE]
nonimmune_meta_aligned <- nonimmune_meta[match(common_cells, cell_id)]
norm_factor <- as.numeric(nonimmune_meta_aligned$nCount_RNA)
norm_factor[norm_factor <= 0 | is.na(norm_factor)] <- NA_real_
nonimmune_marker_log <- log1p(t(t(nonimmune_marker_counts) / norm_factor * 10000))

expr_rows <- list()
for (cl in sort(unique(nonimmune_meta_aligned$cluster))) {
  idx <- which(nonimmune_meta_aligned$cluster == cl)
  for (gene in rownames(nonimmune_marker_log)) {
    values <- nonimmune_marker_log[gene, idx]
    raw_values <- nonimmune_marker_counts[gene, idx]
    expr_rows[[length(expr_rows) + 1]] <- data.table(
      cluster = cl,
      gene_symbol = gene,
      mean_log_norm = mean(values, na.rm = TRUE),
      pct_expr = mean(raw_values > 0, na.rm = TRUE)
    )
  }
}
nonimmune_expr <- rbindlist(expr_rows)
nonimmune_expr <- merge(nonimmune_expr, unique(marker_programs), by = "gene_symbol", allow.cartesian = TRUE)
setorder(nonimmune_expr, cluster, program, gene_symbol)
write_tsv(nonimmune_expr, file.path(upload_dir, "gse225857_nonimmune_marker_expression_by_cluster.tsv"))

plot_expr <- nonimmune_expr[cluster %in% c(fibroblast_clusters, "E01_endothelial_SELP", "E02_endothelial_DLL4", "E03_endothelial_NOTCH3", "Tu01_AREG", "Tu02_DEFA5", "Tu03_SRRM2")]
plot_expr[, cluster := factor(cluster, levels = rev(unique(plot_expr$cluster)))]
p_dot <- ggplot(plot_expr, aes(x = gene_symbol, y = cluster)) +
  geom_point(aes(size = pct_expr, color = mean_log_norm)) +
  scale_color_gradient(low = "#D8D8D8", high = "#2B4C7E") +
  scale_size(range = c(0.2, 5.5), limits = c(0, 1)) +
  facet_grid(. ~ program, scales = "free_x", space = "free_x") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    strip.text.x = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  ) +
  labs(title = "GSE225857 nonimmune marker expression by author cluster", x = NULL, y = NULL, color = "Mean log norm", size = "Fraction")
ggsave(file.path(upload_dir, "gse225857_nonimmune_marker_dotplot.png"), p_dot, width = 14, height = 6.8, dpi = 220, bg = "white")

safe_feature_symbols <- function(feature_df) {
  if (ncol(feature_df) >= 2) {
    return(as.character(feature_df[[2]]))
  }
  as.character(feature_df[[1]])
}

score_programs <- function(log_mat, programs) {
  scores <- list()
  for (program_name in sort(unique(programs$program))) {
    genes <- intersect(programs[program == program_name, gene_symbol], rownames(log_mat))
    if (length(genes) == 0) {
      scores[[program_name]] <- rep(NA_real_, ncol(log_mat))
    } else {
      scores[[program_name]] <- Matrix::colMeans(log_mat[genes, , drop = FALSE])
    }
  }
  as.data.table(scores)
}

read_positions <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  pos <- as.data.table(read.csv(con, header = FALSE))
  if (ncol(pos) < 6) {
    stop("Unexpected tissue position column count: ", path)
  }
  setnames(pos, names(pos)[1:6], c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"))
  pos[, barcode := as.character(barcode)]
  pos
}

spatial_scores <- list()
spatial_summary_rows <- list()
spatial_coverage_rows <- list()
for (i in seq_len(nrow(sample_table))) {
  gsm <- sample_table$gsm[i]
  spatial_sample <- sample_table$spatial_sample[i]
  prefix <- paste0(gsm, "_", spatial_sample)
  matrix_path <- file.path(work_dir, paste0(prefix, ".matrix.mtx.gz"))
  feature_path <- file.path(work_dir, paste0(prefix, ".features.tsv.gz"))
  barcode_path <- file.path(work_dir, paste0(prefix, ".barcodes.tsv.gz"))
  position_path <- file.path(work_dir, paste0(prefix, "_tissue_positions_list.csv.gz"))
  mat <- Matrix::readMM(gzfile(matrix_path, open = "rt"))
  features <- read_gz_table(feature_path, header = FALSE)
  barcodes <- readLines(gzfile(barcode_path, open = "rt"))
  colnames(mat) <- barcodes
  symbols <- safe_feature_symbols(features)
  rownames(mat) <- make.unique(symbols)
  target_present <- target_genes[target_genes %in% symbols]
  spatial_coverage_rows[[length(spatial_coverage_rows) + 1]] <- data.table(
    spatial_sample = spatial_sample,
    tissue_site = sample_table$tissue_site[i],
    lesion_type = sample_table$lesion_type[i],
    gene_symbol = target_genes,
    present = target_genes %in% symbols
  )
  keep_rows <- rownames(mat) %in% target_present
  marker_counts <- mat[keep_rows, , drop = FALSE]
  marker_symbols <- sub("[.][0-9]+$", "", rownames(marker_counts))
  if (nrow(marker_counts) > 0 && any(duplicated(marker_symbols))) {
    collapsed <- rowsum(as.matrix(marker_counts), group = marker_symbols)
    marker_counts <- Matrix::Matrix(collapsed, sparse = TRUE)
  } else {
    rownames(marker_counts) <- marker_symbols
  }
  total_counts <- Matrix::colSums(mat)
  total_counts[total_counts <= 0] <- NA_real_
  marker_log <- log1p(t(t(marker_counts) / total_counts * 10000))
  scores <- score_programs(marker_log, marker_programs)
  scores[, barcode := colnames(mat)]
  positions <- read_positions(position_path)
  scores <- merge(scores, positions, by = "barcode", all.x = TRUE)
  scores[, `:=`(
    gsm = gsm,
    spatial_sample = spatial_sample,
    tissue_site = sample_table$tissue_site[i],
    lesion_type = sample_table$lesion_type[i],
    nCount_spatial = as.numeric(total_counts),
    nFeature_spatial = as.numeric(Matrix::colSums(mat > 0))
  )]
  spatial_scores[[length(spatial_scores) + 1]] <- scores
  spatial_summary_rows[[length(spatial_summary_rows) + 1]] <- data.table(
    gsm = gsm,
    spatial_sample = spatial_sample,
    tissue_site = sample_table$tissue_site[i],
    lesion_type = sample_table$lesion_type[i],
    features = nrow(mat),
    spots = ncol(mat),
    nonzero = length(mat@x),
    marker_genes_present = length(target_present)
  )
  rm(mat, marker_counts, marker_log, scores)
  gc()
}

spatial_score_dt <- rbindlist(spatial_scores, fill = TRUE)
program_cols <- sort(unique(marker_programs$program))
setcolorder(
  spatial_score_dt,
  c(
    "gsm", "spatial_sample", "tissue_site", "lesion_type", "barcode",
    "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres",
    "nCount_spatial", "nFeature_spatial", program_cols
  )
)
write_tsv(spatial_score_dt, file.path(upload_dir, "gse225857_spatial_program_scores_wide.tsv"))
write_tsv(rbindlist(spatial_summary_rows), file.path(upload_dir, "gse225857_spatial_sample_summary.tsv"))
spatial_gene_coverage <- merge(marker_programs, rbindlist(spatial_coverage_rows), by = "gene_symbol", allow.cartesian = TRUE)
setorder(spatial_gene_coverage, spatial_sample, program, gene_symbol)
write_tsv(spatial_gene_coverage, file.path(upload_dir, "gse225857_spatial_marker_gene_coverage.tsv"))

summary_long <- melt(
  spatial_score_dt,
  id.vars = c("spatial_sample", "tissue_site", "lesion_type"),
  measure.vars = program_cols,
  variable.name = "program",
  value.name = "score"
)
program_summary <- summary_long[, metric_summary(score), by = .(spatial_sample, tissue_site, lesion_type, program)]
setorder(program_summary, spatial_sample, program)
write_tsv(program_summary, file.path(upload_dir, "gse225857_spatial_program_score_summary.tsv"))

heat <- program_summary[, .(spatial_sample, program, mean)]
p_heat <- ggplot(heat, aes(x = spatial_sample, y = program, fill = mean)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#F2F2F2", high = "#7A2E2E") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 8)) +
  labs(title = "GSE225857 spatial program mean score", x = NULL, y = NULL, fill = "Mean score")
ggsave(file.path(upload_dir, "gse225857_spatial_program_score_heatmap.png"), p_heat, width = 7, height = 5, dpi = 220, bg = "white")

plot_spatial_score <- function(dt, score_col, out_name, title_text) {
  p <- ggplot(dt, aes(x = pxl_col_in_fullres, y = -pxl_row_in_fullres, color = .data[[score_col]])) +
    geom_point(size = 0.32, alpha = 0.92) +
    facet_wrap(~ spatial_sample, ncol = 3) +
    scale_color_gradient(low = "#E6E6E6", high = "#8A1538", na.value = "#F5F5F5") +
    coord_equal() +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(size = 9),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    ) +
    labs(title = title_text, color = "Score")
  ggsave(file.path(upload_dir, out_name), p, width = 8.5, height = 6.6, dpi = 240, bg = "white")
}

plot_spatial_score(spatial_score_dt, "fibroblast_ecm", "gse225857_spatial_fibroblast_ecm_score.png", "GSE225857 spatial fibroblast ECM score")
plot_spatial_score(spatial_score_dt, "f02_mcam", "gse225857_spatial_f02_mcam_score.png", "GSE225857 spatial F02 MCAM score")
plot_spatial_score(spatial_score_dt, "f03_cxcl14", "gse225857_spatial_f03_cxcl14_score.png", "GSE225857 spatial F03 CXCL14 score")
plot_spatial_score(spatial_score_dt, "epithelial_tumor", "gse225857_spatial_epithelial_tumor_score.png", "GSE225857 spatial epithelial tumor score")

status_text <- c(
  "# GSE225857 Nonimmune and Spatial Mapping",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  paste0("Raw tar: `", raw_tar, "`"),
  paste0("Nonimmune metadata cells: ", nrow(nonimmune_meta)),
  paste0("Nonimmune marker genes requested: ", length(target_genes)),
  paste0("Nonimmune marker genes present: ", sum(nonimmune_gene_coverage$present)),
  paste0("Spatial samples parsed: ", nrow(sample_table)),
  paste0("Spatial spot rows written: ", nrow(spatial_score_dt)),
  paste0("Spatial program count: ", length(program_cols)),
  "",
  "This step used author nonimmune clusters and predefined marker programs.",
  "It did not recluster nonimmune cells, perform deconvolution, infer histology regions, or create final cell labels.",
  "No raw files, dense count matrices, spatial images, or RDS objects should be committed to GitHub."
)
writeLines(status_text, file.path(upload_dir, "STATUS.md"))
