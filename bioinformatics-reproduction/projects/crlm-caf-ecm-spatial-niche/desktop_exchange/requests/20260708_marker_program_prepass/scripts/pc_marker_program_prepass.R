options(stringsAsFactors = FALSE)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Required R package missing: Matrix")
}

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
  script_path <- normalizePath("bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_marker_program_prepass/scripts/pc_marker_program_prepass.R", winslash = "/", mustWork = TRUE)
}

request_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- find_repo_root(request_dir)
project_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche")
raw_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "raw_core_geo")
working_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_working_qc_sparse_objects_20260708")
score_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_marker_program_prepass_20260708")
upload_dir <- file.path(project_dir, "desktop_exchange", "uploads", "20260708_marker_program_prepass")
dir.create(score_object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

read_gz_table_base <- function(path, header = TRUE) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, header = header, check.names = FALSE, sep = "\t")
}

safe_feature_symbols <- function(feature_df) {
  if (ncol(feature_df) >= 2) {
    return(as.character(feature_df[[2]]))
  }
  as.character(feature_df[[1]])
}

file_size_or_na <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  as.numeric(file.info(path)$size)
}

programs <- list(
  fibroblast_ecm = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PRELP"),
  caf_activation = c("FAP", "ACTA2", "POSTN", "THBS2", "MMP2", "MMP11"),
  inflammatory_caf = c("CXCL14", "CFD", "C3", "IL6", "CXCL12", "SFRP1"),
  pericyte_like = c("MCAM", "RGS5", "PDGFRB"),
  endothelial = c("PECAM1", "VWF", "KDR"),
  epithelial_tumor = c("EPCAM", "KRT8", "KRT18"),
  immune_context = c("PTPRC", "CD3D", "LYZ")
)

sparse_log_norm_scores <- function(counts, symbols, genes) {
  idx <- which(symbols %in% genes)
  present <- symbols[idx]
  if (length(idx) == 0) {
    return(rep(NA_real_, ncol(counts)))
  }
  lib <- as.numeric(Matrix::colSums(counts))
  lib[lib <= 0] <- 1
  sub <- counts[idx, , drop = FALSE]
  sub <- as(sub, "dgCMatrix")
  if (length(sub@x) > 0) {
    sub@x <- log1p(sub@x / rep.int(lib, diff(sub@p)) * 10000)
  }
  as.numeric(Matrix::colSums(sub)) / length(genes)
}

summary_stats <- function(x) {
  x <- as.numeric(x)
  data.frame(
    n = length(x),
    q1 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)),
    median = stats::median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    q3 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE))
  )
}

plot_stacked_fraction <- function(tab, group_col, program_col, count_col, path, title_text) {
  if (nrow(tab) == 0) {
    return(invisible(FALSE))
  }
  groups <- sort(unique(as.character(tab[[group_col]])))
  program_names <- names(programs)
  mat <- matrix(0, nrow = length(program_names), ncol = length(groups), dimnames = list(program_names, groups))
  for (i in seq_len(nrow(tab))) {
    g <- as.character(tab[[group_col]][i])
    p <- as.character(tab[[program_col]][i])
    if (p %in% rownames(mat) && g %in% colnames(mat)) {
      mat[p, g] <- mat[p, g] + as.numeric(tab[[count_col]][i])
    }
  }
  denom <- colSums(mat)
  frac <- sweep(mat, 2, ifelse(denom > 0, denom, 1), "/")
  png(path, width = 2200, height = 1100, res = 150)
  par(mar = c(10, 5, 3, 9), xpd = FALSE)
  cols <- c("#4C78A8", "#F58518", "#E45756", "#72B7B2", "#54A24B", "#B279A2", "#9D755D")
  barplot(frac, col = cols, border = NA, las = 2, ylab = "Fraction", main = title_text, ylim = c(0, 1))
  par(xpd = TRUE)
  legend("topright", inset = c(-0.13, 0), legend = rownames(frac), fill = cols, bty = "n", cex = 0.75)
  dev.off()
  invisible(TRUE)
}

status_lines <- c()

gse225857_tar <- file.path(raw_dir, "GSE225857", "GSE225857_RAW.tar")
gse225857_temp <- tempfile("gse225857_marker_program_")
dir.create(gse225857_temp)
utils::untar(gse225857_tar, files = c("GSM7058755_non_immune_meta.txt.gz"), exdir = gse225857_temp)
nonimmune_meta <- read_gz_table_base(file.path(gse225857_temp, "GSM7058755_non_immune_meta.txt.gz"), header = TRUE)
cluster_values <- as.character(nonimmune_meta[["cluster"]])
source_prefix <- sub("[0-9_].*$", "", cluster_values)
f_prefix <- source_prefix == "F"
f_prefix_table <- aggregate(
  rep(1L, sum(f_prefix)),
  by = list(
    cluster = cluster_values[f_prefix],
    organs = as.character(nonimmune_meta[["organs"]][f_prefix]),
    patients_organ = as.character(nonimmune_meta[["patients_organ"]][f_prefix])
  ),
  FUN = sum
)
names(f_prefix_table)[4] <- "cell_count"
f_prefix_table <- f_prefix_table[order(f_prefix_table$cluster, f_prefix_table$organs, f_prefix_table$patients_organ), ]
write_tsv(f_prefix_table, file.path(upload_dir, "gse225857_f_prefix_source_clusters_by_organ.tsv"))
status_lines <- c(status_lines, sprintf("GSE225857 F-prefix source-cluster rows: %s", nrow(f_prefix_table)))

object_paths <- sort(list.files(working_object_dir, pattern = "_working_qc_sparse_matrix[.]rds$", full.names = TRUE))
if (length(object_paths) == 0) {
  stop("No working QC sparse RDS objects found in ", working_object_dir)
}

coverage_rows <- list()
score_summary_rows <- list()
top_count_rows <- list()
stromal_count_rows <- list()
prepass_manifest_rows <- list()

for (object_path in object_paths) {
  obj <- readRDS(object_path)
  counts <- obj$counts
  features <- obj$features
  metadata <- obj$metadata
  symbols <- safe_feature_symbols(features)
  dataset_id <- as.character(obj$dataset_id)
  object_label <- sub("[.]rds$", "", basename(object_path))

  sample_label <- if ("sample_label" %in% names(metadata)) as.character(metadata$sample_label) else rep("all", nrow(metadata))
  if (!"sample_label" %in% names(metadata) && "barcode_suffix" %in% names(metadata)) {
    sample_label <- as.character(metadata$barcode_suffix)
  }
  tissue_site <- if ("tissue_site" %in% names(metadata)) as.character(metadata$tissue_site) else rep("", nrow(metadata))
  if (all(tissue_site == "") && !is.null(obj$tissue_site)) {
    tissue_site <- rep(as.character(obj$tissue_site), nrow(metadata))
  }

  score_list <- list()
  for (program_name in names(programs)) {
    genes <- programs[[program_name]]
    present_genes <- intersect(genes, symbols)
    coverage_rows[[length(coverage_rows) + 1]] <- data.frame(
      dataset_id = dataset_id,
      object_label = object_label,
      marker_program = program_name,
      requested_genes = paste(genes, collapse = ","),
      present_genes = paste(present_genes, collapse = ","),
      present_gene_count = length(present_genes),
      requested_gene_count = length(genes)
    )
    score_list[[program_name]] <- sparse_log_norm_scores(counts, symbols, genes)
  }

  score_mat <- do.call(cbind, score_list)
  colnames(score_mat) <- names(score_list)
  top_program <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]

  group_df <- data.frame(
    dataset_id = dataset_id,
    object_label = object_label,
    sample_label = sample_label,
    tissue_site = tissue_site,
    top_program = top_program
  )

  unique_groups <- unique(group_df[, c("dataset_id", "object_label", "sample_label", "tissue_site")])
  for (gi in seq_len(nrow(unique_groups))) {
    idx <- which(
      group_df$dataset_id == unique_groups$dataset_id[gi] &
        group_df$object_label == unique_groups$object_label[gi] &
        group_df$sample_label == unique_groups$sample_label[gi] &
        group_df$tissue_site == unique_groups$tissue_site[gi]
    )
    for (program_name in colnames(score_mat)) {
      s <- summary_stats(score_mat[idx, program_name])
      score_summary_rows[[length(score_summary_rows) + 1]] <- cbind(unique_groups[gi, ], marker_program = program_name, s)
    }
    tab <- table(group_df$top_program[idx])
    for (program_name in names(tab)) {
      top_count_rows[[length(top_count_rows) + 1]] <- data.frame(
        unique_groups[gi, ],
        top_marker_program = program_name,
        cell_count = as.integer(tab[[program_name]]),
        total_cells = length(idx),
        fraction = as.integer(tab[[program_name]]) / length(idx)
      )
    }
    stromal_programs <- c("fibroblast_ecm", "caf_activation", "inflammatory_caf", "pericyte_like")
    stromal_count <- sum(group_df$top_program[idx] %in% stromal_programs)
    stromal_count_rows[[length(stromal_count_rows) + 1]] <- data.frame(
      unique_groups[gi, ],
      stromal_review_program_count = stromal_count,
      total_cells = length(idx),
      stromal_review_program_fraction = stromal_count / length(idx)
    )
  }

  prepass_path <- file.path(score_object_dir, paste0(object_label, "_marker_program_prepass.rds"))
  saveRDS(
    list(
      dataset_id = dataset_id,
      object_label = object_label,
      sample_label = sample_label,
      tissue_site = tissue_site,
      score_programs = programs,
      score_matrix = score_mat,
      top_marker_program = top_program
    ),
    prepass_path
  )
  prepass_manifest_rows[[length(prepass_manifest_rows) + 1]] <- data.frame(
    dataset_id = dataset_id,
    object_label = object_label,
    object_path = prepass_path,
    file_size_bytes = file_size_or_na(prepass_path),
    cells = nrow(score_mat),
    marker_programs = ncol(score_mat)
  )

  rm(obj, counts, score_mat)
  gc()
}

coverage <- do.call(rbind, coverage_rows)
score_summary <- do.call(rbind, score_summary_rows)
top_counts <- do.call(rbind, top_count_rows)
stromal_counts <- do.call(rbind, stromal_count_rows)
prepass_manifest <- do.call(rbind, prepass_manifest_rows)

write_tsv(coverage, file.path(upload_dir, "marker_program_gene_coverage.tsv"))
write_tsv(score_summary, file.path(upload_dir, "marker_program_score_summary_by_group.tsv"))
write_tsv(top_counts, file.path(upload_dir, "marker_program_top_counts_by_group.tsv"))
write_tsv(stromal_counts, file.path(upload_dir, "stromal_review_counts_by_group.tsv"))
write_tsv(prepass_manifest, file.path(upload_dir, "marker_program_prepass_object_manifest.tsv"))

gse178_top <- top_counts[top_counts$dataset_id == "GSE178318", , drop = FALSE]
plot_stacked_fraction(
  gse178_top,
  group_col = "sample_label",
  program_col = "top_marker_program",
  count_col = "cell_count",
  path = file.path(upload_dir, "gse178318_marker_program_top_counts_by_suffix.png"),
  title_text = "GSE178318 marker-program prepass by barcode suffix"
)

gse245_top <- top_counts[top_counts$dataset_id == "GSE245552", , drop = FALSE]
if (nrow(gse245_top) > 0) {
  tissue_top <- aggregate(
    gse245_top$cell_count,
    by = list(tissue_site = gse245_top$tissue_site, top_marker_program = gse245_top$top_marker_program),
    FUN = sum
  )
  names(tissue_top)[3] <- "cell_count"
  plot_stacked_fraction(
    tissue_top,
    group_col = "tissue_site",
    program_col = "top_marker_program",
    count_col = "cell_count",
    path = file.path(upload_dir, "gse245552_marker_program_top_counts_by_tissue.png"),
    title_text = "GSE245552 marker-program prepass by tissue"
  )
}

status_lines <- c(
  status_lines,
  sprintf("Working objects scored: %s", length(object_paths)),
  sprintf("Marker program coverage rows: %s", nrow(coverage)),
  sprintf("Top-count rows: %s", nrow(top_counts)),
  sprintf("Stromal review count rows: %s", nrow(stromal_counts)),
  sprintf("PC-local marker-program prepass RDS rows: %s", nrow(prepass_manifest))
)

status_text <- c(
  "# CRLM CAF/ECM marker-program prepass",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  paste0("Working object directory: `", working_object_dir, "`"),
  paste0("PC-local score object directory: `", score_object_dir, "`"),
  "",
  status_lines,
  "",
  "No final cell-type annotation, integration, clustering, marker testing, or publication figure generation was performed.",
  "Per-cell score RDS objects are PC-local and should not be committed to GitHub."
)
writeLines(status_text, file.path(upload_dir, "STATUS.md"))
