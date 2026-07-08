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
  script_path <- normalizePath("bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_stromal_first_pass_review/scripts/pc_stromal_first_pass_review.R", winslash = "/", mustWork = TRUE)
}

request_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- find_repo_root(request_dir)
project_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche")
working_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_working_qc_sparse_objects_20260708")
score_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_marker_program_prepass_20260708")
review_object_dir <- file.path(project_dir, "desktop_exchange", "requests", "20260708_core_geo_download", "derived_stromal_first_pass_review_20260708")
upload_dir <- file.path(project_dir, "desktop_exchange", "uploads", "20260708_stromal_first_pass_review")
dir.create(review_object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

file_size_or_na <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  as.numeric(file.info(path)$size)
}

safe_feature_symbols <- function(feature_df) {
  if (ncol(feature_df) >= 2) {
    return(as.character(feature_df[[2]]))
  }
  as.character(feature_df[[1]])
}

recycle_field <- function(x, n, field_name) {
  if (is.null(x)) {
    stop("Required field missing: ", field_name)
  }
  if (length(x) == n) {
    return(as.character(x))
  }
  if (length(x) == 1) {
    return(rep(as.character(x), n))
  }
  stop("Field length mismatch for ", field_name, ": expected 1 or ", n, ", got ", length(x))
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

review_class_colors <- c(
  caf_ecm_review = "#4C78A8",
  pericyte_review = "#72B7B2",
  endothelial_review = "#54A24B",
  epithelial_review = "#B279A2",
  immune_review = "#9D755D",
  unmapped_review = "#BAB0AC"
)

program_to_review_class <- c(
  fibroblast_ecm = "caf_ecm_review",
  caf_activation = "caf_ecm_review",
  inflammatory_caf = "caf_ecm_review",
  pericyte_like = "pericyte_review",
  endothelial = "endothelial_review",
  epithelial_tumor = "epithelial_review",
  immune_context = "immune_review"
)

marker_sets <- list(
  fibroblast_ecm = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PRELP"),
  caf_activation = c("FAP", "ACTA2", "POSTN", "THBS2", "MMP2", "MMP11"),
  inflammatory_caf = c("CXCL14", "CFD", "C3", "IL6", "CXCL12", "SFRP1"),
  pericyte_like = c("MCAM", "RGS5", "PDGFRB", "CSPG4"),
  endothelial = c("PECAM1", "VWF", "KDR"),
  epithelial_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  immune_context = c("PTPRC", "CD3D", "LYZ"),
  cycling = c("MKI67", "TOP2A")
)

marker_table <- do.call(rbind, lapply(names(marker_sets), function(marker_group) {
  data.frame(marker_group = marker_group, gene_symbol = marker_sets[[marker_group]])
}))

plot_stacked_fraction <- function(tab, group_col, class_col, count_col, path, title_text, width = 2200, height = 1100) {
  if (nrow(tab) == 0) {
    return(invisible(FALSE))
  }
  groups <- sort(unique(as.character(tab[[group_col]])))
  classes <- names(review_class_colors)
  mat <- matrix(0, nrow = length(classes), ncol = length(groups), dimnames = list(classes, groups))
  for (i in seq_len(nrow(tab))) {
    group_value <- as.character(tab[[group_col]][i])
    class_value <- as.character(tab[[class_col]][i])
    if (group_value %in% colnames(mat) && class_value %in% rownames(mat)) {
      mat[class_value, group_value] <- mat[class_value, group_value] + as.numeric(tab[[count_col]][i])
    }
  }
  denom <- colSums(mat)
  frac <- sweep(mat, 2, ifelse(denom > 0, denom, 1), "/")
  png(path, width = width, height = height, res = 150)
  par(mar = c(10, 5, 3, 9), xpd = FALSE)
  barplot(frac, col = review_class_colors[rownames(frac)], border = NA, las = 2, ylab = "Fraction", main = title_text, ylim = c(0, 1))
  par(xpd = TRUE)
  legend("topright", inset = c(-0.13, 0), legend = rownames(frac), fill = review_class_colors[rownames(frac)], bty = "n", cex = 0.75)
  dev.off()
  invisible(TRUE)
}

compute_gene_expression_rows <- function(counts, symbols, group_frame, marker_table) {
  lib <- as.numeric(Matrix::colSums(counts))
  lib[lib <= 0] <- 1
  group_cols <- c("dataset_id", "object_label", "sample_label", "tissue_site", "review_class", "top_marker_program")
  groups <- unique(group_frame[, group_cols, drop = FALSE])
  rows <- list()
  for (gi in seq_len(nrow(groups))) {
    idx <- which(
      group_frame$dataset_id == groups$dataset_id[gi] &
        group_frame$object_label == groups$object_label[gi] &
        group_frame$sample_label == groups$sample_label[gi] &
        group_frame$tissue_site == groups$tissue_site[gi] &
        group_frame$review_class == groups$review_class[gi] &
        group_frame$top_marker_program == groups$top_marker_program[gi]
    )
    for (mi in seq_len(nrow(marker_table))) {
      gene <- marker_table$gene_symbol[mi]
      gene_idx <- which(symbols == gene)
      if (length(gene_idx) == 0) {
        rows[[length(rows) + 1]] <- data.frame(
          groups[gi, ],
          marker_group = marker_table$marker_group[mi],
          gene_symbol = gene,
          present = FALSE,
          cell_count = length(idx),
          pct_expr = NA_real_,
          mean_log_norm = NA_real_
        )
      } else {
        raw_counts <- as.numeric(Matrix::colSums(counts[gene_idx, idx, drop = FALSE]))
        norm_expr <- log1p(raw_counts / lib[idx] * 10000)
        rows[[length(rows) + 1]] <- data.frame(
          groups[gi, ],
          marker_group = marker_table$marker_group[mi],
          gene_symbol = gene,
          present = TRUE,
          cell_count = length(idx),
          pct_expr = mean(raw_counts > 0),
          mean_log_norm = mean(norm_expr)
        )
      }
    }
  }
  do.call(rbind, rows)
}

weighted_expression_by_dataset_class <- function(expr) {
  expr_present <- expr[expr$present & !is.na(expr$mean_log_norm) & !is.na(expr$pct_expr), , drop = FALSE]
  if (nrow(expr_present) == 0) {
    return(data.frame())
  }
  expr_present$weighted_mean_log_norm <- expr_present$mean_log_norm * expr_present$cell_count
  expr_present$weighted_pct_expr <- expr_present$pct_expr * expr_present$cell_count
  agg <- aggregate(
    cbind(weighted_mean_log_norm, weighted_pct_expr, cell_count) ~ dataset_id + review_class + marker_group + gene_symbol,
    data = expr_present,
    FUN = sum
  )
  agg$mean_log_norm <- agg$weighted_mean_log_norm / agg$cell_count
  agg$pct_expr <- agg$weighted_pct_expr / agg$cell_count
  agg <- agg[, c("dataset_id", "review_class", "marker_group", "gene_symbol", "cell_count", "pct_expr", "mean_log_norm")]
  agg[order(agg$dataset_id, agg$review_class, agg$marker_group, agg$gene_symbol), ]
}

plot_expression_heatmap <- function(tab, path) {
  if (nrow(tab) == 0) {
    return(invisible(FALSE))
  }
  gene_order <- unique(marker_table$gene_symbol)
  col_order <- sort(unique(paste(tab$dataset_id, tab$review_class, sep = ":")))
  mat <- matrix(NA_real_, nrow = length(gene_order), ncol = length(col_order), dimnames = list(gene_order, col_order))
  for (i in seq_len(nrow(tab))) {
    col_name <- paste(tab$dataset_id[i], tab$review_class[i], sep = ":")
    mat[tab$gene_symbol[i], col_name] <- tab$mean_log_norm[i]
  }
  mat[is.na(mat)] <- 0
  png(path, width = 1800, height = 1700, res = 150)
  par(mar = c(11, 9, 4, 2))
  z <- t(mat[nrow(mat):1, , drop = FALSE])
  image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = z,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "Dataset review-class marker expression",
    col = colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(100)
  )
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.65)
  axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.75)
  box()
  dev.off()
  invisible(TRUE)
}

score_paths <- sort(list.files(score_object_dir, pattern = "_marker_program_prepass[.]rds$", full.names = TRUE))
if (length(score_paths) == 0) {
  stop("No marker-program prepass RDS objects found in ", score_object_dir)
}

count_rows <- list()
margin_rows <- list()
rank_rows <- list()
expression_rows <- list()
manifest_rows <- list()

for (score_path in score_paths) {
  score_obj <- readRDS(score_path)
  required_fields <- c("dataset_id", "object_label", "sample_label", "tissue_site", "score_matrix", "top_marker_program")
  missing_fields <- setdiff(required_fields, names(score_obj))
  if (length(missing_fields) > 0) {
    stop("Score object is missing fields: ", paste(missing_fields, collapse = ", "), " in ", score_path)
  }

  score_mat <- as.matrix(score_obj$score_matrix)
  n_cells <- nrow(score_mat)
  top_marker_program <- recycle_field(score_obj$top_marker_program, n_cells, "top_marker_program")
  sample_label <- recycle_field(score_obj$sample_label, n_cells, "sample_label")
  tissue_site <- recycle_field(score_obj$tissue_site, n_cells, "tissue_site")
  dataset_id <- as.character(score_obj$dataset_id)
  object_label <- as.character(score_obj$object_label)

  if (!all(top_marker_program %in% names(program_to_review_class))) {
    unknown_programs <- sort(unique(top_marker_program[!top_marker_program %in% names(program_to_review_class)]))
    stop("Unknown top marker programs in ", object_label, ": ", paste(unknown_programs, collapse = ", "))
  }

  sorted_scores <- t(apply(score_mat, 1, sort, decreasing = TRUE))
  top_score <- sorted_scores[, 1]
  second_score <- sorted_scores[, 2]
  score_margin <- top_score - second_score
  review_class <- unname(program_to_review_class[top_marker_program])
  review_class[is.na(review_class)] <- "unmapped_review"

  group_frame <- data.frame(
    dataset_id = dataset_id,
    object_label = object_label,
    sample_label = sample_label,
    tissue_site = tissue_site,
    review_class = review_class,
    top_marker_program = top_marker_program,
    top_score = top_score,
    second_score = second_score,
    score_margin = score_margin
  )

  working_path <- file.path(working_object_dir, paste0(object_label, ".rds"))
  if (!file.exists(working_path)) {
    stop("Working object not found for ", object_label, ": ", working_path)
  }
  working_obj <- readRDS(working_path)
  if (!all(c("counts", "features") %in% names(working_obj))) {
    stop("Working object is missing counts or features: ", working_path)
  }
  counts <- working_obj$counts
  if (ncol(counts) != n_cells) {
    stop("Cell count mismatch for ", object_label, ": score rows ", n_cells, ", count columns ", ncol(counts))
  }
  symbols <- safe_feature_symbols(working_obj$features)

  unique_groups <- unique(group_frame[, c("dataset_id", "object_label", "sample_label", "tissue_site"), drop = FALSE])
  for (gi in seq_len(nrow(unique_groups))) {
    group_idx <- which(
      group_frame$dataset_id == unique_groups$dataset_id[gi] &
        group_frame$object_label == unique_groups$object_label[gi] &
        group_frame$sample_label == unique_groups$sample_label[gi] &
        group_frame$tissue_site == unique_groups$tissue_site[gi]
    )
    total_cells <- length(group_idx)
    class_table <- table(group_frame$review_class[group_idx])
    caf_ecm_count <- if ("caf_ecm_review" %in% names(class_table)) as.integer(class_table[["caf_ecm_review"]]) else 0L
    pericyte_count <- if ("pericyte_review" %in% names(class_table)) as.integer(class_table[["pericyte_review"]]) else 0L
    endothelial_count <- if ("endothelial_review" %in% names(class_table)) as.integer(class_table[["endothelial_review"]]) else 0L
    epithelial_count <- if ("epithelial_review" %in% names(class_table)) as.integer(class_table[["epithelial_review"]]) else 0L
    immune_count <- if ("immune_review" %in% names(class_table)) as.integer(class_table[["immune_review"]]) else 0L

    for (class_name in names(class_table)) {
      class_count <- as.integer(class_table[[class_name]])
      count_rows[[length(count_rows) + 1]] <- data.frame(
        unique_groups[gi, ],
        review_class = class_name,
        cell_count = class_count,
        total_cells = total_cells,
        fraction = class_count / total_cells
      )
      class_idx <- group_idx[group_frame$review_class[group_idx] == class_name]
      margin_summary <- summary_stats(group_frame$score_margin[class_idx])
      top_score_summary <- summary_stats(group_frame$top_score[class_idx])
      names(margin_summary) <- paste0("score_margin_", names(margin_summary))
      names(top_score_summary) <- paste0("top_score_", names(top_score_summary))
      margin_rows[[length(margin_rows) + 1]] <- cbind(
        unique_groups[gi, ],
        review_class = class_name,
        margin_summary,
        top_score_summary
      )
    }

    rank_rows[[length(rank_rows) + 1]] <- data.frame(
      unique_groups[gi, ],
      total_cells = total_cells,
      caf_ecm_review_cells = caf_ecm_count,
      caf_ecm_review_fraction = caf_ecm_count / total_cells,
      pericyte_review_cells = pericyte_count,
      endothelial_review_cells = endothelial_count,
      epithelial_review_cells = epithelial_count,
      immune_review_cells = immune_count
    )
  }

  expression_rows[[length(expression_rows) + 1]] <- compute_gene_expression_rows(counts, symbols, group_frame, marker_table)

  review_object_path <- file.path(review_object_dir, paste0(object_label, "_stromal_first_pass_review.rds"))
  saveRDS(
    list(
      dataset_id = dataset_id,
      object_label = object_label,
      review_classes = program_to_review_class,
      review_table = group_frame
    ),
    review_object_path
  )
  manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
    dataset_id = dataset_id,
    object_label = object_label,
    object_path = review_object_path,
    file_size_bytes = file_size_or_na(review_object_path),
    cells = n_cells
  )

  rm(score_obj, score_mat, sorted_scores, working_obj, counts)
  gc()
}

review_counts <- do.call(rbind, count_rows)
margin_summary <- do.call(rbind, margin_rows)
sample_rank <- do.call(rbind, rank_rows)
marker_expression <- do.call(rbind, expression_rows)
object_manifest <- do.call(rbind, manifest_rows)
dataset_expression <- weighted_expression_by_dataset_class(marker_expression)

review_counts <- review_counts[order(review_counts$dataset_id, review_counts$sample_label, review_counts$review_class), ]
margin_summary <- margin_summary[order(margin_summary$dataset_id, margin_summary$sample_label, margin_summary$review_class), ]
sample_rank <- sample_rank[order(sample_rank$dataset_id, -sample_rank$caf_ecm_review_fraction, sample_rank$sample_label), ]
marker_expression <- marker_expression[order(marker_expression$dataset_id, marker_expression$sample_label, marker_expression$review_class, marker_expression$marker_group, marker_expression$gene_symbol), ]

write_tsv(review_counts, file.path(upload_dir, "review_class_counts_by_group.tsv"))
write_tsv(margin_summary, file.path(upload_dir, "review_class_score_margin_summary.tsv"))
write_tsv(sample_rank, file.path(upload_dir, "sample_caf_ecm_review_rank.tsv"))
write_tsv(marker_expression, file.path(upload_dir, "marker_expression_by_review_class.tsv"))
write_tsv(dataset_expression, file.path(upload_dir, "dataset_review_class_marker_expression.tsv"))
write_tsv(object_manifest, file.path(upload_dir, "stromal_first_pass_object_manifest.tsv"))

gse178_counts <- review_counts[review_counts$dataset_id == "GSE178318", , drop = FALSE]
plot_stacked_fraction(
  gse178_counts,
  group_col = "sample_label",
  class_col = "review_class",
  count_col = "cell_count",
  path = file.path(upload_dir, "gse178318_review_class_fraction_by_suffix.png"),
  title_text = "GSE178318 review-class fraction by barcode suffix",
  width = 2300,
  height = 1150
)

gse245_counts <- review_counts[review_counts$dataset_id == "GSE245552", , drop = FALSE]
if (nrow(gse245_counts) > 0) {
  tissue_counts <- aggregate(
    gse245_counts$cell_count,
    by = list(tissue_site = gse245_counts$tissue_site, review_class = gse245_counts$review_class),
    FUN = sum
  )
  names(tissue_counts)[3] <- "cell_count"
  plot_stacked_fraction(
    tissue_counts,
    group_col = "tissue_site",
    class_col = "review_class",
    count_col = "cell_count",
    path = file.path(upload_dir, "gse245552_review_class_fraction_by_tissue.png"),
    title_text = "GSE245552 review-class fraction by tissue",
    width = 2200,
    height = 1100
  )
  plot_stacked_fraction(
    gse245_counts,
    group_col = "sample_label",
    class_col = "review_class",
    count_col = "cell_count",
    path = file.path(upload_dir, "gse245552_review_class_fraction_by_sample.png"),
    title_text = "GSE245552 review-class fraction by sample",
    width = 3000,
    height = 1300
  )
}

plot_expression_heatmap(dataset_expression, file.path(upload_dir, "dataset_review_class_marker_expression_heatmap.png"))

status_text <- c(
  "# CRLM CAF/ECM stromal first-pass review",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  paste0("Working object directory: `", working_object_dir, "`"),
  paste0("Marker-program score object directory: `", score_object_dir, "`"),
  paste0("PC-local review object directory: `", review_object_dir, "`"),
  "",
  paste0("Score objects reviewed: ", length(score_paths)),
  paste0("Review count rows: ", nrow(review_counts)),
  paste0("Score-margin summary rows: ", nrow(margin_summary)),
  paste0("Sample CAF/ECM review-rank rows: ", nrow(sample_rank)),
  paste0("Marker expression rows: ", nrow(marker_expression)),
  paste0("Dataset-level marker expression rows: ", nrow(dataset_expression)),
  paste0("PC-local review RDS rows: ", nrow(object_manifest)),
  "",
  "No final annotation, integration, clustering, UMAP, marker testing, or publication figure generation was performed.",
  "PC-local RDS objects should not be committed to GitHub."
)
writeLines(status_text, file.path(upload_dir, "STATUS.md"))
