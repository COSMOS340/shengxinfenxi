options(stringsAsFactors = FALSE)

cran_repos <- getOption("repos")
if (is.null(cran_repos) || is.na(cran_repos[["CRAN"]]) || cran_repos[["CRAN"]] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

install_enabled <- identical(Sys.getenv("INSTALL_MISSING_R_PACKAGES", unset = "1"), "1")

ensure_package <- function(pkg, source = "cran", required = TRUE) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(TRUE)
  }
  if (install_enabled) {
    message("Installing missing package: ", pkg)
    if (identical(source, "bioc")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  available <- requireNamespace(pkg, quietly = TRUE)
  if (!available && required) {
    stop("Missing required package after install attempt: ", pkg)
  }
  available
}

required_packages <- c("data.table", "Matrix", "Seurat", "ggplot2", "patchwork", "scales")
invisible(vapply(required_packages, ensure_package, logical(1), source = "cran", required = TRUE))

optional_package_status <- data.frame(
  package = c("SingleCellExperiment", "scDblFinder", "BiocParallel", "FNN", "harmony"),
  source = c("bioc", "bioc", "bioc", "cran", "cran"),
  available = FALSE,
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(optional_package_status))) {
  optional_package_status$available[[i]] <- ensure_package(
    optional_package_status$package[[i]],
    source = optional_package_status$source[[i]],
    required = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    stop("Run this script with Rscript so --file is available.")
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

safe_write_png <- function(plot, path, width, height, dpi = 300) {
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

entropy_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  p <- as.numeric(table(x)) / length(x)
  -sum(p * log(p))
}

norm_entropy_value <- function(x) {
  x <- x[!is.na(x)]
  label_count <- length(unique(x))
  if (length(x) == 0 || label_count <= 1) {
    return(0)
  }
  entropy_value(x) / log(label_count)
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

get_data_layer <- function(obj) {
  args <- names(formals(Seurat::GetAssayData))
  if ("layer" %in% args) {
    return(Seurat::GetAssayData(obj, assay = DefaultAssay(obj), layer = "data"))
  }
  Seurat::GetAssayData(obj, assay = DefaultAssay(obj), slot = "data")
}

script_path <- get_script_path()
request_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
output_dir <- file.path(request_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

previous_request_dir <- normalizePath(
  file.path(request_dir, "..", "20260706_gse164522_full_mac_mono_analysis"),
  mustWork = FALSE
)
default_rds_path <- file.path(
  previous_request_dir,
  "local_objects_not_for_github",
  "gse164522_full_mac_mono_seurat.rds"
)
rds_path <- Sys.getenv("GSE164522_MAC_MONO_RDS", unset = default_rds_path)
rds_path <- normalizePath(rds_path, mustWork = TRUE)

cluster_col <- Sys.getenv("CLUSTER_COL", unset = "RNA_snn_res.0.6")
marker_set_path <- file.path(request_dir, "inputs", "lineage_marker_sets.tsv")
cleaning_set_path <- file.path(request_dir, "inputs", "gse164522_full_mac_mono_cleaning_sets.tsv")
marker_review_path <- file.path(request_dir, "inputs", "gse164522_full_mac_mono_cluster_marker_boundary_review.tsv")
if (!file.exists(marker_set_path)) {
  stop("Missing marker set input: ", marker_set_path)
}
if (!file.exists(cleaning_set_path)) {
  stop("Missing cleaning set input: ", cleaning_set_path)
}
if (!file.exists(marker_review_path)) {
  stop("Missing marker review input: ", marker_review_path)
}

message("Reading Seurat object: ", rds_path)
obj <- readRDS(rds_path)
meta <- as.data.table(obj[[]], keep.rownames = "expression_cell_id")
if (!cluster_col %in% names(meta)) {
  stop("Expected cluster column not found: ", cluster_col, ". Metadata columns: ", paste(names(meta), collapse = ", "))
}
meta[, cluster := as.character(get(cluster_col))]
Idents(obj) <- cluster_col

cluster_levels <- sort(unique(meta$cluster))
cluster_order_numeric <- suppressWarnings(as.numeric(cluster_levels))
if (!any(is.na(cluster_order_numeric))) {
  cluster_levels <- cluster_levels[order(cluster_order_numeric)]
}
cluster_base <- data.table(cluster = cluster_levels)

cleaning_sets <- fread(cleaning_set_path)
marker_review <- fread(marker_review_path)
for (dt_name in c("cleaning_sets", "marker_review")) {
  dt <- get(dt_name)
  if (!"cluster" %in% names(dt)) {
    stop("Input table lacks cluster column: ", dt_name)
  }
  dt[, cluster := as.character(cluster)]
  assign(dt_name, dt)
}

qc_cols <- intersect(
  c("n_genes", "n_counts", "percent_mito", "percent_mt_calculated", "nFeature_RNA", "nCount_RNA"),
  names(meta)
)
qc_summary_list <- list()
for (col in qc_cols) {
  one <- meta[, .(
    median = median(get(col), na.rm = TRUE),
    mean = mean(get(col), na.rm = TRUE),
    q10 = as.numeric(quantile(get(col), 0.10, na.rm = TRUE)),
    q90 = as.numeric(quantile(get(col), 0.90, na.rm = TRUE))
  ), by = cluster]
  one[, metric := col]
  qc_summary_list[[col]] <- one
}
qc_summary_long <- rbindlist(qc_summary_list, use.names = TRUE, fill = TRUE)
fwrite(qc_summary_long, file.path(output_dir, "gse164522_audit_cluster_qc_summary.tsv"), sep = "\t")

dominance_for <- function(dt, column) {
  if (!column %in% names(dt)) {
    return(NULL)
  }
  tab <- dt[, .N, by = .(cluster, label = get(column))]
  totals <- tab[, .(cluster_cells = sum(N)), by = cluster]
  tab <- merge(tab, totals, by = "cluster", all.x = TRUE)
  tab[, fraction := N / cluster_cells]
  max_rows <- tab[order(cluster, -fraction, label), .SD[1], by = cluster]
  ent <- dt[, .(
    entropy = entropy_value(get(column)),
    norm_entropy = norm_entropy_value(get(column)),
    label_count = uniqueN(get(column), na.rm = TRUE)
  ), by = cluster]
  out <- merge(max_rows[, .(cluster, label, fraction)], ent, by = "cluster", all = TRUE)
  setnames(
    out,
    c("label", "fraction", "entropy", "norm_entropy", "label_count"),
    paste0(c("max_", "max_", "", "norm_", ""), column, c("", "_fraction", "_entropy", "_entropy", "_count"))
  )
  out
}

composition_list <- lapply(c("sample", "patient", "tissue", "celltype_sub", "celltype_major", "expression_group"), function(x) {
  dominance_for(meta, x)
})
composition_summary <- Reduce(function(x, y) merge(x, y, by = "cluster", all = TRUE), c(list(cluster_base), composition_list[!vapply(composition_list, is.null, logical(1))]))
fwrite(composition_summary, file.path(output_dir, "gse164522_audit_cluster_composition_summary.tsv"), sep = "\t")

marker_sets <- fread(marker_set_path)
required_marker_cols <- c("score_set", "gene")
missing_marker_cols <- setdiff(required_marker_cols, names(marker_sets))
if (length(missing_marker_cols) > 0) {
  stop("Marker set columns missing: ", paste(missing_marker_cols, collapse = ", "))
}
marker_sets[, score_set := sanitize_name(score_set)]
gene_universe <- rownames(obj)
marker_coverage <- marker_sets[, .(
  input_genes = .N,
  present_genes = sum(gene %in% gene_universe),
  missing_genes = paste(sort(setdiff(gene, gene_universe)), collapse = ","),
  present_gene_list = paste(sort(intersect(gene, gene_universe)), collapse = ",")
), by = score_set]
fwrite(marker_coverage, file.path(output_dir, "gse164522_audit_marker_set_coverage.tsv"), sep = "\t")

log_data <- get_data_layer(obj)
lineage_score_dt <- data.table(expression_cell_id = colnames(obj), cluster = meta$cluster[match(colnames(obj), meta$expression_cell_id)])
for (set_name in marker_coverage$score_set) {
  genes <- marker_sets[score_set == set_name & gene %in% rownames(log_data), gene]
  score_col <- paste0("score_", set_name)
  if (length(genes) == 0) {
    lineage_score_dt[, (score_col) := NA_real_]
  } else {
    lineage_score_dt[, (score_col) := Matrix::colMeans(log_data[genes, , drop = FALSE])]
  }
  z_col <- paste0(score_col, "_z")
  values <- lineage_score_dt[[score_col]]
  lineage_score_dt[, (z_col) := as.numeric(scale(values))]
}
fwrite(lineage_score_dt, file.path(output_dir, "gse164522_audit_per_cell_lineage_scores.tsv.gz"), sep = "\t")

score_z_cols <- grep("^score_.*_z$", names(lineage_score_dt), value = TRUE)
score_raw_cols <- setdiff(grep("^score_", names(lineage_score_dt), value = TRUE), score_z_cols)
score_summary_raw <- melt(
  lineage_score_dt,
  id.vars = c("expression_cell_id", "cluster"),
  measure.vars = score_raw_cols,
  variable.name = "score_name",
  value.name = "score"
)
score_summary_raw[, score_set := sub("^score_", "", score_name)]
score_summary_z <- melt(
  lineage_score_dt,
  id.vars = c("expression_cell_id", "cluster"),
  measure.vars = score_z_cols,
  variable.name = "score_name",
  value.name = "score_z"
)
score_summary_z[, score_set := sub("^score_", "", sub("_z$", "", score_name))]
score_long <- merge(
  score_summary_raw[, .(expression_cell_id, cluster, score_set, score)],
  score_summary_z[, .(expression_cell_id, cluster, score_set, score_z)],
  by = c("expression_cell_id", "cluster", "score_set"),
  all = TRUE
)
cluster_lineage_summary <- score_long[, .(
  median_score = median(score, na.rm = TRUE),
  mean_score = mean(score, na.rm = TRUE),
  median_score_z = median(score_z, na.rm = TRUE),
  mean_score_z = mean(score_z, na.rm = TRUE)
), by = .(cluster, score_set)]
fwrite(cluster_lineage_summary, file.path(output_dir, "gse164522_audit_cluster_lineage_score_summary.tsv"), sep = "\t")

myeloid_by_cell <- score_long[score_set == "Myeloid_Mac_Mono", .(expression_cell_id, myeloid_score_z = score_z)]
score_long[, score_z_for_rank := fifelse(is.na(score_z), -Inf, score_z)]
non_myeloid_by_cell <- score_long[score_set != "Myeloid_Mac_Mono", .SD[which.max(score_z_for_rank)], by = expression_cell_id]
non_myeloid_by_cell[!is.finite(score_z_for_rank), `:=`(score_z, NA_real_, score_set = NA_character_)]
setnames(non_myeloid_by_cell, c("score_set", "score_z"), c("top_non_myeloid_score_set", "top_non_myeloid_score_z"))
lineage_conflict <- merge(myeloid_by_cell, non_myeloid_by_cell[, .(expression_cell_id, top_non_myeloid_score_set, top_non_myeloid_score_z)], by = "expression_cell_id", all = TRUE)
lineage_conflict <- merge(lineage_conflict, meta[, .(expression_cell_id, cluster)], by = "expression_cell_id", all.x = TRUE)
lineage_conflict[, non_myeloid_exceeds_myeloid := !is.na(top_non_myeloid_score_z) & !is.na(myeloid_score_z) & top_non_myeloid_score_z > myeloid_score_z]
mode_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  names(sort(table(x), decreasing = TRUE))[1]
}
lineage_conflict_summary <- lineage_conflict[, .(
  top_non_myeloid_fraction = mean(non_myeloid_exceeds_myeloid, na.rm = TRUE),
  top_non_myeloid_set_mode = mode_or_na(top_non_myeloid_score_set),
  median_myeloid_score_z = median(myeloid_score_z, na.rm = TRUE),
  median_top_non_myeloid_score_z = median(top_non_myeloid_score_z, na.rm = TRUE)
), by = cluster]

doublet_status <- "scDblFinder not run"
doublet_cluster_summary <- data.table(cluster = cluster_base$cluster)
doublet_per_cell_path <- file.path(output_dir, "gse164522_audit_scDblFinder_per_cell.tsv.gz")
doublet_cluster_path <- file.path(output_dir, "gse164522_audit_scDblFinder_cluster_summary.tsv")
if (all(optional_package_status$available[match(c("SingleCellExperiment", "scDblFinder"), optional_package_status$package)])) {
  doublet_status <- tryCatch({
    set.seed(20260706)
    sce <- Seurat::as.SingleCellExperiment(obj, assay = DefaultAssay(obj))
    cluster_labels <- factor(meta$cluster[match(colnames(sce), meta$expression_cell_id)])
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      sce <- scDblFinder::scDblFinder(sce, clusters = cluster_labels, BPPARAM = BiocParallel::SerialParam())
    } else {
      sce <- scDblFinder::scDblFinder(sce, clusters = cluster_labels)
    }
    cd <- as.data.table(as.data.frame(SummarizedExperiment::colData(sce)), keep.rownames = "expression_cell_id")
    cd[, cluster := as.character(cluster_labels[match(expression_cell_id, colnames(sce))])]
    scdbl_cols <- grep("^scDblFinder", names(cd), value = TRUE)
    fwrite(cd[, c("expression_cell_id", "cluster", scdbl_cols), with = FALSE], doublet_per_cell_path, sep = "\t")
    class_col <- "scDblFinder.class"
    score_col <- "scDblFinder.score"
    if (!all(c(class_col, score_col) %in% names(cd))) {
      stop("scDblFinder output lacks expected columns: scDblFinder.class and scDblFinder.score")
    }
    doublet_cluster_summary <<- cd[, .(
      scDblFinder_cells = .N,
      scDblFinder_doublet_cells = sum(get(class_col) == "doublet", na.rm = TRUE),
      scDblFinder_doublet_fraction = mean(get(class_col) == "doublet", na.rm = TRUE),
      scDblFinder_median_score = median(get(score_col), na.rm = TRUE),
      scDblFinder_q90_score = as.numeric(quantile(get(score_col), 0.90, na.rm = TRUE))
    ), by = cluster]
    fwrite(doublet_cluster_summary, doublet_cluster_path, sep = "\t")
    "scDblFinder completed"
  }, error = function(e) {
    doublet_cluster_summary <<- data.table(cluster = cluster_base$cluster)
    paste0("scDblFinder failed: ", conditionMessage(e))
  })
} else {
  missing_doublet <- optional_package_status$package[
    optional_package_status$package %in% c("SingleCellExperiment", "scDblFinder") & !optional_package_status$available
  ]
  doublet_status <- paste0("scDblFinder skipped because packages are unavailable: ", paste(missing_doublet, collapse = ", "))
}

knn_status <- "FNN not run"
knn_summary <- data.table(cluster = cluster_base$cluster)
if (requireNamespace("FNN", quietly = TRUE) && "pca" %in% Reductions(obj)) {
  knn_status <- tryCatch({
    pca <- Embeddings(obj, "pca")
    dims_use <- seq_len(min(30, ncol(pca)))
    k_use <- min(30, nrow(pca) - 1)
    nn <- FNN::get.knn(pca[, dims_use, drop = FALSE], k = k_use)
    label_meta <- meta[match(rownames(pca), expression_cell_id)]
    local_entropy_for <- function(labels) {
      vapply(seq_len(nrow(nn$nn.index)), function(i) norm_entropy_value(labels[nn$nn.index[i, ]]), numeric(1))
    }
    per_cell_knn <- data.table(
      expression_cell_id = rownames(pca),
      cluster = label_meta$cluster
    )
    for (field in intersect(c("sample", "patient", "tissue", "expression_group"), names(label_meta))) {
      per_cell_knn[, paste0("knn_", field, "_norm_entropy") := local_entropy_for(label_meta[[field]])]
    }
    fwrite(per_cell_knn, file.path(output_dir, "gse164522_audit_per_cell_knn_batch_mixing.tsv.gz"), sep = "\t")
    knn_cols <- grep("^knn_.*_norm_entropy$", names(per_cell_knn), value = TRUE)
    knn_summary <<- per_cell_knn[, lapply(.SD, median, na.rm = TRUE), by = cluster, .SDcols = knn_cols]
    fwrite(knn_summary, file.path(output_dir, "gse164522_audit_cluster_knn_batch_mixing.tsv"), sep = "\t")
    "FNN nearest-neighbor batch mixing completed"
  }, error = function(e) {
    knn_summary <<- data.table(cluster = cluster_base$cluster)
    paste0("FNN nearest-neighbor batch mixing failed: ", conditionMessage(e))
  })
}

integrated <- Reduce(
  function(x, y) merge(x, y, by = "cluster", all = TRUE),
  list(
    cluster_base,
    marker_review,
    cleaning_sets[, setdiff(names(cleaning_sets), intersect(names(cleaning_sets), setdiff(names(marker_review), "cluster"))), with = FALSE],
    dcast(qc_summary_long, cluster ~ metric, value.var = "median"),
    composition_summary,
    lineage_conflict_summary,
    doublet_cluster_summary,
    knn_summary
  )
)

if (!"review_label" %in% names(integrated) && "review_label.x" %in% names(integrated)) {
  setnames(integrated, "review_label.x", "review_label")
}
if (!"initial_action" %in% names(integrated) && "initial_action.x" %in% names(integrated)) {
  setnames(integrated, "initial_action.x", "initial_action")
}

integrated[, sample_dominated_70 := !is.na(max_sample_fraction) & max_sample_fraction >= 0.70]
integrated[, patient_dominated_70 := !is.na(max_patient_fraction) & max_patient_fraction >= 0.70]
integrated[, tissue_dominated_90 := !is.na(max_tissue_fraction) & max_tissue_fraction >= 0.90]
integrated[, author_subtype_dominated_80 := !is.na(max_celltype_sub_fraction) & max_celltype_sub_fraction >= 0.80]
if ("scDblFinder_doublet_fraction" %in% names(integrated)) {
  integrated[, doublet_enriched_15 := !is.na(scDblFinder_doublet_fraction) & scDblFinder_doublet_fraction >= 0.15]
} else {
  integrated[, scDblFinder_doublet_fraction := NA_real_]
  integrated[, scDblFinder_doublet_cells := NA_real_]
  integrated[, scDblFinder_median_score := NA_real_]
  integrated[, scDblFinder_q90_score := NA_real_]
  integrated[, doublet_enriched_15 := NA]
}
integrated[, non_myeloid_score_enriched := !is.na(top_non_myeloid_fraction) & top_non_myeloid_fraction >= 0.50]
integrated[, lineage_boundary_flag := !is.na(review_label) & review_label != "mac_mono_supported"]
integrated[, small_cluster_flag := !is.na(cell_count) & cell_count < 100]
integrated[, doublet_enriched_for_status := fifelse(is.na(doublet_enriched_15), FALSE, doublet_enriched_15)]
integrated[, audit_status := fifelse(
  lineage_boundary_flag & doublet_enriched_for_status,
  "lineage_boundary_and_doublet_enriched_review",
  fifelse(
    lineage_boundary_flag & non_myeloid_score_enriched,
    "lineage_boundary_and_module_score_review",
    fifelse(
      lineage_boundary_flag,
      "lineage_boundary_review",
      fifelse(
        sample_dominated_70 | patient_dominated_70 | tissue_dominated_90,
        "batch_or_sample_review",
        fifelse(small_cluster_flag, "small_cluster_review", "retain_for_now")
      )
    )
  )
)]

fwrite(integrated, file.path(output_dir, "gse164522_audit_cluster_integrated_review.tsv"), sep = "\t")
fwrite(
  integrated[, .(
    cluster,
    cell_count,
    review_label,
    initial_action,
    audit_status,
    lineage_boundary_flag,
    non_myeloid_score_enriched,
    doublet_enriched_15,
    sample_dominated_70,
    patient_dominated_70,
    tissue_dominated_90,
    author_subtype_dominated_80,
    top_non_myeloid_set_mode,
    top_non_myeloid_fraction,
    scDblFinder_doublet_fraction,
    max_sample,
    max_sample_fraction,
    max_patient,
    max_patient_fraction,
    max_tissue,
    max_tissue_fraction,
    max_celltype_sub,
    max_celltype_sub_fraction
  )],
  file.path(output_dir, "gse164522_audit_cluster_review_status.tsv"),
  sep = "\t"
)

umap_dt <- NULL
if ("umap" %in% Reductions(obj)) {
  umap_dt <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "expression_cell_id")
  umap_cols <- setdiff(names(umap_dt), "expression_cell_id")
  if (length(umap_cols) >= 2) {
    setnames(umap_dt, umap_cols[1:2], c("umap_1", "umap_2"))
    umap_dt <- merge(umap_dt, meta[, .(expression_cell_id, cluster, tissue, sample, patient, celltype_sub)], by = "expression_cell_id", all.x = TRUE)
    status_cols <- integrated[, .(cluster, review_label, initial_action, audit_status)]
    umap_dt <- merge(umap_dt, status_cols, by = "cluster", all.x = TRUE)
    base_umap <- function(color_col, title) {
      ggplot(umap_dt, aes(x = umap_1, y = umap_2, color = get(color_col))) +
        geom_point(size = 0.16, alpha = 0.75, stroke = 0) +
        coord_equal() +
        theme_classic(base_size = 9) +
        theme(
          legend.title = element_blank(),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 11)
        ) +
        labs(title = title, x = "UMAP 1", y = "UMAP 2")
    }
    safe_write_png(base_umap("initial_action", "Initial marker action"), file.path(output_dir, "gse164522_audit_umap_initial_action.png"), 6.5, 5.2)
    safe_write_png(base_umap("review_label", "Marker review label"), file.path(output_dir, "gse164522_audit_umap_marker_review.png"), 7.2, 5.2)
    safe_write_png(base_umap("audit_status", "Integrated audit status"), file.path(output_dir, "gse164522_audit_umap_status.png"), 8.5, 5.5)
    safe_write_png(base_umap("tissue", "Tissue"), file.path(output_dir, "gse164522_audit_umap_tissue.png"), 7.0, 5.2)
  }
}

harmony_status <- "Harmony not run"
if (requireNamespace("harmony", quietly = TRUE) && "pca" %in% Reductions(obj) && "sample" %in% names(meta)) {
  harmony_status <- tryCatch({
    obj_harmony <- obj
    pca_harmony <- Embeddings(obj_harmony, "pca")
    dims_harmony <- seq_len(min(30, ncol(pca_harmony)))
    obj_harmony <- harmony::RunHarmony(
      object = obj_harmony,
      group.by.vars = "sample",
      reduction.use = "pca",
      dims.use = dims_harmony,
      assay.use = DefaultAssay(obj_harmony),
      verbose = FALSE
    )
    obj_harmony <- RunUMAP(
      obj_harmony,
      reduction = "harmony",
      dims = dims_harmony,
      reduction.name = "umap_harmony_sample",
      seed.use = 20260706,
      verbose = FALSE
    )
    h_dt <- as.data.table(Embeddings(obj_harmony, "umap_harmony_sample"), keep.rownames = "expression_cell_id")
    h_cols <- setdiff(names(h_dt), "expression_cell_id")
    setnames(h_dt, h_cols[1:2], c("umap_harmony_1", "umap_harmony_2"))
    h_dt <- merge(h_dt, meta[, .(expression_cell_id, cluster, tissue, sample, patient, celltype_sub)], by = "expression_cell_id", all.x = TRUE)
    h_dt <- merge(h_dt, integrated[, .(cluster, review_label, initial_action, audit_status)], by = "cluster", all.x = TRUE)
    fwrite(h_dt, file.path(output_dir, "gse164522_audit_harmony_sample_umap_coordinates.tsv.gz"), sep = "\t")
    harmony_umap <- function(color_col, title) {
      ggplot(h_dt, aes(x = umap_harmony_1, y = umap_harmony_2, color = get(color_col))) +
        geom_point(size = 0.16, alpha = 0.75, stroke = 0) +
        coord_equal() +
        theme_classic(base_size = 9) +
        theme(
          legend.title = element_blank(),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 11)
        ) +
        labs(title = title, x = "Harmony UMAP 1", y = "Harmony UMAP 2")
    }
    safe_write_png(harmony_umap("cluster", "Harmony by sample: cluster"), file.path(output_dir, "gse164522_audit_harmony_sample_umap_cluster.png"), 7.0, 5.2)
    safe_write_png(harmony_umap("tissue", "Harmony by sample: tissue"), file.path(output_dir, "gse164522_audit_harmony_sample_umap_tissue.png"), 7.0, 5.2)
    safe_write_png(harmony_umap("audit_status", "Harmony by sample: audit status"), file.path(output_dir, "gse164522_audit_harmony_sample_umap_status.png"), 8.5, 5.5)
    "Harmony by sample completed"
  }, error = function(e) {
    paste0("Harmony by sample failed: ", conditionMessage(e))
  })
}

plot_bar_fraction <- function(dt, value_col, title, ylab, path) {
  if (!value_col %in% names(dt)) {
    return(invisible(NULL))
  }
  p <- ggplot(dt, aes(x = factor(cluster, levels = cluster_levels), y = get(value_col), fill = review_label)) +
    geom_col(width = 0.78) +
    geom_hline(yintercept = 0.70, linetype = "dashed", linewidth = 0.3) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11)
    ) +
    labs(title = title, y = ylab)
  safe_write_png(p, path, 9, 5)
}

plot_bar_fraction(integrated, "max_sample_fraction", "Dominant sample fraction by cluster", "Max sample fraction", file.path(output_dir, "gse164522_audit_cluster_max_sample_fraction.png"))
plot_bar_fraction(integrated, "max_patient_fraction", "Dominant patient fraction by cluster", "Max patient fraction", file.path(output_dir, "gse164522_audit_cluster_max_patient_fraction.png"))
plot_bar_fraction(integrated, "max_celltype_sub_fraction", "Dominant author subtype fraction by cluster", "Max author subtype fraction", file.path(output_dir, "gse164522_audit_cluster_author_subtype_fraction.png"))
if ("scDblFinder_doublet_fraction" %in% names(integrated)) {
  plot_bar_fraction(integrated, "scDblFinder_doublet_fraction", "scDblFinder doublet fraction by cluster", "Doublet fraction", file.path(output_dir, "gse164522_audit_cluster_scDblFinder_doublet_fraction.png"))
}

if ("nFeature_RNA" %in% names(integrated)) {
  p_qc <- ggplot(integrated, aes(x = factor(cluster, levels = cluster_levels), y = nFeature_RNA, fill = review_label)) +
    geom_col(width = 0.78) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11)
    ) +
    labs(title = "Median nFeature_RNA by cluster", y = "Median nFeature_RNA")
  safe_write_png(p_qc, file.path(output_dir, "gse164522_audit_cluster_median_nFeature_RNA.png"), 9, 5)
}

lineage_tile <- copy(cluster_lineage_summary)
lineage_tile[, cluster := factor(cluster, levels = cluster_levels)]
lineage_tile[, score_set := factor(score_set, levels = unique(marker_coverage$score_set))]
p_lineage <- ggplot(lineage_tile, aes(x = score_set, y = cluster, fill = median_score_z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b91c1c", midpoint = 0, na.value = "grey90") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 11)
  ) +
  labs(title = "Cluster median lineage score z", fill = "Median z")
safe_write_png(p_lineage, file.path(output_dir, "gse164522_audit_lineage_score_heatmap.png"), 8, 6)

status_long <- melt(
  integrated[, .(
    cluster,
    lineage_boundary_flag,
    non_myeloid_score_enriched,
    doublet_enriched_15,
    sample_dominated_70,
    patient_dominated_70,
    tissue_dominated_90,
    author_subtype_dominated_80,
    small_cluster_flag
  )],
  id.vars = "cluster",
  variable.name = "flag",
  value.name = "flag_value"
)
status_long[, cluster := factor(cluster, levels = cluster_levels)]
p_flags <- ggplot(status_long, aes(x = flag, y = cluster, fill = flag_value)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("TRUE" = "#b91c1c", "FALSE" = "#e5e7eb"), na.value = "#f3f4f6") +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 11)
  ) +
  labs(title = "Cluster audit flags", fill = "Flag")
safe_write_png(p_flags, file.path(output_dir, "gse164522_audit_cluster_flag_matrix.png"), 8, 6)

writeLines(
  c(
    "# GSE164522 Macrophage/Monocyte doublet and batch audit",
    "",
    paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    paste0("Seurat object: ", rds_path),
    paste0("Cells: ", ncol(obj)),
    paste0("Genes: ", nrow(obj)),
    paste0("Cluster column: ", cluster_col),
    paste0("Clusters: ", length(cluster_levels)),
    "",
    "No clusters were removed in this audit.",
    "",
    "Package availability:",
    paste0(optional_package_status$package, ": ", optional_package_status$available),
    "",
    paste0("Doublet status: ", doublet_status),
    paste0("Nearest-neighbor batch mixing status: ", knn_status),
    paste0("Harmony status: ", harmony_status),
    "",
    "Review thresholds written as boolean flags:",
    "sample_dominated_70: max_sample_fraction >= 0.70",
    "patient_dominated_70: max_patient_fraction >= 0.70",
    "tissue_dominated_90: max_tissue_fraction >= 0.90",
    "author_subtype_dominated_80: max_celltype_sub_fraction >= 0.80",
    "doublet_enriched_15: scDblFinder_doublet_fraction >= 0.15",
    "non_myeloid_score_enriched: top_non_myeloid_fraction >= 0.50",
    "",
    "Main tables:",
    "gse164522_audit_cluster_integrated_review.tsv",
    "gse164522_audit_cluster_review_status.tsv"
  ),
  file.path(output_dir, "STATUS.md")
)

sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
