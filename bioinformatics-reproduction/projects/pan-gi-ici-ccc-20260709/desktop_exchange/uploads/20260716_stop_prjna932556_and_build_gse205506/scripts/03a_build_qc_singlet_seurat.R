#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(scuttle)
  library(mclust)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 03a_build_qc_singlet_seurat.R <gse205506_root> <output_dir> <seed>",
    call. = FALSE
  )
}

gse_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
seed <- as.integer(args[[3]])
analysis_dir <- file.path(gse_root, "r_analysis")
per_sample_dir <- file.path(analysis_dir, "per_sample_filtered")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(per_sample_dir, recursive = TRUE, showWarnings = FALSE)

matrix_audit_path <- file.path(output_dir, "gse205506_matrix_structure_audit.tsv")
geo_metadata_path <- file.path(output_dir, "gse205506_geo_sample_metadata.tsv")
response_mapping_path <- file.path(output_dir, "gse205506_response_mapping.tsv")
response_join_summary_path <- file.path(output_dir, "gse205506_response_join_summary.tsv")

required_paths <- c(matrix_audit_path, geo_metadata_path, response_mapping_path, response_join_summary_path)
if (any(!file.exists(required_paths))) {
  stop("Missing required audit files: ", paste(required_paths[!file.exists(required_paths)], collapse = "; "), call. = FALSE)
}

write_tsv <- function(x, filename) {
  fwrite(x, file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

matrix_audit <- fread(matrix_audit_path, sep = "\t", header = TRUE, na.strings = NULL)
geo_metadata <- fread(geo_metadata_path, sep = "\t", header = TRUE, na.strings = NULL)
response_mapping <- fread(response_mapping_path, sep = "\t", header = TRUE, na.strings = NULL)
response_join_summary <- fread(response_join_summary_path, sep = "\t", header = TRUE, na.strings = NULL)

if (nrow(matrix_audit) != 40L || any(!matrix_audit$dimension_match) || any(!matrix_audit$metadata_match)) {
  stop("Matrix structure audit did not pass for all 40 samples", call. = FALSE)
}
if (any(!response_join_summary$pass)) {
  stop("Response mapping audit did not pass", call. = FALSE)
}

expected_geo_treatments <- c("untreated", "anti-PD-1", "Anti-PD-1+celecoxib")
if (!setequal(unique(geo_metadata$treatment), expected_geo_treatments)) {
  stop("Unexpected exact GEO treatment values: ", paste(unique(geo_metadata$treatment), collapse = ", "), call. = FALSE)
}

mt_mixture_audit <- function(percent_mt) {
  bounded <- pmin(pmax(percent_mt, 0.01), 99.99)
  transformed <- qlogis(bounded / 100)
  fit <- tryCatch(
    Mclust(transformed, G = 1:3, modelNames = "V", verbose = FALSE),
    error = function(e) NULL
  )

  fallback_flag <- isOutlier(percent_mt, type = "higher", nmads = 3)
  fallback_threshold <- as.numeric(attr(fallback_flag, "thresholds")[["higher"]])
  if (!is.finite(fallback_threshold)) fallback_threshold <- max(percent_mt, na.rm = TRUE)

  if (is.null(fit)) {
    return(list(
      method = "upper_3_MAD_mclust_failed",
      fit = NULL,
      high_component = NA_integer_,
      posterior_high = as.numeric(percent_mt > fallback_threshold),
      primary_threshold = fallback_threshold,
      component_means_percent_mt = NA_real_,
      component_gap_percent_mt = NA_real_
    ))
  }

  means <- as.numeric(fit$parameters$mean)

  if (fit$G < 2L) {
    return(list(
      method = "upper_3_MAD_no_separated_mixture",
      fit = fit,
      high_component = NA_integer_,
      posterior_high = as.numeric(percent_mt > fallback_threshold),
      primary_threshold = fallback_threshold,
      component_means_percent_mt = plogis(means) * 100,
      component_gap_percent_mt = NA_real_
    ))
  }

  order_means <- order(means)
  high_component <- order_means[[length(order_means)]]
  next_component <- order_means[[length(order_means) - 1L]]
  component_means_percent <- plogis(means) * 100
  component_gap <- component_means_percent[[high_component]] - component_means_percent[[next_component]]

  if (!is.finite(component_gap) || component_gap < 10) {
    return(list(
      method = "upper_3_MAD_mixture_gap_below_10_percent",
      fit = fit,
      high_component = NA_integer_,
      posterior_high = as.numeric(percent_mt > fallback_threshold),
      primary_threshold = fallback_threshold,
      component_means_percent_mt = component_means_percent,
      component_gap_percent_mt = component_gap
    ))
  }

  posterior_high <- fit$z[, high_component]
  primary_high <- posterior_high >= 0.5
  primary_threshold <- if (any(primary_high)) min(percent_mt[primary_high], na.rm = TRUE) else max(percent_mt, na.rm = TRUE)
  list(
    method = paste0("mclust_G", fit$G, "_posterior_high_component"),
    fit = fit,
    high_component = high_component,
    posterior_high = posterior_high,
    primary_threshold = primary_threshold,
    component_means_percent_mt = component_means_percent,
    component_gap_percent_mt = component_gap
  )
}

threshold_value <- function(flag, side) {
  thresholds <- attr(flag, "thresholds")
  if (is.null(thresholds) || length(thresholds) == 0L) return(NA_real_)
  if (is.matrix(thresholds)) {
    if (!(side %in% rownames(thresholds)) || ncol(thresholds) != 1L) return(NA_real_)
    value <- as.numeric(thresholds[side, 1L])
  } else {
    if (!(side %in% names(thresholds))) return(NA_real_)
    value <- as.numeric(thresholds[[side]])
  }
  if (!is.finite(value)) NA_real_ else value
}

sample_paths <- character(nrow(matrix_audit))
qc_rows <- vector("list", nrow(matrix_audit))
doublet_rows <- vector("list", nrow(matrix_audit))
ambient_rows <- vector("list", nrow(matrix_audit))
sensitivity_rows <- list()
feature_reference <- NULL

for (i in seq_len(nrow(matrix_audit))) {
  row <- matrix_audit[i]
  accession <- row$geo_accession
  features <- fread(row$features_path, header = FALSE, sep = "\t", showProgress = FALSE)
  barcodes <- fread(row$barcodes_path, header = FALSE, sep = "\t", showProgress = FALSE)
  if (ncol(features) != 3L || ncol(barcodes) != 1L) {
    stop("Unexpected feature or barcode structure for ", accession, call. = FALSE)
  }
  setnames(features, c("gene_id", "gene_symbol", "feature_type"))
  setnames(barcodes, "barcode_source")

  if (is.null(feature_reference)) {
    feature_reference <- copy(features)
    feature_reference[, seurat_feature_name := make.unique(gene_symbol)]
    write_tsv(feature_reference, "gse205506_feature_reference.tsv")
  } else if (!identical(
    as.data.frame(features),
    as.data.frame(feature_reference[, .(gene_id, gene_symbol, feature_type)])
  )) {
    stop("Feature table differs from the first verified sample for ", accession, call. = FALSE)
  }

  con <- gzfile(row$matrix_path, open = "rt")
  counts <- tryCatch(readMM(con), finally = close(con))
  counts <- as(counts, "CsparseMatrix")
  rownames(counts) <- feature_reference$seurat_feature_name
  colnames(counts) <- paste0(accession, "_", barcodes$barcode_source)

  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)
  mt_index <- grepl("^MT-", feature_reference$gene_symbol)
  percent_mt <- 100 * Matrix::colSums(counts[mt_index, , drop = FALSE]) / pmax(n_count, 1)

  mt_audit <- mt_mixture_audit(percent_mt)
  if (startsWith(mt_audit$method, "mclust_")) {
    mt_primary_flag <- mt_audit$posterior_high >= 0.5
  } else {
    mt_primary_flag <- percent_mt > mt_audit$primary_threshold
  }
  provisional_low_mt <- !mt_primary_flag & n_count > 0 & n_feature > 0
  if (sum(provisional_low_mt) < 200L) {
    stop("Fewer than 200 provisional low-mitochondrial cells for ", accession, call. = FALSE)
  }

  count_low_reference <- isOutlier(n_count[provisional_low_mt], type = "lower", nmads = 3, log = TRUE)
  feature_low_reference <- isOutlier(n_feature[provisional_low_mt], type = "lower", nmads = 3, log = TRUE)
  count_lower <- threshold_value(count_low_reference, "lower")
  feature_lower <- threshold_value(feature_low_reference, "lower")
  if (!is.finite(count_lower)) count_lower <- min(n_count[provisional_low_mt])
  if (!is.finite(feature_lower)) feature_lower <- min(n_feature[provisional_low_mt])

  count_low <- n_count < count_lower
  feature_low <- n_feature < feature_lower
  pre_qc_keep <- n_count > 0 & n_feature > 0 & !mt_primary_flag & !count_low & !feature_low
  if (sum(pre_qc_keep) < 200L) {
    stop("Fewer than 200 cells passed primary QC for ", accession, call. = FALSE)
  }

  sce <- SingleCellExperiment(list(counts = counts[, pre_qc_keep, drop = FALSE]))
  set.seed(seed + i)
  sce <- scDblFinder(sce, verbose = FALSE)
  dbl_class <- as.character(colData(sce)[["scDblFinder.class"]])
  dbl_score <- as.numeric(colData(sce)[["scDblFinder.score"]])
  singlet_keep <- dbl_class == "singlet"
  kept_cell_names <- colnames(sce)[singlet_keep]

  geo_row <- geo_metadata[geo_accession == accession]
  if (nrow(geo_row) != 1L) stop("GEO metadata row mismatch for ", accession, call. = FALSE)
  response_row <- response_mapping[subject == geo_row$subject]
  if (nrow(response_row) != 1L) stop("Response mapping row mismatch for ", accession, call. = FALSE)

  timepoint <- switch(
    geo_row$treatment,
    "untreated" = "pre-treatment",
    "anti-PD-1" = "post-treatment",
    "Anti-PD-1+celecoxib" = "post-treatment",
    stop("Unmapped exact GEO treatment value for ", accession, call. = FALSE)
  )

  object_counts <- counts[, kept_cell_names, drop = FALSE]
  object <- CreateSeuratObject(
    counts = object_counts,
    project = accession,
    assay = "RNA",
    min.cells = 0,
    min.features = 0
  )
  kept_index <- match(kept_cell_names, colnames(counts))
  score_index <- match(kept_cell_names, colnames(sce))
  object[["nCount_RNA_source_qc"]] <- n_count[kept_index]
  object[["nFeature_RNA_source_qc"]] <- n_feature[kept_index]
  object[["percent.mt"]] <- percent_mt[kept_index]
  object[["scDblFinder.score"]] <- dbl_score[score_index]
  object[["scDblFinder.class"]] <- dbl_class[score_index]
  object[["geo_accession"]] <- accession
  object[["geo_title"]] <- geo_row$title
  object[["geo_source_name_ch1"]] <- geo_row$source_name_ch1
  object[["geo_subject"]] <- geo_row$subject
  object[["geo_cell_type_source"]] <- geo_row$cell_type_source
  object[["geo_tissue"]] <- geo_row$tissue
  object[["geo_genotype"]] <- geo_row$genotype
  object[["geo_treatment"]] <- geo_row$treatment
  object[["geo_description"]] <- geo_row$description
  object[["derived_timepoint_from_exact_geo_treatment"]] <- timepoint
  object[["table_s1_response"]] <- response_row$response
  object[["table_s1_treatment"]] <- response_row$treatment_table_s1
  object[["table_s1_tumor_anatomical_location"]] <- response_row$tumor_anatomical_location
  object[["table_s1_mismatch_repair_defective_protein"]] <- response_row$mismatch_repair_defective_protein
  object[["table_s1_microsatellite_status"]] <- response_row$microsatellite_status
  object[["table_s1_stage_cTMN"]] <- response_row$stage_cTMN
  object[["table_s1_sex"]] <- response_row$sex
  object[["table_s1_age"]] <- response_row$age

  sample_rds <- file.path(per_sample_dir, paste0(accession, ".rds"))
  saveRDS(object, sample_rds, compress = FALSE)
  sample_paths[[i]] <- sample_rds

  qc_rows[[i]] <- data.table(
    geo_accession = accession,
    subject = geo_row$subject,
    genotype = geo_row$genotype,
    treatment = geo_row$treatment,
    timepoint = timepoint,
    response = response_row$response,
    cells_input = ncol(counts),
    cells_zero_count = sum(n_count == 0),
    cells_high_mitochondrial = sum(mt_primary_flag),
    cells_low_count = sum(count_low),
    cells_low_feature = sum(feature_low),
    cells_pass_primary_qc = sum(pre_qc_keep),
    cells_scDblFinder_doublet = sum(dbl_class == "doublet"),
    cells_final_singlet = ncol(object),
    final_retention_fraction = ncol(object) / ncol(counts),
    count_lower_threshold = count_lower,
    feature_lower_threshold = feature_lower,
    percent_mt_primary_threshold = mt_audit$primary_threshold,
    percent_mt_method = mt_audit$method,
    percent_mt_component_means = paste(round(mt_audit$component_means_percent_mt, 4), collapse = ";"),
    percent_mt_component_gap = mt_audit$component_gap_percent_mt,
    final_nCount_RNA_median = median(object$nCount_RNA_source_qc),
    final_nFeature_RNA_median = median(object$nFeature_RNA_source_qc),
    final_percent_mt_median = median(object$percent.mt),
    per_sample_rds = normalizePath(sample_rds, winslash = "/", mustWork = TRUE),
    per_sample_rds_bytes = file.info(sample_rds)$size
  )

  doublet_rows[[i]] <- data.table(
    geo_accession = accession,
    cells_entered_scDblFinder = length(dbl_class),
    singlets = sum(dbl_class == "singlet"),
    doublets = sum(dbl_class == "doublet"),
    doublet_fraction = mean(dbl_class == "doublet"),
    scDblFinder_threshold = as.numeric(metadata(sce)[["scDblFinder.threshold"]]),
    score_median = median(dbl_score),
    score_q95 = as.numeric(quantile(dbl_score, 0.95, type = 8)),
    status = "completed_per_sample"
  )

  ambient_rows[[i]] <- data.table(
    geo_accession = accession,
    filtered_matrix_available = TRUE,
    unfiltered_droplet_matrix_available = FALSE,
    soup_profile_estimable = FALSE,
    correction_applied = FALSE,
    reason = "GEO archive contains filtered matrix triplets only; no unfiltered droplet matrix was provided for contamination estimation"
  )

  for (scenario in c("permissive", "primary", "stringent")) {
    posterior_cutoff <- c(permissive = 0.9, primary = 0.5, stringent = 0.1)[[scenario]]
    nmads <- c(permissive = 3.5, primary = 3.0, stringent = 2.5)[[scenario]]
    if (startsWith(mt_audit$method, "mclust_")) {
      mt_flag_scenario <- mt_audit$posterior_high >= posterior_cutoff
      mt_threshold_scenario <- if (any(mt_flag_scenario)) min(percent_mt[mt_flag_scenario]) else max(percent_mt)
    } else {
      mt_ref_flag <- isOutlier(percent_mt, type = "higher", nmads = nmads)
      mt_threshold_scenario <- threshold_value(mt_ref_flag, "higher")
      if (!is.finite(mt_threshold_scenario)) mt_threshold_scenario <- max(percent_mt)
      mt_flag_scenario <- percent_mt > mt_threshold_scenario
    }
    low_mt_scenario <- !mt_flag_scenario & n_count > 0 & n_feature > 0
    if (!any(low_mt_scenario)) {
      sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
        geo_accession = accession,
        scenario = scenario,
        mt_high_component_posterior_cutoff = if (startsWith(mt_audit$method, "mclust_")) posterior_cutoff else NA_real_,
        nmads = nmads,
        count_lower_threshold = NA_real_,
        feature_lower_threshold = NA_real_,
        percent_mt_threshold = mt_threshold_scenario,
        cells_input = ncol(counts),
        cells_pass_qc_before_doublet = 0L,
        retention_fraction_before_doublet = 0
      )
      next
    }
    count_ref <- isOutlier(n_count[low_mt_scenario], type = "lower", nmads = nmads, log = TRUE)
    feature_ref <- isOutlier(n_feature[low_mt_scenario], type = "lower", nmads = nmads, log = TRUE)
    count_threshold <- threshold_value(count_ref, "lower")
    feature_threshold <- threshold_value(feature_ref, "lower")
    if (!is.finite(count_threshold)) count_threshold <- min(n_count[low_mt_scenario])
    if (!is.finite(feature_threshold)) feature_threshold <- min(n_feature[low_mt_scenario])
    scenario_keep <- low_mt_scenario & n_count >= count_threshold & n_feature >= feature_threshold
    sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
      geo_accession = accession,
      scenario = scenario,
      mt_high_component_posterior_cutoff = if (startsWith(mt_audit$method, "mclust_")) posterior_cutoff else NA_real_,
      nmads = nmads,
      count_lower_threshold = count_threshold,
      feature_lower_threshold = feature_threshold,
      percent_mt_threshold = mt_threshold_scenario,
      cells_input = ncol(counts),
      cells_pass_qc_before_doublet = sum(scenario_keep),
      retention_fraction_before_doublet = mean(scenario_keep)
    )
  }

  rm(counts, object_counts, object, sce)
  invisible(gc())
  message(sprintf("[%d/40] %s final singlets: %d", i, accession, qc_rows[[i]]$cells_final_singlet))
}

qc_summary <- rbindlist(qc_rows, use.names = TRUE, fill = TRUE)
doublet_audit <- rbindlist(doublet_rows, use.names = TRUE, fill = TRUE)
ambient_audit <- rbindlist(ambient_rows, use.names = TRUE, fill = TRUE)
qc_sensitivity <- rbindlist(sensitivity_rows, use.names = TRUE, fill = TRUE)
setorder(qc_summary, geo_accession)
setorder(doublet_audit, geo_accession)
setorder(ambient_audit, geo_accession)
setorder(qc_sensitivity, geo_accession, -nmads)

write_tsv(qc_summary, "gse205506_qc_filter_summary.tsv")
write_tsv(doublet_audit, "gse205506_doublet_audit.tsv")
write_tsv(ambient_audit, "gse205506_environment_rna_audit.tsv")
write_tsv(qc_sensitivity, "gse205506_qc_filter_sensitivity_audit.tsv")

merged <- readRDS(sample_paths[[1]])
for (i in seq.int(2L, length(sample_paths))) {
  next_object <- readRDS(sample_paths[[i]])
  merged <- merge(
    x = merged,
    y = next_object,
    merge.data = FALSE,
    merge.dr = FALSE
  )
  rm(next_object)
  invisible(gc())
  message(sprintf("[merge %d/40] %s", i, basename(sample_paths[[i]])))
}
merged$orig.ident <- merged$geo_accession
merged@misc$qc_method <- list(
  percent_mt = "per-sample mclust high-component posterior >= 0.5 when separated by at least 10 percentage points; otherwise upper 3 MAD",
  lower_counts = "per-sample lower 3 MAD on log-transformed counts among provisional low-mitochondrial cells",
  lower_features = "per-sample lower 3 MAD on log-transformed detected features among provisional low-mitochondrial cells",
  doublets = "per-sample scDblFinder 1.20.2 classification"
)
merged@misc$response_mapping_source <- normalizePath(response_mapping_path, winslash = "/", mustWork = TRUE)

merged_rds <- file.path(analysis_dir, "gse205506_seurat_qc_singlets.rds")
saveRDS(merged, merged_rds, compress = FALSE)
merged_summary <- data.table(
  object_path = normalizePath(merged_rds, winslash = "/", mustWork = TRUE),
  object_size_bytes = file.info(merged_rds)$size,
  features = nrow(merged),
  cells = ncol(merged),
  samples = uniqueN(merged$geo_accession),
  subjects = uniqueN(merged$geo_subject),
  pCR_cells = sum(merged$table_s1_response == "pCR"),
  non_pCR_cells = sum(merged$table_s1_response == "non-pCR"),
  pre_treatment_cells = sum(merged$derived_timepoint_from_exact_geo_treatment == "pre-treatment"),
  post_treatment_cells = sum(merged$derived_timepoint_from_exact_geo_treatment == "post-treatment")
)
write_tsv(merged_summary, "gse205506_qc_singlet_object_summary.tsv")

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_qc_singlet_object.txt"))
message("QC singlet Seurat object completed: ", merged_rds)
