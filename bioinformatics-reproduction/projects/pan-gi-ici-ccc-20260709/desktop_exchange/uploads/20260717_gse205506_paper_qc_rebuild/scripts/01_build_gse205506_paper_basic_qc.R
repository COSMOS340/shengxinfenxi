#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
})

set.seed(340)
gse_root <- "F:/pan-gi-ici-ccc-20260709/GSE205506"
project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
previous_output <- file.path(project_dir, "desktop_exchange/uploads/20260716_stop_prjna932556_and_build_gse205506")
output_dir <- file.path(project_dir, "desktop_exchange/uploads/20260717_gse205506_paper_qc_rebuild")
analysis_dir <- file.path(gse_root, "r_analysis_paper_qc_20260717")
basic_object_path <- file.path(analysis_dir, "gse205506_article_basic_qc_primary.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

log_step <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
  flush.console()
}

required_previous <- c(
  "gse205506_matrix_structure_audit.tsv",
  "gse205506_geo_sample_metadata.tsv",
  "gse205506_response_mapping.tsv",
  "gse205506_response_join_summary.tsv",
  "gse205506_author_marker_reference_top50.tsv"
)
missing_previous <- required_previous[!file.exists(file.path(previous_output, required_previous))]
if (length(missing_previous)) stop("Missing accepted prior audit files: ", paste(missing_previous, collapse = "; "), call. = FALSE)

matrix_audit <- fread(file.path(previous_output, "gse205506_matrix_structure_audit.tsv"))
geo_metadata <- fread(file.path(previous_output, "gse205506_geo_sample_metadata.tsv"))
response_mapping <- fread(file.path(previous_output, "gse205506_response_mapping.tsv"))
response_join_summary <- fread(file.path(previous_output, "gse205506_response_join_summary.tsv"))
if (nrow(matrix_audit) != 40L || any(!matrix_audit$dimension_match) || any(!matrix_audit$metadata_match)) {
  stop("Accepted matrix audit does not pass for all 40 samples", call. = FALSE)
}
if (any(!response_join_summary$pass)) stop("Accepted response join audit contains a failure", call. = FALSE)

previous_rds <- list.files(file.path(gse_root, "r_analysis"), pattern = "\\.rds$", full.names = TRUE)
if (length(previous_rds) != 6L) stop("Expected six prior version-1 RDS files, observed ", length(previous_rds), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("R package digest is required", call. = FALSE)
rejected_inventory <- rbindlist(lapply(previous_rds, function(path) {
  data.table(
    object_name = basename(path),
    absolute_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    size_bytes = as.numeric(file.info(path)$size),
    sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE),
    status = "rejected_qc_version_1_preserved",
    reason = "Previous per-sample mitochondrial filtering did not reproduce the source-paper QC; file preserved and not overwritten"
  )
}))
write_tsv(rejected_inventory, "gse205506_rejected_qc_v1_object_inventory.tsv")

expected_geo_treatments <- c("untreated", "anti-PD-1", "Anti-PD-1+celecoxib")
if (!setequal(unique(geo_metadata$treatment), expected_geo_treatments)) {
  stop("Unexpected exact GEO treatment values", call. = FALSE)
}

feature_reference <- NULL
basic_counts <- vector("list", nrow(matrix_audit))
cell_audits <- vector("list", nrow(matrix_audit))
sample_stage_rows <- vector("list", nrow(matrix_audit))
per_sample_gene_rows <- vector("list", nrow(matrix_audit))
global_detected_cells <- NULL

log_step("Pass 1: reading all 40 author matrices and applying source-paper numerical cell bounds")
for (i in seq_len(nrow(matrix_audit))) {
  row <- matrix_audit[i]
  accession <- row$geo_accession
  features <- fread(row$features_path, header = FALSE, sep = "\t", showProgress = FALSE)
  barcodes <- fread(row$barcodes_path, header = FALSE, sep = "\t", showProgress = FALSE)
  if (ncol(features) != 3L || ncol(barcodes) != 1L) stop("Unexpected feature/barcode structure for ", accession, call. = FALSE)
  setnames(features, c("gene_id", "gene_symbol", "feature_type"))
  setnames(barcodes, "barcode_source")

  if (is.null(feature_reference)) {
    feature_reference <- copy(features)
    feature_reference[, seurat_feature_name := make.unique(gene_symbol)]
    global_detected_cells <- numeric(nrow(feature_reference))
  } else if (!identical(as.data.frame(features), as.data.frame(feature_reference[, .(gene_id, gene_symbol, feature_type)]))) {
    stop("Feature reference differs for ", accession, call. = FALSE)
  }

  con <- gzfile(row$matrix_path, open = "rt")
  counts <- tryCatch(readMM(con), finally = close(con))
  counts <- as(counts, "CsparseMatrix")
  rownames(counts) <- feature_reference$seurat_feature_name
  cell_ids <- paste0(accession, "_", barcodes$barcode_source)
  colnames(counts) <- cell_ids

  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)
  mt_index <- grepl("^MT-", feature_reference$gene_symbol)
  percent_mt <- 100 * Matrix::colSums(counts[mt_index, , drop = FALSE]) / pmax(n_count, 1)
  pass_feature_lower <- n_feature >= 500
  pass_feature_upper <- n_feature <= 5000
  pass_count_lower <- n_count >= 400
  pass_count_upper <- n_count <= 25000
  basic_pass <- pass_feature_lower & pass_feature_upper & pass_count_lower & pass_count_upper
  failure_reason <- vapply(seq_along(cell_ids), function(j) {
    reasons <- c()
    if (!pass_feature_lower[[j]]) reasons <- c(reasons, "nFeature_RNA_below_500")
    if (!pass_feature_upper[[j]]) reasons <- c(reasons, "nFeature_RNA_above_5000")
    if (!pass_count_lower[[j]]) reasons <- c(reasons, "nCount_RNA_below_400")
    if (!pass_count_upper[[j]]) reasons <- c(reasons, "nCount_RNA_above_25000")
    if (!length(reasons)) "PASS" else paste(reasons, collapse = ";")
  }, character(1))

  geo_row <- geo_metadata[geo_accession == accession]
  if (nrow(geo_row) != 1L) stop("GEO metadata row mismatch for ", accession, call. = FALSE)
  response_row <- response_mapping[subject == geo_row$subject]
  if (nrow(response_row) != 1L) stop("Response mapping row mismatch for ", accession, call. = FALSE)
  timepoint <- switch(
    geo_row$treatment,
    "untreated" = "pre-treatment",
    "anti-PD-1" = "post-treatment",
    "Anti-PD-1+celecoxib" = "post-treatment",
    stop("Unmapped GEO treatment value", call. = FALSE)
  )

  global_detected_cells <- global_detected_cells + Matrix::rowSums(counts > 0)
  per_sample_gene_rows[[i]] <- data.table(
    geo_accession = accession,
    genes_total = nrow(counts),
    genes_detected_in_at_least_3_cells_per_sample = sum(Matrix::rowSums(counts > 0) >= 3),
    interpretation = "sensitivity only; formal route applies the >=3-cell rule after merging all samples"
  )
  basic_counts[[i]] <- counts[, basic_pass, drop = FALSE]

  cell_audits[[i]] <- data.table(
    cell_id = cell_ids,
    barcode_source = barcodes$barcode_source,
    geo_accession = accession,
    geo_subject = geo_row$subject,
    geo_tissue = geo_row$tissue,
    geo_treatment = geo_row$treatment,
    derived_timepoint_from_exact_geo_treatment = timepoint,
    table_s1_response = response_row$response,
    nFeature_RNA_before_gene_filter = as.numeric(n_feature),
    nCount_RNA_before_gene_filter = as.numeric(n_count),
    percent_mt_before_gene_filter = as.numeric(percent_mt),
    pass_nFeature_RNA_ge_500 = pass_feature_lower,
    pass_nFeature_RNA_le_5000 = pass_feature_upper,
    pass_nCount_RNA_ge_400 = pass_count_lower,
    pass_nCount_RNA_le_25000 = pass_count_upper,
    paper_basic_qc_pass = basic_pass,
    paper_basic_qc_failure_reason = failure_reason,
    scDblFinder_class_sensitivity = NA_character_,
    scDblFinder_score_sensitivity = NA_real_,
    primary_inclusion_after_basic_qc = basic_pass,
    first_round_broad_compartment = NA_character_,
    mt_model_group = NA_character_,
    mt_model_center = NA_real_,
    mt_model_sigma = NA_real_,
    mt_model_p_upper = NA_real_,
    mt_model_p_bonferroni = NA_real_,
    compartment_mt_qc_pass = NA,
    formal_final_inclusion = FALSE
  )
  sample_stage_rows[[i]] <- data.table(
    geo_accession = accession,
    geo_subject = geo_row$subject,
    geo_tissue = geo_row$tissue,
    geo_treatment = geo_row$treatment,
    timepoint = timepoint,
    response = response_row$response,
    cells_input = ncol(counts),
    cells_paper_basic_qc_pass = sum(basic_pass),
    cells_paper_basic_qc_fail = sum(!basic_pass),
    basic_qc_retention_fraction = mean(basic_pass)
  )
  rm(counts)
  invisible(gc())
  log_step(sprintf("[%d/40] %s input=%d basic_pass=%d", i, accession, length(cell_ids), sum(basic_pass)))
}

merged_gene_keep <- global_detected_cells >= 3
if (sum(merged_gene_keep) < 1000L) stop("Merged >=3-cell gene rule retained fewer than 1000 genes", call. = FALSE)
gene_sensitivity <- rbindlist(per_sample_gene_rows)
gene_sensitivity <- rbind(
  gene_sensitivity,
  data.table(
    geo_accession = "MERGED_FORMAL_ROUTE",
    genes_total = length(global_detected_cells),
    genes_detected_in_at_least_3_cells_per_sample = sum(merged_gene_keep),
    interpretation = "formal route; genes detected in at least three cells across the merged raw object"
  ),
  fill = TRUE
)
write_tsv(gene_sensitivity, "gse205506_paper_qc_gene_filter_sensitivity.tsv")
feature_reference[, `:=`(
  detected_cells_merged_raw = global_detected_cells,
  formal_merged_gene_filter_pass = merged_gene_keep
)]
write_tsv(feature_reference, "gse205506_paper_qc_feature_audit.tsv")

log_step("Pass 2: per-sample scDblFinder sensitivity and Seurat object construction")
sample_objects <- vector("list", length(basic_counts))
doublet_rows <- vector("list", length(basic_counts))
for (i in seq_along(basic_counts)) {
  accession <- matrix_audit$geo_accession[[i]]
  counts <- basic_counts[[i]][merged_gene_keep, , drop = FALSE]
  if (ncol(counts) < 200L) stop("Fewer than 200 article-QC-pass cells for ", accession, call. = FALSE)
  sce <- SingleCellExperiment(list(counts = counts))
  set.seed(340 + i)
  sce <- scDblFinder(sce, verbose = FALSE)
  dbl_class <- as.character(colData(sce)[["scDblFinder.class"]])
  dbl_score <- as.numeric(colData(sce)[["scDblFinder.score"]])
  audit_index <- match(colnames(counts), cell_audits[[i]]$cell_id)
  if (anyNA(audit_index)) stop("Cell audit join failed for ", accession, call. = FALSE)
  cell_audits[[i]][audit_index, `:=`(
    scDblFinder_class_sensitivity = dbl_class,
    scDblFinder_score_sensitivity = dbl_score
  )]

  geo_row <- geo_metadata[geo_accession == accession]
  response_row <- response_mapping[subject == geo_row$subject]
  timepoint <- if (geo_row$treatment == "untreated") "pre-treatment" else "post-treatment"
  meta <- data.frame(
    row.names = colnames(counts),
    source_nCount_RNA_before_gene_filter = cell_audits[[i]][audit_index, nCount_RNA_before_gene_filter],
    source_nFeature_RNA_before_gene_filter = cell_audits[[i]][audit_index, nFeature_RNA_before_gene_filter],
    percent.mt = cell_audits[[i]][audit_index, percent_mt_before_gene_filter],
    scDblFinder.class.sensitivity = dbl_class,
    scDblFinder.score.sensitivity = dbl_score,
    geo_accession = accession,
    geo_title = geo_row$title,
    geo_source_name_ch1 = geo_row$source_name_ch1,
    geo_subject = geo_row$subject,
    geo_cell_type_source = geo_row$cell_type_source,
    geo_tissue = geo_row$tissue,
    geo_genotype = geo_row$genotype,
    geo_treatment = geo_row$treatment,
    geo_description = geo_row$description,
    derived_timepoint_from_exact_geo_treatment = timepoint,
    table_s1_response = response_row$response,
    table_s1_treatment = response_row$treatment_table_s1,
    table_s1_tumor_anatomical_location = response_row$tumor_anatomical_location,
    table_s1_mismatch_repair_defective_protein = response_row$mismatch_repair_defective_protein,
    table_s1_microsatellite_status = response_row$microsatellite_status,
    table_s1_stage_cTMN = response_row$stage_cTMN,
    table_s1_sex = response_row$sex,
    table_s1_age = response_row$age,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  object <- CreateSeuratObject(counts = counts, project = accession, assay = "RNA", min.cells = 0, min.features = 0, meta.data = meta)
  object$orig.ident <- object$geo_accession
  sample_objects[[i]] <- object
  doublet_rows[[i]] <- data.table(
    geo_accession = accession,
    cells_paper_basic_qc_pass = length(dbl_class),
    scDblFinder_singlets = sum(dbl_class == "singlet"),
    scDblFinder_doublets = sum(dbl_class == "doublet"),
    scDblFinder_doublet_fraction = mean(dbl_class == "doublet"),
    primary_object_exclusion_by_scDblFinder = 0L,
    route = "secondary sensitivity only"
  )
  rm(sce, counts, object)
  basic_counts[i] <- list(NULL)
  invisible(gc())
  log_step(sprintf("[scDblFinder %d/40] %s doublets=%d sensitivity_only", i, accession, sum(dbl_class == "doublet")))
}

cell_audit <- rbindlist(cell_audits, use.names = TRUE, fill = TRUE)
fwrite(cell_audit, file.path(output_dir, "gse205506_paper_qc_cell_audit.tsv.gz"), sep = "\t", quote = FALSE, na = "", compress = "gzip")
doublet_sensitivity <- rbindlist(doublet_rows)
write_tsv(doublet_sensitivity, "gse205506_paper_qc_doublet_sensitivity.tsv")

sample_counts <- rbindlist(sample_stage_rows)
sample_counts <- merge(sample_counts, doublet_sensitivity[, .(geo_accession, scDblFinder_singlets, scDblFinder_doublets)], by = "geo_accession", all.x = TRUE, sort = FALSE)
write_tsv(sample_counts, "gse205506_paper_qc_sample_counts.tsv")
stage_counts <- data.table(
  stage = c("raw_author_matrices", "paper_basic_qc_pass", "paper_basic_qc_fail", "scDblFinder_singlet_sensitivity", "scDblFinder_doublet_sensitivity"),
  cells = c(
    sum(sample_counts$cells_input),
    sum(sample_counts$cells_paper_basic_qc_pass),
    sum(sample_counts$cells_paper_basic_qc_fail),
    sum(sample_counts$scDblFinder_singlets),
    sum(sample_counts$scDblFinder_doublets)
  ),
  primary_route = c(TRUE, TRUE, TRUE, FALSE, FALSE),
  notes = c(
    "all barcodes from 40 author matrices",
    "500<=nFeature_RNA<=5000 and 400<=nCount_RNA<=25000; scDblFinder not used for exclusion",
    "failed one or more source-paper numerical bounds",
    "secondary sensitivity only",
    "secondary sensitivity only"
  )
)
write_tsv(stage_counts, "gse205506_paper_qc_stage_counts.tsv")

log_step("Merging 40 paper-basic-QC Seurat objects")
merged <- sample_objects[[1L]]
for (i in 2:length(sample_objects)) {
  merged <- merge(merged, sample_objects[[i]], merge.data = FALSE, merge.dr = FALSE)
  sample_objects[i] <- list(NULL)
  invisible(gc())
  log_step(sprintf("[merge %d/40]", i))
}
merged <- JoinLayers(merged, assay = "RNA")
merged@misc$paper_basic_qc <- list(
  gene_filter = "features detected in at least three cells across merged raw object",
  cell_filter = "500<=source nFeature_RNA<=5000 and 400<=source nCount_RNA<=25000",
  doublet_primary_rule = "no computational doublet exclusion; upper gene and UMI limits are the source-paper primary rule",
  scDblFinder = "per-sample secondary sensitivity only",
  seed = 340
)
saveRDS(merged, basic_object_path, compress = FALSE)

basic_summary <- data.table(
  object_path = normalizePath(basic_object_path, winslash = "/", mustWork = TRUE),
  size_bytes = as.numeric(file.info(basic_object_path)$size),
  sha256 = digest::digest(file = basic_object_path, algo = "sha256", serialize = FALSE),
  features = nrow(merged),
  cells = ncol(merged),
  samples = uniqueN(merged$geo_accession),
  subjects = uniqueN(merged$geo_subject),
  scDblFinder_doublets_retained_in_primary_object = sum(merged$scDblFinder.class.sensitivity == "doublet"),
  delivery = "local_only"
)
write_tsv(basic_summary, "gse205506_paper_qc_basic_object_summary.tsv")

parameters <- data.table(
  parameter = c(
    "seed", "merged_gene_min_detected_cells", "nFeature_RNA_lower", "nFeature_RNA_upper",
    "nCount_RNA_lower", "nCount_RNA_upper", "primary_doublet_rule", "secondary_doublet_method"
  ),
  value = c(340, 3, 500, 5000, 400, 25000, "article gene/UMI bounds; no scDblFinder exclusion", paste0("scDblFinder ", packageVersion("scDblFinder")))
)
write_tsv(parameters, "gse205506_paper_qc_parameters.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_gse205506_paper_basic_qc.txt"))
log_step("GSE205506 paper basic QC primary object complete")
