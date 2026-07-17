#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 01_gse205506_inventory_basic_qc.R <gse205506_root> <output_dir> <raw_tar_sha256>",
    call. = FALSE
  )
}

gse_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
raw_tar_sha256 <- args[[3]]

raw_tar <- file.path(gse_root, "raw", "GSE205506_RAW.tar")
extracted_dir <- file.path(gse_root, "extracted")
supplement_dir <- file.path(gse_root, "supplement")
soft_file <- file.path(supplement_dir, "GSE205506_family.soft.gz")

required_paths <- c(raw_tar, extracted_dir, supplement_dir, soft_file)
missing_paths <- required_paths[!file.exists(required_paths) & !dir.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required paths: ", paste(missing_paths, collapse = "; "), call. = FALSE)
}

write_tsv <- function(x, name) {
  fwrite(x, file.path(output_dir, name), sep = "\t", quote = FALSE, na = "")
}

read_gzip_lines <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

parse_soft_samples <- function(path) {
  lines <- read_gzip_lines(path)
  sample_start <- grep("^\\^SAMPLE = ", lines)
  if (length(sample_start) == 0L) {
    stop("No ^SAMPLE blocks found in ", path, call. = FALSE)
  }

  sample_end <- c(sample_start[-1L] - 1L, length(lines))
  rows <- lapply(seq_along(sample_start), function(i) {
    block <- lines[sample_start[[i]]:sample_end[[i]]]
    exact_value <- function(key) {
      hit <- block[startsWith(block, paste0(key, " = "))]
      if (length(hit) == 0L) return(NA_character_)
      if (length(hit) > 1L) {
        stop("Duplicate exact key ", key, " in sample block ", block[[1]], call. = FALSE)
      }
      sub(paste0("^", key, " = "), "", hit)
    }

    characteristics <- block[startsWith(block, "!Sample_characteristics_ch1 = ")]
    characteristics <- sub("^!Sample_characteristics_ch1 = ", "", characteristics)
    characteristic_map <- list()
    for (entry in characteristics) {
      split_at <- regexpr(": ", entry, fixed = TRUE)[[1]]
      if (split_at < 1L) {
        stop("Malformed exact Sample_characteristics_ch1 entry: ", entry, call. = FALSE)
      }
      field <- substr(entry, 1L, split_at - 1L)
      value <- substr(entry, split_at + 2L, nchar(entry))
      if (!is.null(characteristic_map[[field]])) {
        stop("Duplicate characteristic field ", field, " in sample block ", block[[1]], call. = FALSE)
      }
      characteristic_map[[field]] <- value
    }

    required_characteristics <- c("subject", "cell type", "tissue", "genotype", "treatment")
    absent <- required_characteristics[!required_characteristics %in% names(characteristic_map)]
    if (length(absent) > 0L) {
      stop(
        "Missing exact characteristic fields in sample block ", block[[1]], ": ",
        paste(absent, collapse = ", "), call. = FALSE
      )
    }

    data.table(
      geo_accession = exact_value("!Sample_geo_accession"),
      title = exact_value("!Sample_title"),
      source_name_ch1 = exact_value("!Sample_source_name_ch1"),
      subject = characteristic_map[["subject"]],
      cell_type_source = characteristic_map[["cell type"]],
      tissue = characteristic_map[["tissue"]],
      genotype = characteristic_map[["genotype"]],
      treatment = characteristic_map[["treatment"]],
      description = exact_value("!Sample_description")
    )
  })

  metadata <- rbindlist(rows, use.names = TRUE)
  if (anyNA(metadata$geo_accession) || anyDuplicated(metadata$geo_accession)) {
    stop("GEO sample accessions are missing or duplicated in SOFT metadata", call. = FALSE)
  }
  metadata[]
}

parse_matrix_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  header <- readLines(con, n = 1L, warn = FALSE)
  if (!identical(header, "%%MatrixMarket matrix coordinate integer general")) {
    stop("Unexpected Matrix Market header in ", path, ": ", header, call. = FALSE)
  }
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) stop("Missing matrix dimension line in ", path, call. = FALSE)
    if (!startsWith(line, "%")) break
  }
  values <- strsplit(trimws(line), "[[:space:]]+")[[1]]
  if (length(values) != 3L || anyNA(suppressWarnings(as.numeric(values)))) {
    stop("Malformed matrix dimension line in ", path, ": ", line, call. = FALSE)
  }
  setNames(as.numeric(values), c("n_features", "n_barcodes", "nnz"))
}

count_gzip_lines <- function(path) {
  length(read_gzip_lines(path))
}

archive_files <- list.files(extracted_dir, full.names = TRUE, recursive = FALSE)
archive_info <- file.info(archive_files)
archive_names <- basename(archive_files)

role_patterns <- c(
  barcodes = "_barcodes.tsv.gz$",
  features = "_features.tsv.gz$",
  matrix = "_matrix.mtx.gz$"
)

role <- rep(NA_character_, length(archive_names))
for (role_name in names(role_patterns)) {
  role[grepl(role_patterns[[role_name]], archive_names)] <- role_name
}
if (length(archive_names) != 120L || anyNA(role)) {
  stop(
    "Expected 120 observed archive members matching the three observed suffixes; found ",
    length(archive_names), " total and ", sum(is.na(role)), " unmatched", call. = FALSE
  )
}

geo_accession <- sub("^((GSM[0-9]+)).*$", "\\1", archive_names)
if (any(!grepl("^GSM[0-9]+_", archive_names))) {
  stop("At least one archive member does not begin with an observed GSM accession and underscore", call. = FALSE)
}

sample_stem <- archive_names
for (role_name in names(role_patterns)) {
  hit <- role == role_name
  sample_stem[hit] <- sub(role_patterns[[role_name]], "", archive_names[hit])
}

tar_inventory <- data.table(
  tar_member_name = archive_names,
  extracted_absolute_path = normalizePath(archive_files, winslash = "/", mustWork = TRUE),
  member_size_bytes = as.numeric(archive_info$size),
  geo_accession = geo_accession,
  sample_stem = sample_stem,
  file_role = role,
  matrix_format = fifelse(role == "matrix", "Matrix Market MTX gzip", "TSV gzip"),
  extraction_verified = TRUE
)
setorder(tar_inventory, geo_accession, file_role)
write_tsv(tar_inventory, "gse205506_tar_inventory.tsv")

geo_metadata <- parse_soft_samples(soft_file)
if (nrow(geo_metadata) != 40L) {
  stop("Expected 40 GEO sample blocks; found ", nrow(geo_metadata), call. = FALSE)
}
write_tsv(geo_metadata, "gse205506_geo_sample_metadata.tsv")

triplet_counts <- dcast(tar_inventory, geo_accession + sample_stem ~ file_role, length)
if (nrow(triplet_counts) != 40L || any(triplet_counts[, .(barcodes, features, matrix)] != 1L)) {
  stop("The observed archive members do not form exactly 40 one-to-one matrix triplets", call. = FALSE)
}

matrix_rows <- vector("list", nrow(triplet_counts))
qc_rows <- vector("list", nrow(triplet_counts))
sensitivity_rows <- list()

for (i in seq_len(nrow(triplet_counts))) {
  accession <- triplet_counts$geo_accession[[i]]
  stem <- triplet_counts$sample_stem[[i]]
  current <- tar_inventory[geo_accession == accession & sample_stem == stem]
  matrix_path <- current[file_role == "matrix", extracted_absolute_path]
  features_path <- current[file_role == "features", extracted_absolute_path]
  barcodes_path <- current[file_role == "barcodes", extracted_absolute_path]

  dimensions <- parse_matrix_header(matrix_path)
  feature_rows <- count_gzip_lines(features_path)
  barcode_rows <- count_gzip_lines(barcodes_path)
  metadata_rows <- geo_metadata[geo_accession == accession, .N]

  features <- fread(features_path, header = FALSE, sep = "\t", showProgress = FALSE)
  barcodes <- fread(barcodes_path, header = FALSE, sep = "\t", showProgress = FALSE)
  if (ncol(features) != 3L || ncol(barcodes) != 1L) {
    stop("Unexpected feature or barcode column count for ", accession, call. = FALSE)
  }

  con <- gzfile(matrix_path, open = "rt")
  counts <- tryCatch(readMM(con), finally = close(con))
  counts <- as(counts, "dgCMatrix")
  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)
  mt_index <- grepl("^MT-", features[[2]])
  percent_mt <- if (any(mt_index)) {
    100 * Matrix::colSums(counts[mt_index, , drop = FALSE]) / pmax(n_count, 1)
  } else {
    rep(NA_real_, length(n_count))
  }

  dimension_match <- identical(as.numeric(dim(counts)), as.numeric(dimensions[c("n_features", "n_barcodes")])) &&
    feature_rows == dimensions[["n_features"]] &&
    barcode_rows == dimensions[["n_barcodes"]]

  matrix_rows[[i]] <- data.table(
    geo_accession = accession,
    sample_stem = stem,
    matrix_path = matrix_path,
    features_path = features_path,
    barcodes_path = barcodes_path,
    matrix_n_features = dimensions[["n_features"]],
    feature_rows = feature_rows,
    matrix_n_barcodes = dimensions[["n_barcodes"]],
    barcode_rows = barcode_rows,
    matrix_nnz = dimensions[["nnz"]],
    metadata_rows = metadata_rows,
    dimension_match = dimension_match,
    metadata_match = metadata_rows == 1L
  )

  metric_summary <- function(x, prefix) {
    values <- c(
      min = min(x, na.rm = TRUE),
      q01 = unname(quantile(x, 0.01, na.rm = TRUE, type = 8)),
      q05 = unname(quantile(x, 0.05, na.rm = TRUE, type = 8)),
      median = median(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      q95 = unname(quantile(x, 0.95, na.rm = TRUE, type = 8)),
      q99 = unname(quantile(x, 0.99, na.rm = TRUE, type = 8)),
      max = max(x, na.rm = TRUE),
      mad = mad(x, constant = 1, na.rm = TRUE)
    )
    setNames(as.list(values), paste0(prefix, "_", names(values)))
  }

  qc_rows[[i]] <- as.data.table(c(
    list(
      geo_accession = accession,
      subject = geo_metadata[geo_accession == accession, subject],
      genotype = geo_metadata[geo_accession == accession, genotype],
      treatment = geo_metadata[geo_accession == accession, treatment],
      cells = length(n_count),
      mitochondrial_gene_count = sum(mt_index)
    ),
    metric_summary(n_count, "nCount_RNA"),
    metric_summary(n_feature, "nFeature_RNA"),
    metric_summary(percent_mt, "percent_mt")
  ))

  for (mad_multiplier in c(6, 5, 4)) {
    scenario <- c(`6` = "permissive", `5` = "primary", `4` = "stringent")[[as.character(mad_multiplier)]]
    lower_feature <- max(min(n_feature), median(n_feature) - mad_multiplier * mad(n_feature, constant = 1))
    upper_feature <- min(max(n_feature), median(n_feature) + mad_multiplier * mad(n_feature, constant = 1))
    upper_count <- min(max(n_count), median(n_count) + mad_multiplier * mad(n_count, constant = 1))
    upper_mt <- min(max(percent_mt, na.rm = TRUE), median(percent_mt, na.rm = TRUE) + mad_multiplier * mad(percent_mt, constant = 1, na.rm = TRUE))
    keep <- n_feature >= lower_feature & n_feature <= upper_feature & n_count <= upper_count & percent_mt <= upper_mt

    sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
      geo_accession = accession,
      scenario = scenario,
      mad_multiplier = mad_multiplier,
      nFeature_RNA_lower = lower_feature,
      nFeature_RNA_upper = upper_feature,
      nCount_RNA_upper = upper_count,
      percent_mt_upper = upper_mt,
      cells_total = length(keep),
      cells_retained = sum(keep),
      retention_fraction = mean(keep)
    )
  }

  rm(counts, n_count, n_feature, percent_mt)
  invisible(gc())
  message(sprintf("[%d/40] audited %s", i, accession))
}

matrix_audit <- rbindlist(matrix_rows, use.names = TRUE)
qc_summary <- rbindlist(qc_rows, use.names = TRUE, fill = TRUE)
qc_sensitivity <- rbindlist(sensitivity_rows, use.names = TRUE)
setorder(matrix_audit, geo_accession)
setorder(qc_summary, geo_accession)
setorder(qc_sensitivity, geo_accession, -mad_multiplier)

if (any(!matrix_audit$dimension_match) || any(!matrix_audit$metadata_match)) {
  stop("At least one matrix failed dimension or GEO metadata matching", call. = FALSE)
}

write_tsv(matrix_audit, "gse205506_matrix_structure_audit.tsv")
write_tsv(qc_summary, "gse205506_qc_summary.tsv")
write_tsv(qc_sensitivity, "gse205506_qc_sensitivity_audit.tsv")

inventory_paths <- c(
  raw_tar,
  list.files(supplement_dir, full.names = TRUE),
  tar_inventory$extracted_absolute_path
)
inventory_info <- file.info(inventory_paths)
input_inventory <- data.table(
  absolute_path = normalizePath(inventory_paths, winslash = "/", mustWork = TRUE),
  filename = basename(inventory_paths),
  size_bytes = as.numeric(inventory_info$size),
  source = fifelse(
    normalizePath(inventory_paths, winslash = "/", mustWork = TRUE) == normalizePath(raw_tar, winslash = "/", mustWork = TRUE),
    "GEO GSE205506 supplementary file",
    fifelse(startsWith(normalizePath(inventory_paths, winslash = "/", mustWork = TRUE), normalizePath(extracted_dir, winslash = "/", mustWork = TRUE)),
      "Extracted from GSE205506_RAW.tar",
      "GEO family metadata or publication-related local file"
    )
  ),
  sha256 = fifelse(basename(inventory_paths) == "GSE205506_RAW.tar", raw_tar_sha256, NA_character_)
)
setorder(input_inventory, absolute_path)
write_tsv(input_inventory, "gse205506_input_inventory.tsv")

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_inventory_basic_qc.txt"))
message("Inventory and basic QC completed")
