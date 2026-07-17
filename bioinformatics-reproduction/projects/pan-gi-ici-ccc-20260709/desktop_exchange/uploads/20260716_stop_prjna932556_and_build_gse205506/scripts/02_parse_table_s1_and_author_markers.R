#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: 02_parse_table_s1_and_author_markers.R <gse205506_root> <output_dir> <mmc2_sha256> <mmc3_sha256>",
    call. = FALSE
  )
}

gse_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
mmc2_sha256 <- args[[3]]
mmc3_sha256 <- args[[4]]
mmc2_path <- file.path(gse_root, "supplement", "mmc2.xlsx")
mmc3_path <- file.path(gse_root, "supplement", "mmc3.xlsx")
geo_metadata_path <- file.path(output_dir, "gse205506_geo_sample_metadata.tsv")
input_inventory_path <- file.path(output_dir, "gse205506_input_inventory.tsv")

required_paths <- c(mmc2_path, mmc3_path, geo_metadata_path, input_inventory_path)
if (any(!file.exists(required_paths))) {
  stop("Missing required files: ", paste(required_paths[!file.exists(required_paths)], collapse = "; "), call. = FALSE)
}

write_tsv <- function(x, filename) {
  fwrite(x, file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

input_inventory <- fread(input_inventory_path, sep = "\t", header = TRUE, na.strings = NULL)
input_inventory <- input_inventory[!filename %chin% c("mmc2.xlsx", "mmc3.xlsx")]
supplement_inventory <- data.table(
  absolute_path = c(
    normalizePath(mmc2_path, winslash = "/", mustWork = TRUE),
    normalizePath(mmc3_path, winslash = "/", mustWork = TRUE)
  ),
  filename = c("mmc2.xlsx", "mmc3.xlsx"),
  size_bytes = as.numeric(file.info(c(mmc2_path, mmc3_path))$size),
  source = c(
    "User-provided original Table S1 workbook for DOI 10.1016/j.ccell.2023.04.011 after publisher redirect loop",
    "User-provided author compartment-specific cluster marker workbook"
  ),
  sha256 = c(mmc2_sha256, mmc3_sha256)
)
input_inventory <- rbindlist(list(input_inventory, supplement_inventory), use.names = TRUE, fill = TRUE)
setorder(input_inventory, filename, absolute_path)
write_tsv(input_inventory, "gse205506_input_inventory.tsv")

clean_cell <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

raw_s1 <- read_excel(mmc2_path, sheet = "Sheet1", col_names = FALSE, .name_repair = "minimal")
if (!identical(dim(raw_s1), c(24L, 12L))) {
  stop("Unexpected mmc2.xlsx Sheet1 dimensions: ", paste(dim(raw_s1), collapse = " x "), call. = FALSE)
}

expected_title <- "Supplementary Table 1. Clinical characteristics of 19 d-MMR/MSI-H CRC patients, related to Figure 1."
if (!identical(clean_cell(raw_s1[[1]][[1]]), expected_title)) {
  stop("Unexpected Table S1 title", call. = FALSE)
}

expected_main_headers <- c(
  "Treatment", "Patient ID", "Tumor anatomical location", "Mismatch repair defective protein",
  "Microsatellite status", "Stage(cTMN)", "Sex", "Age", "Tumor", "", "Normal", "Treatment response#"
)
expected_sub_headers <- c("", "", "", "", "", "", "", "", "pre-treatment", "Post-treatment", "Post-treatment", "")
observed_main_headers <- unname(vapply(raw_s1[2, ], clean_cell, character(1)))
observed_sub_headers <- unname(vapply(raw_s1[3, ], clean_cell, character(1)))
if (!identical(observed_main_headers, expected_main_headers)) {
  stop("Table S1 main headers differ from the verified workbook structure", call. = FALSE)
}
if (!identical(observed_sub_headers, expected_sub_headers)) {
  stop("Table S1 subheaders differ from the verified workbook structure", call. = FALSE)
}

s1 <- as.data.table(raw_s1[4:22, ])
setnames(s1, c(
  "treatment_table_s1", "subject", "tumor_anatomical_location", "mismatch_repair_defective_protein",
  "microsatellite_status", "stage_cTMN", "sex", "age", "tumor_pre_treatment",
  "tumor_post_treatment", "normal_post_treatment", "response"
))
s1[, names(s1) := lapply(.SD, clean_cell)]

last_treatment <- ""
for (i in seq_len(nrow(s1))) {
  if (nzchar(s1$treatment_table_s1[[i]])) {
    last_treatment <- s1$treatment_table_s1[[i]]
  } else {
    s1$treatment_table_s1[[i]] <- last_treatment
  }
}

if (nrow(s1) != 19L || any(!nzchar(s1$subject)) || anyDuplicated(s1$subject)) {
  stop("Table S1 patient identifiers are incomplete or duplicated", call. = FALSE)
}
if (!setequal(unique(s1$response), c("pCR", "non-pCR"))) {
  stop("Unexpected Table S1 response values: ", paste(unique(s1$response), collapse = ", "), call. = FALSE)
}
response_counts <- s1[, .N, by = response]
if (response_counts[response == "pCR", N] != 15L || response_counts[response == "non-pCR", N] != 4L) {
  stop("Table S1 response counts are not 15 pCR and 4 non-pCR", call. = FALSE)
}

s1[, age := as.integer(age)]
s1[, `:=`(
  source_file = normalizePath(mmc2_path, winslash = "/", mustWork = TRUE),
  source_sheet = "Sheet1",
  source_sha256 = mmc2_sha256,
  source_url = "https://www.cell.com/cms/10.1016/j.ccell.2023.04.011/attachment/c9110253-41ec-458a-bc9d-c6be5e3f9b21/mmc2.xlsx",
  source_title = expected_title
)]
setcolorder(s1, c(
  "subject", "response", "treatment_table_s1", "tumor_anatomical_location",
  "mismatch_repair_defective_protein", "microsatellite_status", "stage_cTMN", "sex", "age",
  "tumor_pre_treatment", "tumor_post_treatment", "normal_post_treatment",
  "source_file", "source_sheet", "source_sha256", "source_url", "source_title"
))
s1[, patient_order__ := as.integer(sub("^P", "", subject))]
setorder(s1, patient_order__)
s1[, patient_order__ := NULL]
write_tsv(s1, "gse205506_response_mapping.tsv")

geo_metadata <- fread(geo_metadata_path, sep = "\t", header = TRUE, na.strings = NULL)
if (!all(c("geo_accession", "subject") %in% names(geo_metadata))) {
  stop("GEO metadata lacks exact geo_accession or subject columns", call. = FALSE)
}

mapping_subject_counts <- s1[, .(response_mapping_rows = .N), by = subject]
join_audit <- merge(
  geo_metadata[, .(geo_accession, subject_geo = subject, title, tissue, genotype, treatment)],
  s1[, .(subject_geo = subject, response, treatment_table_s1)],
  by = "subject_geo",
  all.x = TRUE,
  allow.cartesian = FALSE
)
join_audit <- merge(
  join_audit,
  mapping_subject_counts,
  by.x = "subject_geo",
  by.y = "subject",
  all.x = TRUE,
  allow.cartesian = FALSE
)
join_audit[, `:=`(
  match_status = fifelse(is.na(response), "unmatched", "matched"),
  response_blank = is.na(response) | !nzchar(response),
  response_mapping_duplicate = is.na(response_mapping_rows) | response_mapping_rows != 1L,
  many_to_many = FALSE
)]
setorder(join_audit, geo_accession)
write_tsv(join_audit, "gse205506_response_join_audit.tsv")

unmatched_mapping_subjects <- setdiff(s1$subject, unique(geo_metadata$subject))
join_summary <- data.table(
  metric = c(
    "table_s1_patients", "table_s1_pCR", "table_s1_non_pCR", "geo_samples", "geo_subjects",
    "geo_samples_matched", "geo_samples_unmatched", "mapping_subjects_unmatched_to_geo",
    "mapping_subject_duplicates", "many_to_many_connections", "response_field_blanks"
  ),
  value = c(
    nrow(s1), s1[response == "pCR", .N], s1[response == "non-pCR", .N], nrow(geo_metadata),
    uniqueN(geo_metadata$subject), join_audit[match_status == "matched", .N],
    join_audit[match_status == "unmatched", .N], length(unmatched_mapping_subjects),
    s1[, sum(duplicated(subject))], join_audit[many_to_many == TRUE, .N],
    join_audit[response_blank == TRUE, .N]
  ),
  pass = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
)
join_summary[metric == "geo_samples_unmatched", pass := value == 0]
join_summary[metric == "mapping_subjects_unmatched_to_geo", pass := value == 0]
join_summary[metric == "mapping_subject_duplicates", pass := value == 0]
join_summary[metric == "many_to_many_connections", pass := value == 0]
join_summary[metric == "response_field_blanks", pass := value == 0]
write_tsv(join_summary, "gse205506_response_join_summary.tsv")
if (any(!join_summary$pass)) {
  stop("Response mapping failed one or more exact join checks", call. = FALSE)
}

expected_marker_sheets <- c(
  "T_I_NK cell markers", "B cell markers", "Myeloid cell markers",
  "Endothelial cell markers", "Fibroblast markers"
)
observed_marker_sheets <- excel_sheets(mmc3_path)
if (!identical(observed_marker_sheets, expected_marker_sheets)) {
  stop("mmc3.xlsx sheet names differ from the verified workbook", call. = FALSE)
}

expected_marker_headers <- c(
  "Gene_name", "pvalue", "avg_log2FC", "pct exp in cluster",
  "pct exp in other cluster", "p_val_adj", "cluster"
)
marker_tables <- lapply(observed_marker_sheets, function(sheet_name) {
  title_row <- read_excel(mmc3_path, sheet = sheet_name, col_names = FALSE, n_max = 1L, .name_repair = "minimal")
  title <- clean_cell(title_row[[1]][[1]])
  markers <- as.data.table(read_excel(mmc3_path, sheet = sheet_name, skip = 1L, .name_repair = "minimal"))
  if (!identical(names(markers), expected_marker_headers)) {
    stop("Unexpected marker headers in sheet ", sheet_name, call. = FALSE)
  }
  markers[, gene_source := clean_cell(Gene_name)]
  markers[, gene := sub('^"(.*)"$', "\\1", gene_source)]
  if (any(grepl('^"|"$', markers$gene))) {
    stop("Unresolved outer quote in marker gene field on sheet ", sheet_name, call. = FALSE)
  }
  markers[, `:=`(
    compartment_sheet = sheet_name,
    source_title = title,
    source_file = normalizePath(mmc3_path, winslash = "/", mustWork = TRUE),
    source_sha256 = mmc3_sha256
  )]
  markers
})
author_markers <- rbindlist(marker_tables, use.names = TRUE, fill = TRUE)
author_markers[, avg_log2FC := as.numeric(avg_log2FC)]
author_markers[, p_val_adj := as.numeric(p_val_adj)]
if (any(!nzchar(author_markers$gene)) || any(!nzchar(author_markers$cluster))) {
  stop("Author marker table has blank gene or cluster fields", call. = FALSE)
}

setorder(author_markers, compartment_sheet, cluster, -avg_log2FC, p_val_adj, gene)
author_top10 <- author_markers[, {
  top_rows <- head(.SD, 10L)
  top_rows[, rank := seq_len(.N)]
  top_rows
}, by = .(compartment_sheet, cluster)]
author_top10 <- author_top10[, .(
  compartment_sheet, cluster, rank, gene, gene_source,
  avg_log2FC, pvalue, p_val_adj, pct_exp_in_cluster = `pct exp in cluster`,
  pct_exp_in_other_cluster = `pct exp in other cluster`, source_title, source_file, source_sha256
)]
write_tsv(author_top10, "gse205506_author_marker_reference_top10.tsv")

author_top50 <- author_markers[, {
  top_rows <- head(.SD, 50L)
  top_rows[, rank := seq_len(.N)]
  top_rows
}, by = .(compartment_sheet, cluster)]
author_top50 <- author_top50[, .(
  compartment_sheet, cluster, rank, gene, gene_source,
  avg_log2FC, pvalue, p_val_adj, pct_exp_in_cluster = `pct exp in cluster`,
  pct_exp_in_other_cluster = `pct exp in other cluster`, source_title, source_file, source_sha256
)]
write_tsv(author_top50, "gse205506_author_marker_reference_top50.tsv")

author_cluster_summary <- author_markers[, .(
  marker_rows = .N,
  positive_log2fc_rows = sum(avg_log2FC > 0, na.rm = TRUE),
  top10_genes = paste(head(gene, 10L), collapse = ";")
), by = .(compartment_sheet, cluster)]
setorder(author_cluster_summary, compartment_sheet, cluster)
write_tsv(author_cluster_summary, "gse205506_author_marker_cluster_summary.tsv")

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_response_and_author_markers.txt"))
message("Table S1 mapping and author marker reference parsing completed")
