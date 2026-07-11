suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
exchange_dir <- file.path(project_dir, "desktop_exchange")
prior_dir <- file.path(exchange_dir, "uploads", "20260709_response_object_construction")
out_dir <- file.path(exchange_dir, "uploads", "20260711_gse235863_timepoint_prjna932556")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

soft_gz <- file.path(project_dir, "00_metadata", "source_metadata", "GSE235863_family.soft.gz")
runinfo_path <- file.path(project_dir, "00_metadata", "source_metadata", "PRJNA932556_sra_runinfo.csv")
response_map_path <- file.path(project_dir, "02_input", "response_metadata", "PRJNA932556_sample_response_mapping.tsv")

cd45_major_path <- file.path(prior_dir, "gse235863_cd45_response_celltype_counts.tsv")
cd45_sub_path <- file.path(prior_dir, "gse235863_cd45_response_subcluster_counts.tsv")
cd8_sub_path <- file.path(prior_dir, "gse235863_cd8_response_subcluster_counts.tsv")
cd8_doublet_path <- file.path(prior_dir, "gse235863_cd8_doublet_sensitivity_counts.tsv")

write_tsv <- function(x, path) {
  fwrite(x, path, sep = "\t", na = "", quote = FALSE)
}

derive_timepoint <- function(sample) {
  out <- rep(NA_character_, length(sample))
  out[grepl("(^|-)pre-", sample)] <- "pre-treatment"
  out[grepl("(^|-)post-", sample)] <- "post-treatment"
  out[grepl("(^|-)af-", sample)] <- "post-treatment"
  out
}

derive_tissue <- function(tissue_code) {
  fifelse(tissue_code == "P", "blood",
    fifelse(tissue_code == "T", "liver tumor",
      fifelse(tissue_code == "N", "non-tumor/normal", NA_character_)))
}

derive_title_timepoint <- function(title) {
  out <- rep(NA_character_, length(title))
  out[grepl("^Pre-treatment", title)] <- "pre-treatment"
  out[grepl("^Post-treatment", title)] <- "post-treatment"
  out
}

derive_response_group <- function(raw_response) {
  fifelse(grepl("non-responder|\\bNR\\b|\\bSD\\b|\\bPD\\b", raw_response, ignore.case = TRUE),
    "non_responder",
    fifelse(grepl("responder|\\bPR\\b|\\bCR\\b", raw_response, ignore.case = TRUE),
      "responder", NA_character_))
}

derive_major_from_subcluster <- function(sub_cluster) {
  prefix <- sub("_.*$", "", sub_cluster)
  fifelse(prefix == "CD4", "CD4T",
    fifelse(prefix == "CD8", "CD8T",
      fifelse(prefix %in% c("B", "ILC", "Myeloid"), prefix, prefix)))
}

normalize_count_table <- function(dt, label_col, include_major = TRUE) {
  dt <- copy(dt)
  dt[, derived_timepoint := derive_timepoint(sample)]
  dt[, derived_tissue := derive_tissue(tissue)]
  dt[, label_type := label_col]
  if (label_col == "major_cluster") {
    dt[, sub_cluster := major_cluster]
  } else {
    dt[, major_cluster := derive_major_from_subcluster(sub_cluster)]
  }
  dt[, cell_label := get(label_col)]
  totals <- dt[, .(sample_total_cells = sum(cell_count)), by = .(sample)]
  dt <- merge(dt, totals, by = "sample", all.x = TRUE)
  dt[, proportion_within_sample := cell_count / sample_total_cells]
  keep <- c(
    "dataset_id", "sample", "patient", "raw_response", "derived_response_group",
    "tissue", "derived_tissue", "derived_timepoint", "major_cluster", "sub_cluster",
    "label_type", "cell_label", "cell_count", "sample_total_cells",
    "proportion_within_sample"
  )
  dt[, ..keep]
}

parse_soft <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  sample_starts <- grep("^\\^SAMPLE = ", lines)
  out <- vector("list", length(sample_starts))
  for (i in seq_along(sample_starts)) {
    start <- sample_starts[i]
    end <- if (i < length(sample_starts)) sample_starts[i + 1] - 1 else length(lines)
    block <- lines[start:end]
    gsm <- sub("^\\^SAMPLE = ", "", block[1])
    title <- sub("^!Sample_title = ", "", block[grep("^!Sample_title = ", block)[1]])
    chars <- sub("^!Sample_characteristics_ch1 = ", "", block[grep("^!Sample_characteristics_ch1 = ", block)])
    get_char <- function(key) {
      hit <- chars[grepl(paste0("^", key, ": "), chars)]
      if (length(hit) == 0) return(NA_character_)
      sub(paste0("^", key, ": "), "", hit[1])
    }
    out[[i]] <- data.table(
      source_gsm = gsm,
      source_sample_title = title,
      source_tissue_raw = get_char("tissue"),
      source_patient_raw = get_char("patient"),
      source_response_raw = get_char("response"),
      source_cell_type_raw = get_char("cell type"),
      source_library_type_raw = get_char("library type")
    )
  }
  rbindlist(out, fill = TRUE)
}

soft <- parse_soft(soft_gz)
soft[, title_patient := sub(".* patient (P[0-9]+) .*", "\\1", source_sample_title)]
soft[!grepl(" patient P[0-9]+ ", source_sample_title), title_patient := NA_character_]
soft[, title_timepoint := derive_title_timepoint(source_sample_title)]
soft[, derived_response_group := derive_response_group(source_response_raw)]

cd45_major <- fread(cd45_major_path)
cd45_sub <- fread(cd45_sub_path)
cd8_sub <- fread(cd8_sub_path)
cd8_doublet <- fread(cd8_doublet_path)
cd45_major[, dataset_id := "GSE235863_CD45"]
cd45_sub[, dataset_id := "GSE235863_CD45"]
cd8_sub[, dataset_id := "GSE235863_CD8"]
cd8_doublet[, dataset_id := "GSE235863_CD8"]

cd45_major_norm <- normalize_count_table(cd45_major, "major_cluster")
cd45_sub_norm <- normalize_count_table(cd45_sub, "sub_cluster")
cd45_all <- rbindlist(list(cd45_major_norm, cd45_sub_norm), use.names = TRUE)

cd8_sub_norm <- normalize_count_table(cd8_sub, "sub_cluster")
cd8_doublet_norm <- copy(cd8_doublet)
cd8_doublet_norm[, derived_timepoint := derive_timepoint(sample)]
cd8_doublet_norm[, derived_tissue := derive_tissue(tissue)]
cd8_doublet_norm[, major_cluster := derive_major_from_subcluster(sub_cluster)]
cd8_doublet_norm[, cell_label := sub_cluster]

make_mapping <- function(dt, dataset_id) {
  samples <- unique(dt[, .(dataset_id, object_sample = sample, object_patient_raw = patient, tissue, raw_response)])
  samples[, derived_timepoint := derive_timepoint(object_sample)]
  samples[, derived_tissue := derive_tissue(tissue)]
  samples[, derived_patient := object_patient_raw]
  samples[, derived_response_group := derive_response_group(raw_response)]
  samples[, object_tissue_code := tissue]
  samples[, source_gsm := NA_character_]
  samples[, source_sample_title := NA_character_]
  samples[, source_patient_raw := NA_character_]
  samples[, source_response_raw := NA_character_]
  samples[, source_tissue_raw := NA_character_]
  samples[, mapping_status := NA_character_]
  samples[, mapping_evidence := NA_character_]
  samples[, exclusion_reason := NA_character_]

  for (i in seq_len(nrow(samples))) {
    row <- samples[i]
    hits <- soft[
      title_patient == row$derived_patient &
        title_timepoint == row$derived_timepoint &
        source_tissue_raw == row$derived_tissue
    ]
    if (nrow(hits) >= 1) {
      hit <- hits[1]
      samples[i, source_gsm := hit$source_gsm]
      samples[i, source_sample_title := hit$source_sample_title]
      samples[i, source_patient_raw := hit$source_patient_raw]
      samples[i, source_response_raw := hit$source_response_raw]
      samples[i, source_tissue_raw := hit$source_tissue_raw]
      if (!identical(hit$source_patient_raw, hit$title_patient) ||
          !identical(hit$source_patient_raw, row$derived_patient)) {
        samples[i, mapping_status := "unresolved_excluded"]
        samples[i, mapping_evidence := paste0(
          "GEO title patient=", hit$title_patient,
          "; characteristics patient=", hit$source_patient_raw,
          "; object patient=", row$derived_patient
        )]
        samples[i, exclusion_reason := "conflicting_source_fields_unresolved"]
      } else {
        samples[i, mapping_status := "exact_source_match"]
        samples[i, mapping_evidence := paste0(
          "matched GEO title, tissue, timepoint, and characteristics patient: ",
          hit$source_gsm
        )]
      }
    } else if (!is.na(row$derived_timepoint) && !is.na(row$derived_tissue)) {
      samples[i, mapping_status := "object_only_timepoint_tissue_resolved"]
      samples[i, mapping_evidence := "object sample name and object tissue code resolved timepoint and tissue; no exact GEO sample block matched"]
    } else {
      samples[i, mapping_status := "unresolved_excluded"]
      samples[i, mapping_evidence := "object sample name or tissue code did not resolve to requested pre/post blood/liver tumor strata"]
      samples[i, exclusion_reason := "unresolved_tissue_or_timepoint"]
    }
  }

  samples[, .(
    dataset_id, object_sample, object_patient_raw, source_gsm, source_sample_title,
    source_patient_raw, source_response_raw, source_tissue_raw, derived_patient,
    derived_response_group, derived_tissue, derived_timepoint, mapping_status,
    mapping_evidence, exclusion_reason
  )]
}

mapping <- rbindlist(list(
  make_mapping(cd45_major, "GSE235863_CD45"),
  make_mapping(cd8_sub, "GSE235863_CD8")
), use.names = TRUE, fill = TRUE)
write_tsv(mapping, file.path(out_dir, "gse235863_sample_mapping_audit.tsv"))

valid_samples <- mapping[mapping_status %in% c("exact_source_match", "object_only_timepoint_tissue_resolved") &
  derived_tissue %in% c("blood", "liver tumor") &
  derived_timepoint %in% c("pre-treatment", "post-treatment"),
  .(dataset_id, sample = object_sample, mapping_status)]

cd45_valid <- merge(cd45_all, valid_samples[dataset_id == "GSE235863_CD45"], by = c("dataset_id", "sample"))
cd45_sample_level <- cd45_valid[order(patient, derived_tissue, derived_timepoint, label_type, cell_label)]
write_tsv(cd45_sample_level, file.path(out_dir, "gse235863_cd45_sample_level_composition.tsv"))

cd45_patient_timepoint <- cd45_sample_level[, .(
  cell_count = sum(cell_count),
  total_cells = sum(sample_total_cells),
  proportion = sum(cell_count) / sum(sample_total_cells),
  contributing_samples = paste(sort(unique(sample)), collapse = ";")
), by = .(
  dataset_id, patient, raw_response, derived_response_group, derived_tissue,
  derived_timepoint, major_cluster, sub_cluster, label_type, cell_label
)]
setorder(cd45_patient_timepoint, patient, derived_tissue, derived_timepoint, label_type, cell_label)
write_tsv(cd45_patient_timepoint, file.path(out_dir, "gse235863_cd45_patient_timepoint_composition.tsv"))
write_tsv(
  cd45_patient_timepoint[label_type == "sub_cluster"],
  file.path(out_dir, "gse235863_cd45_patient_timepoint_subcluster.tsv")
)

cd8_valid_samples <- valid_samples[dataset_id == "GSE235863_CD8"]
cd8_sensitivity <- merge(cd8_doublet_norm, cd8_valid_samples, by = c("dataset_id", "sample"))
cd8_sensitivity <- cd8_sensitivity[, .(
  dataset_id, sample, patient, raw_response, derived_response_group,
  tissue, derived_tissue, derived_timepoint, doublet_filter_state,
  major_cluster, sub_cluster, predicted_doublets, cell_count,
  proportion_within_group, total_cells_before_filter,
  retained_cells_after_filter, excluded_predicted_doublets
)]
setorder(cd8_sensitivity, patient, derived_tissue, derived_timepoint, doublet_filter_state, sub_cluster)
write_tsv(cd8_sensitivity, file.path(out_dir, "gse235863_cd8_sample_level_doublet_sensitivity.tsv"))

cd8_patient <- cd8_sensitivity[, .(
  cell_count = sum(cell_count),
  total_cells_before_filter = max(total_cells_before_filter, na.rm = TRUE),
  retained_cells_after_filter = max(retained_cells_after_filter, na.rm = TRUE),
  excluded_predicted_doublets = max(excluded_predicted_doublets, na.rm = TRUE),
  contributing_samples = paste(sort(unique(sample)), collapse = ";")
), by = .(
  dataset_id, patient, raw_response, derived_response_group, derived_tissue,
  derived_timepoint, doublet_filter_state, major_cluster, sub_cluster
)]
cd8_patient_totals <- cd8_patient[, .(state_total_cells = sum(cell_count)),
  by = .(dataset_id, patient, derived_tissue, derived_timepoint, doublet_filter_state)]
cd8_patient <- merge(cd8_patient, cd8_patient_totals,
  by = c("dataset_id", "patient", "derived_tissue", "derived_timepoint", "doublet_filter_state"))
cd8_patient[, proportion := cell_count / state_total_cells]
setorder(cd8_patient, patient, derived_tissue, derived_timepoint, doublet_filter_state, sub_cluster)
write_tsv(cd8_patient, file.path(out_dir, "gse235863_cd8_patient_timepoint_doublet_sensitivity.tsv"))

make_stratum_counts <- function(dt) {
  single <- dt[label_type == "major_cluster", .(
    n_patients_total = uniqueN(patient),
    n_responder = uniqueN(patient[derived_response_group == "responder"]),
    n_non_responder = uniqueN(patient[derived_response_group == "non_responder"]),
    n_samples = uniqueN(contributing_samples)
  ), by = .(derived_tissue, derived_timepoint)]
  single[, stratum := paste(derived_timepoint, derived_tissue)]
  paired <- dt[label_type == "major_cluster",
    .(has_pre = any(derived_timepoint == "pre-treatment"),
      has_post = any(derived_timepoint == "post-treatment"),
      derived_response_group = first(derived_response_group)),
    by = .(derived_tissue, patient)
  ][has_pre & has_post, .(
    n_patients_total = .N,
    n_responder = sum(derived_response_group == "responder"),
    n_non_responder = sum(derived_response_group == "non_responder"),
    n_samples = NA_integer_
  ), by = .(derived_tissue)]
  paired[, derived_timepoint := "paired pre/post"]
  paired[, stratum := paste("paired pre/post within", derived_tissue)]
  rbindlist(list(
    single[, .(stratum, derived_tissue, derived_timepoint, n_patients_total, n_responder, n_non_responder, n_samples)],
    paired[, .(stratum, derived_tissue, derived_timepoint, n_patients_total, n_responder, n_non_responder, n_samples)]
  ), use.names = TRUE, fill = TRUE)
}

strata <- make_stratum_counts(cd45_patient_timepoint)
write_tsv(strata, file.path(out_dir, "gse235863_analysis_stratum_counts.tsv"))

make_stats <- function(dt) {
  base <- dt[, .(
    dataset_id, patient, derived_response_group, derived_tissue, derived_timepoint,
    major_cluster, sub_cluster, label_type, cell_label, proportion
  )]
  single <- base[, {
    r <- proportion[derived_response_group == "responder"]
    n <- proportion[derived_response_group == "non_responder"]
    nr <- uniqueN(patient[derived_response_group == "responder"])
    nn <- uniqueN(patient[derived_response_group == "non_responder"])
    status <- if (nr < 3 || nn < 3) "descriptive_only_group_n_below_3" else "inferential_test_completed"
    p <- if (status == "inferential_test_completed") wilcox.test(r, n, exact = FALSE)$p.value else NA_real_
    .(
      n_patients_responder = nr,
      n_patients_non_responder = nn,
      median_responder = if (length(r)) median(r) else NA_real_,
      median_non_responder = if (length(n)) median(n) else NA_real_,
      median_difference_responder_minus_non_responder = if (length(r) && length(n)) median(r) - median(n) else NA_real_,
      exact_test_name = if (status == "inferential_test_completed") "wilcoxon_rank_sum_exact_false" else "none_descriptive_only",
      raw_p = p,
      analysis_status = status
    )
  }, by = .(stratum = paste(derived_timepoint, derived_tissue), derived_tissue, derived_timepoint, major_cluster, sub_cluster, label_type, cell_label)]

  wide <- dcast(base, patient + derived_response_group + derived_tissue + major_cluster + sub_cluster + label_type + cell_label ~ derived_timepoint,
    value.var = "proportion")
  if (!("pre-treatment" %in% names(wide))) wide[, `pre-treatment` := NA_real_]
  if (!("post-treatment" %in% names(wide))) wide[, `post-treatment` := NA_real_]
  wide <- wide[!is.na(`pre-treatment`) & !is.na(`post-treatment`)]
  wide[, delta_post_minus_pre := `post-treatment` - `pre-treatment`]
  paired <- wide[, {
    r <- delta_post_minus_pre[derived_response_group == "responder"]
    n <- delta_post_minus_pre[derived_response_group == "non_responder"]
    nr <- uniqueN(patient[derived_response_group == "responder"])
    nn <- uniqueN(patient[derived_response_group == "non_responder"])
    status <- if (nr < 3 || nn < 3) "descriptive_only_group_n_below_3" else "inferential_test_completed"
    p <- if (status == "inferential_test_completed") wilcox.test(r, n, exact = FALSE)$p.value else NA_real_
    .(
      n_patients_responder = nr,
      n_patients_non_responder = nn,
      median_responder = if (length(r)) median(r) else NA_real_,
      median_non_responder = if (length(n)) median(n) else NA_real_,
      median_difference_responder_minus_non_responder = if (length(r) && length(n)) median(r) - median(n) else NA_real_,
      exact_test_name = if (status == "inferential_test_completed") "wilcoxon_rank_sum_exact_false" else "none_descriptive_only",
      raw_p = p,
      analysis_status = status
    )
  }, by = .(stratum = paste("paired pre/post within", derived_tissue), derived_tissue, major_cluster, sub_cluster, label_type, cell_label)]
  paired[, derived_timepoint := "paired pre/post"]
  out <- rbindlist(list(single, paired), use.names = TRUE, fill = TRUE)
  out[, adjusted_p := ifelse(is.na(raw_p), NA_real_, p.adjust(raw_p, method = "BH"))]
  setcolorder(out, c(
    "stratum", "derived_tissue", "derived_timepoint", "major_cluster", "sub_cluster",
    "label_type", "cell_label", "n_patients_responder", "n_patients_non_responder",
    "median_responder", "median_non_responder",
    "median_difference_responder_minus_non_responder", "exact_test_name",
    "raw_p", "adjusted_p", "analysis_status"
  ))
  setorder(out, stratum, label_type, cell_label)
  out
}

stats <- make_stats(cd45_patient_timepoint)
write_tsv(stats, file.path(out_dir, "gse235863_patient_level_composition_statistics.tsv"))

plot_cd45 <- cd45_patient_timepoint[label_type == "major_cluster"]
plot_cd45[, stratum := paste(derived_timepoint, derived_tissue)]
p1 <- ggplot(plot_cd45, aes(x = major_cluster, y = proportion, fill = derived_timepoint)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.82, position = position_dodge(width = 0.78)) +
  geom_point(aes(shape = derived_response_group), position = position_jitterdodge(jitter.width = 0.16, dodge.width = 0.78), size = 1.7, alpha = 0.9) +
  facet_grid(derived_tissue ~ derived_response_group) +
  scale_fill_manual(values = c("pre-treatment" = "#0072B2", "post-treatment" = "#D55E00")) +
  scale_shape_manual(values = c("responder" = 16, "non_responder" = 17)) +
  labs(x = NULL, y = "Proportion within patient/timepoint", fill = "Timepoint", shape = "Response") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(out_dir, "gse235863_cd45_timepoint_composition_qc.png"), p1, width = 11, height = 7.2, dpi = 220)

top_cd8 <- cd8_sensitivity[doublet_filter_state == "after_excluding_predicted_doublets",
  .(cells = sum(cell_count)), by = sub_cluster][order(-cells)][1:min(.N, 14), sub_cluster]
plot_cd8 <- cd8_sensitivity[sub_cluster %in% top_cd8]
p2 <- ggplot(plot_cd8, aes(x = sub_cluster, y = proportion_within_group, fill = doublet_filter_state)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68, alpha = 0.9) +
  facet_grid(derived_tissue + derived_timepoint ~ derived_response_group, scales = "free_y") +
  scale_fill_manual(values = c("before_doublet_exclusion" = "#009E73", "after_excluding_predicted_doublets" = "#CC79A7")) +
  labs(x = NULL, y = "Proportion", fill = "Doublet filter") +
  theme_bw(base_size = 10.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(out_dir, "gse235863_cd8_timepoint_doublet_qc.png"), p2, width = 12.5, height = 8.4, dpi = 220)

runinfo <- fread(runinfo_path)
response_map <- fread(response_map_path)
tool_names <- c("prefetch", "fasterq-dump", "fastq-dump", "cellranger", "kallisto", "bustools", "salmon", "STARsolo")
tool_paths <- Sys.which(tool_names)
tool_dt <- data.table(tool = names(tool_paths), path = unname(tool_paths))
tool_dt[, available := nzchar(path)]

pilot <- merge(runinfo, response_map[, .(sample_name, response_group, response_label, run_accession)],
  by.x = "SampleName", by.y = "sample_name", all.x = TRUE)
pilot[, observed_read_files := "not_observed_no_sra_fastq_tool_available"]
pilot[, observed_read_lengths := as.character(avgLength)]
pilot[, observed_pairing := LibraryLayout]
pilot[, observed_headers := "not_observed_no_fastq_generated"]
pilot[, barcode_umi_evidence := "not_validated"]
pilot[, cdna_evidence := "not_validated"]
pilot[, chemistry := "not_validated"]
pilot[, tool := "none_available_on_path"]
pilot[, reference := AssemblyName]
pilot[, decision := "blocked_read_structure_not_validated"]
pilot[, decision_reason := "SRA metadata alone is insufficient to validate single-cell barcode/UMI/cDNA read structure; required SRA or quantification tools were not available on PATH"]
pilot_audit <- pilot[, .(
  Run, SampleName, response_group, response_label, size_MB, avgLength, LibraryLayout,
  spots_with_mates, download_path, observed_read_files, observed_read_lengths,
  observed_pairing, observed_headers, barcode_umi_evidence, cdna_evidence,
  chemistry, tool, reference, decision, decision_reason
)]
write_tsv(pilot_audit, file.path(out_dir, "prjna932556_read_structure_audit.tsv"))

fastq_metrics <- pilot[, .(
  Run, SampleName, fastq_file = NA_character_, read_role = NA_character_,
  n_reads_checked = NA_integer_, read_length_min = NA_integer_,
  read_length_median = NA_integer_, read_length_max = NA_integer_,
  barcode_umi_pattern = "not_available_no_fastq_generated",
  metric_status = "not_computed_read_structure_not_validated"
)]
write_tsv(fastq_metrics, file.path(out_dir, "prjna932556_pilot_fastq_metrics.tsv"))

feasibility <- data.table(
  resource_id = "PRJNA932556",
  samples_requested = nrow(runinfo),
  total_size_MB = sum(as.numeric(runinfo$size_MB)),
  sra_tools_available = any(tool_dt$available[tool_dt$tool %in% c("prefetch", "fasterq-dump", "fastq-dump")]),
  quantification_tools_available = any(tool_dt$available[tool_dt$tool %in% c("cellranger", "kallisto", "bustools", "salmon", "STARsolo")]),
  read_structure_validated = FALSE,
  decision = "blocked_read_structure_not_validated",
  full_processing_status = "not_started_stop_condition_met",
  reason = "pilot could not observe FASTQ read files, headers, barcode/UMI evidence, cDNA read evidence, or chemistry"
)
write_tsv(feasibility, file.path(out_dir, "prjna932556_quantification_feasibility.tsv"))

command_log <- c(
  "PRJNA932556 pilot command log",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Tool availability checked with Sys.which in R:",
  paste(tool_dt$tool, ifelse(tool_dt$available, tool_dt$path, "NOT_FOUND_ON_PATH"), sep = "\t"),
  "",
  "No SRA prefetch/fasterq-dump/fastq-dump executable was available on PATH.",
  "No Cell Ranger, kallisto/bustools, salmon, or STARsolo executable was available on PATH.",
  "No FASTQ headers or read sequences were observed.",
  "Stopped per request because barcode/UMI structure and chemistry were not validated."
)
writeLines(command_log, file.path(out_dir, "prjna932556_pilot_command_log.txt"))

run_status <- data.table(
  stage = c(
    "GSE235863_CD45_timepoint_rebuild",
    "GSE235863_CD8_doublet_timepoint_rebuild",
    "PRJNA932556_pilot_read_structure",
    "PRJNA932556_full_quantification"
  ),
  status = c(
    "completed",
    "completed_with_unresolved_samples_excluded",
    "completed_stop_condition_met",
    "not_started_stop_condition_met"
  ),
  detail = c(
    "CD45 sample, patient-timepoint, subcluster, stratum, statistics, and QC plot outputs written",
    "CD8 doublet sensitivity outputs written after excluding samples with unresolved tissue or timepoint",
    "Pilot audit written; FASTQ/barcode/UMI/cDNA/chemistry not validated",
    "Full six-sample matrices and Seurat QC object were not generated because pilot failed validation"
  )
)
write_tsv(run_status, file.path(out_dir, "run_status.tsv"))

warning_log <- data.table(
  item = c(
    "GSM7510911",
    "CD8 unresolved tissue/timepoint samples",
    "PRJNA932556 full processing outputs"
  ),
  severity = c("high", "medium", "high"),
  message = c(
    "GEO sample title is Pre-treatment blood of patient P27 (R/PR) but characteristics_ch1 patient is P18; affected P27-pre-P sample excluded from analytical mapping",
    paste(mapping[dataset_id == "GSE235863_CD8" & mapping_status == "unresolved_excluded", object_sample], collapse = ";"),
    "Not generated because read structure/chemistry was not validated in pilot"
  )
)
write_tsv(warning_log, file.path(out_dir, "warning_log.tsv"))

all_outputs <- list.files(out_dir, full.names = TRUE, recursive = FALSE)
all_outputs <- all_outputs[basename(all_outputs) != "object_inventory.tsv"]
inventory <- data.table(
  file_name = basename(all_outputs),
  relative_path = file.path("bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260711_gse235863_timepoint_prjna932556", basename(all_outputs)),
  size_bytes = file.info(all_outputs)$size,
  role = fifelse(grepl("^gse235863", basename(all_outputs)), "GSE235863 timepoint rebuild",
    fifelse(grepl("^prjna932556", basename(all_outputs)), "PRJNA932556 pilot audit", "run metadata"))
)
write_tsv(inventory, file.path(out_dir, "object_inventory.tsv"))

status_md <- c(
  "# 20260711 GSE235863 timepoint rebuild and PRJNA932556 pilot",
  "",
  "## Status",
  "",
  "- GSE235863 CD45: completed timepoint-preserving sample and patient-level composition rebuild.",
  "- GSE235863 CD8: completed timepoint-preserving doublet sensitivity tables after excluding unresolved tissue/timepoint samples.",
  "- PRJNA932556: stopped after pilot with `blocked_read_structure_not_validated`; no full count matrices were generated.",
  "",
  "## Key audit decisions",
  "",
  "- `GSM7510911` raw GEO fields were preserved. Its title indicates `P27`, while `characteristics_ch1::patient` is `P18`; `P27-pre-P` was excluded from analytical mapping.",
  "- Pre/post and blood/liver tumor strata were not collapsed.",
  "- SRA metadata alone was not used to infer 10x chemistry or force single-cell quantification.",
  "",
  "## Output directory",
  "",
  "`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260711_gse235863_timepoint_prjna932556`"
)
writeLines(status_md, file.path(out_dir, "STATUS.md"))

session_file <- file.path(out_dir, "session_info.txt")
zz <- file(session_file, open = "wt")
sink(zz)
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("R executable: ", file.path(R.home("bin"), "R"), "\n", sep = "")
cat("Project dir: ", project_dir, "\n\n", sep = "")
sessionInfo()
sink()
close(zz)

files_for_sha <- list.files(out_dir, full.names = TRUE, recursive = FALSE)
files_for_sha <- files_for_sha[basename(files_for_sha) != "SHA256SUMS.txt"]
sha_values <- vapply(files_for_sha, function(path) {
  cmd <- sprintf(
    "Get-FileHash -LiteralPath %s -Algorithm SHA256 | Select-Object -ExpandProperty Hash",
    shQuote(path, type = "cmd")
  )
  out <- system2("powershell", c("-NoProfile", "-Command", cmd), stdout = TRUE, stderr = TRUE)
  tolower(out[1])
}, character(1))
sha_dt <- data.table(sha256 = unname(sha_values), file_name = basename(names(sha_values)))
write_tsv(sha_dt, file.path(out_dir, "SHA256SUMS.txt"))

message("Wrote outputs to: ", out_dir)
