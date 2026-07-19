#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
output_dir <- file.path(project_dir, "desktop_exchange/uploads/20260717_gse205506_paper_qc_rebuild")

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

log_step <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
  flush.console()
}

required_files <- c(
  "gse205506_paper_qc_parameters.tsv",
  "gse205506_paper_qc_cell_audit.tsv.gz",
  "gse205506_paper_qc_stage_counts.tsv",
  "gse205506_paper_qc_sample_counts.tsv",
  "gse205506_paper_qc_compartment_mt_model.tsv",
  "gse205506_paper_qc_removal_fraction_comparison.tsv",
  "gse205506_paper_qc_doublet_sensitivity.tsv",
  "gse205506_paper_qc_broad_cluster_evidence.tsv",
  "gse205506_paper_qc_hybrid_review.tsv",
  "gse205506_paper_qc_final_cell_ids.tsv.gz",
  "gse205506_paper_qc_final_object_summary.tsv",
  "gse205506_paper_qc_final_counts.tsv",
  "gse205506_paper_qc_final_broad_cluster_evidence.tsv",
  "gse205506_paper_qc_batch_knn_audit.tsv",
  "gse205506_paper_qc_batch_pc_eta_squared.tsv",
  "gse205506_paper_qc_final_all_cell_top10_markers.tsv",
  "gse205506_paper_qc_final_all_cell_top50_markers.tsv",
  "gse205506_paper_qc_compartment_object_inventory.tsv",
  "gse205506_paper_qc_compartment_resolution_summary.tsv",
  "gse205506_paper_qc_package_versions.tsv",
  "gse205506_paper_qc_basic_object_summary.tsv",
  "gse205506_paper_qc_first_round_object_inventory.tsv",
  "gse205506_paper_qc_figure_visual_qc.tsv"
)
compartment_slugs <- c("T_I_NK", "B", "Myeloid", "Endothelial", "Fibroblast", "Epithelial")
for (slug in compartment_slugs) {
  required_files <- c(
    required_files,
    paste0("gse205506_paper_qc_", slug, "_top10_markers.tsv"),
    paste0("gse205506_paper_qc_", slug, "_top50_markers.tsv"),
    paste0("gse205506_paper_qc_", slug, "_author_overlap.tsv"),
    paste0("gse205506_paper_qc_", slug, "_annotation_evidence.tsv")
  )
}
missing_files <- required_files[!file.exists(file.path(output_dir, required_files))]
if (length(missing_files)) stop("Required handoff files are missing: ", paste(missing_files, collapse = "; "), call. = FALSE)

figure_stems <- c(
  "gse205506_paper_qc_umap_broad",
  "gse205506_paper_qc_umap_response",
  "gse205506_paper_qc_umap_timepoint",
  "gse205506_paper_qc_umap_sample",
  "gse205506_paper_qc_before_after",
  "gse205506_paper_qc_marker_dotplot_broad",
  paste0("gse205506_paper_qc_umap_", compartment_slugs[1:5]),
  paste0("gse205506_paper_qc_marker_dotplot_", compartment_slugs[1:5])
)
figure_files <- as.vector(outer(figure_stems, c("png", "pdf"), paste, sep = "."))
missing_figures <- figure_files[!file.exists(file.path(output_dir, "figures", figure_files))]
if (length(missing_figures)) stop("Required figures are missing: ", paste(missing_figures, collapse = "; "), call. = FALSE)

visual_qc <- fread(file.path(output_dir, "gse205506_paper_qc_figure_visual_qc.tsv"))
expected_visual_columns <- c("figure", "inspected", "clipping", "text_overlap", "legend_overflow", "point_washout", "color_separation", "result")
if (!all(expected_visual_columns %in% names(visual_qc))) stop("Figure visual-QC schema is incomplete", call. = FALSE)
if (nrow(visual_qc) != length(figure_stems) || !setequal(visual_qc$figure, paste0(figure_stems, ".png"))) {
  stop("Figure visual-QC rows do not match the required PNG set", call. = FALSE)
}
if (!all(visual_qc$inspected %in% TRUE) || !all(visual_qc$result == "PASS")) {
  stop("At least one required PNG has not passed visual inspection", call. = FALSE)
}

basic <- fread(file.path(output_dir, "gse205506_paper_qc_basic_object_summary.tsv"))
first_round <- fread(file.path(output_dir, "gse205506_paper_qc_first_round_object_inventory.tsv"))
final_summary <- fread(file.path(output_dir, "gse205506_paper_qc_final_object_summary.tsv"))
compartment_inventory <- fread(file.path(output_dir, "gse205506_paper_qc_compartment_object_inventory.tsv"))
rejected <- fread(file.path(output_dir, "gse205506_rejected_qc_v1_object_inventory.tsv"))

local_inventory <- rbindlist(list(
  basic[, .(
    object_type = "article_basic_qc_primary",
    compartment = NA_character_,
    absolute_path = object_path,
    size_bytes,
    sha256,
    cells,
    status = "formal_route_local_only",
    notes = "Article numerical basic QC; scDblFinder calls retained in primary object"
  )],
  first_round[, .(
    object_type = fifelse(grepl("first_round_rpca", object_path), "first_round_rpca_broad_classification", "post_compartment_mt_qc_raw"),
    compartment = NA_character_,
    absolute_path = object_path,
    size_bytes,
    sha256,
    cells,
    status = "formal_route_local_only",
    notes = fifelse(grepl("first_round_rpca", object_path), "First-round integrated object used for broad classification and mitochondrial QC", "Raw RNA object restricted to cells passing compartment mitochondrial QC")
  )],
  final_summary[, .(
    object_type = "formal_post_mt_rpca_atlas",
    compartment = NA_character_,
    absolute_path = object_path,
    size_bytes,
    sha256,
    cells,
    status = "formal_route_local_only",
    notes = "Final post-mitochondrial-QC RPCA atlas"
  )],
  compartment_inventory[, .(
    object_type = "compartment_reclustered_object",
    compartment,
    absolute_path = object_path,
    size_bytes,
    sha256,
    cells,
    status = "formal_route_local_only",
    notes = "Compartment-specific Harmony project implementation; no whole cluster deleted"
  )],
  rejected[, .(
    object_type = "rejected_qc_version_1",
    compartment = NA_character_,
    absolute_path,
    size_bytes,
    sha256,
    cells = NA_integer_,
    status,
    notes = reason
  )]
), use.names = TRUE, fill = TRUE)
write_tsv(local_inventory, "gse205506_paper_qc_local_object_inventory.tsv")

stage_counts <- fread(file.path(output_dir, "gse205506_paper_qc_stage_counts.tsv"))
mt_model <- fread(file.path(output_dir, "gse205506_paper_qc_compartment_mt_model.tsv"))
first_evidence <- fread(file.path(output_dir, "gse205506_paper_qc_broad_cluster_evidence.tsv"))
final_evidence <- fread(file.path(output_dir, "gse205506_paper_qc_final_broad_cluster_evidence.tsv"))

raw_cells <- stage_counts[stage == "raw_author_matrices", cells][[1L]]
basic_cells <- stage_counts[stage == "paper_basic_qc_pass", cells][[1L]]
mt_pass_cells <- stage_counts[stage == "compartment_mt_qc_pass", cells][[1L]]
final_cells <- final_summary$cells[[1L]]
if (mt_pass_cells != final_cells) stop("Final atlas cell count does not match the compartment-mt-QC pass count", call. = FALSE)
if (final_summary$samples[[1L]] != 40L || final_summary$subjects[[1L]] != 19L) {
  stop("Final atlas does not preserve all 40 samples and 19 subjects", call. = FALSE)
}
if (uniqueN(first_evidence$assigned_broad_compartment) != 6L || uniqueN(final_evidence$assigned_broad_compartment) != 6L) {
  stop("First-round or final broad evidence does not support exactly six compartments", call. = FALSE)
}

unresolved_rows <- rbindlist(lapply(compartment_slugs[1:5], function(slug) {
  evidence <- fread(file.path(output_dir, paste0("gse205506_paper_qc_", slug, "_annotation_evidence.tsv")))
  data.table(
    compartment = unique(evidence$compartment),
    clusters = nrow(evidence),
    resolved_author_labels = sum(evidence$final_reviewed_label != "Unresolved_review"),
    unresolved_review_clusters = sum(evidence$final_reviewed_label == "Unresolved_review")
  )
}))
write_tsv(unresolved_rows, "gse205506_paper_qc_compartment_review_summary.tsv")

mt_lines <- vapply(seq_len(nrow(mt_model)), function(index) {
  sprintf(
    "- %s: %d/%d removed (%.2f%%); article comparison %.2f%%.",
    mt_model$model_group[[index]],
    mt_model$removed_cells[[index]],
    mt_model$cells[[index]],
    100 * mt_model$removal_fraction[[index]],
    100 * c(Lymphoid = 0.0920, Myeloid = 0.1284, Fibroblast = 0.0811, Endothelial = 0.0850, Epithelial = 0.2975)[mt_model$model_group[[index]]]
  )
}, character(1))

status_lines <- c(
  "# GSE205506 source-paper QC rebuild",
  "",
  "Status: `COMPLETE_GSE205506_PAPER_QC_REBUILD`",
  "",
  "## Formal route",
  "",
  sprintf("- Re-imported all 40 author matrices: %s raw barcodes from 40 samples and 19 subjects.", format(raw_cells, big.mark = ",")),
  sprintf("- Article basic-QC pass: %s cells; scDblFinder was sensitivity-only and excluded zero cells from the primary route.", format(basic_cells, big.mark = ",")),
  sprintf("- Compartment mitochondrial-QC pass and final atlas: %s cells.", format(final_cells, big.mark = ",")),
  sprintf("- Article reported total: 155,397; observed difference: %+d cells. Equality was not forced.", final_cells - 155397L),
  sprintf("- First-round broad clusters: %d; final broad clusters: %d; exactly six broad compartments represented in both.", nrow(first_evidence), nrow(final_evidence)),
  "- No whole broad or subtype cluster was deleted.",
  "",
  "## Compartment mitochondrial QC",
  "",
  mt_lines,
  "",
  "## Compartment review",
  "",
  vapply(seq_len(nrow(unresolved_rows)), function(index) sprintf(
    "- %s: %d clusters; %d author-overlap labels supported; %d retained as `Unresolved_review`.",
    unresolved_rows$compartment[[index]], unresolved_rows$clusters[[index]],
    unresolved_rows$resolved_author_labels[[index]], unresolved_rows$unresolved_review_clusters[[index]]
  ), character(1)),
  "- Epithelial was saved as a separate sixth object and received no mmc3-derived subtype label.",
  "",
  "## Verification",
  "",
  "- All required PNG/PDF figures were generated with the approved bright palette.",
  sprintf("- Visual inspection passed for %d PNG figures; details are in `gse205506_paper_qc_figure_visual_qc.tsv`.", nrow(visual_qc)),
  "- Exact local paths, byte sizes, and SHA-256 hashes are in `gse205506_paper_qc_local_object_inventory.tsv`.",
  "- The prior QC-version-1 RDS files were preserved and remain listed as rejected in `gse205506_rejected_qc_v1_object_inventory.tsv`.",
  ""
)
writeLines(status_lines, file.path(output_dir, "STATUS.md"), useBytes = TRUE)

review_lines <- c(
  "# Mac review request: GSE205506 source-paper QC rebuild",
  "",
  "Please review these items before using this object for CRC Figure 1:",
  "",
  "1. Compare observed mitochondrial-removal fractions with the article values in `gse205506_paper_qc_removal_fraction_comparison.tsv`; parameters were not tuned to the article total.",
  "2. Review first-round and final broad evidence, especially conflicts listed in `gse205506_paper_qc_hybrid_review.tsv` and the two broad-cluster evidence tables.",
  "3. Review every `Unresolved_review` compartment cluster in the six annotation-evidence tables; no author subtype was forced.",
  "4. Review the bright-palette UMAPs and marker dotplots in `figures/`.",
  "5. Confirm the final cell total and sample/subject traceability from `gse205506_paper_qc_final_object_summary.tsv` and `gse205506_paper_qc_final_counts.tsv`.",
  "",
  "The RDS objects are local-only on F: and are fully identified by path, size, and SHA-256 in `gse205506_paper_qc_local_object_inventory.tsv`.",
  ""
)
writeLines(review_lines, file.path(output_dir, "MAC_REVIEW_REQUEST.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_gse205506_handoff_finalization.txt"))
log_step("GSE205506 handoff completion files generated")
