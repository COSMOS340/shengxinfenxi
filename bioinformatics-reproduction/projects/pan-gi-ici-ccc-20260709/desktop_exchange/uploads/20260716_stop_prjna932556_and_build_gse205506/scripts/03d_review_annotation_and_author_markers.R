#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: 03d_review_annotation_and_author_markers.R <gse205506_root> <output_dir>",
    call. = FALSE
  )
}

gse_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(gse_root, "r_analysis")
input_rds <- file.path(analysis_dir, "gse205506_seurat_final_annotated.rds")
reviewed_rds <- file.path(analysis_dir, "gse205506_seurat_final_reviewed.rds")
figure_dir <- file.path(output_dir, "figures")

required_paths <- c(
  input_rds,
  file.path(output_dir, "gse205506_cluster_top10_markers.tsv"),
  file.path(output_dir, "gse205506_cluster_top50_markers.tsv"),
  file.path(output_dir, "gse205506_author_marker_reference_top50.tsv")
)
if (any(!file.exists(required_paths))) {
  stop("Missing required files: ", paste(required_paths[!file.exists(required_paths)], collapse = "; "), call. = FALSE)
}

write_tsv <- function(x, filename) {
  fwrite(x, file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

save_figure <- function(plot, stem, width, height) {
  plot <- plot + theme(plot.margin = margin(10, 14, 10, 10))
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, device = cairo_pdf, bg = "white", limitsize = FALSE
  )
}

review <- data.table(
  cluster = as.character(0:27),
  detailed_annotation = c(
    "Epithelial", "B_cell", "Epithelial_mito_enriched", "CD4_T_cell",
    "Cytotoxic_T_NK", "Cycling_epithelial", "Epithelial_mito_enriched", "Endothelial",
    "Colonocyte", "Myeloid", "T_cell", "Plasma_cell",
    "Erythroid_mito_enriched", "BEST4_OTOP2_epithelial", "Goblet", "Colonocyte_mito_enriched",
    "Fibroblast", "Pericyte_smooth_muscle", "Cycling_immune", "Cycling_epithelial_mito_enriched",
    "Tuft_epithelial", "Immune_mito_enriched", "Lymphatic_endothelial", "Immune_mixed_mito_enriched",
    "Fibroblast_mito_enriched", "Myeloid_IG_enriched", "Epithelial", "Unresolved_mito_enriched"
  ),
  broad_annotation = c(
    "Epithelial", "B cell", "Epithelial", "T/NK",
    "T/NK", "Epithelial", "Epithelial", "Endothelial",
    "Epithelial", "Myeloid", "T/NK", "Plasma",
    "Erythroid", "Epithelial", "Epithelial", "Epithelial",
    "Fibroblast/SMC", "Fibroblast/SMC", "Immune unresolved", "Epithelial",
    "Epithelial", "Immune unresolved", "Endothelial", "Immune unresolved",
    "Fibroblast/SMC", "Myeloid", "Epithelial", "Unresolved"
  ),
  evidence_summary = c(
    "EPCAM;KRT18;TSPAN8;KRT8", "MS4A1;CD79A;CD74;CD37", "epithelial genes after mitochondrial markers", "IL7R;CD2;RORA;KLRB1",
    "CCL5;NKG7;GZMA;CD3D", "MKI67;STMN1;TOP2A;EPCAM", "LGALS4;PIGR;ELF3;EPCAM", "PECAM1;VWF;EGFL7;SPARCL1",
    "FABP1;SLC26A3;CEACAM7;KRT20", "LYZ;TYROBP;FCER1G;LST1", "TRAC;CD3D;CD2;IL32", "JCHAIN;MZB1;XBP1;immunoglobulin genes",
    "ALAS2;HBM;HBD after mitochondrial markers", "BEST4;OTOP2;CA7;CA4", "TFF3;FCGBP;MUC2;SPINK4", "SLC26A3;FABP1;CA1 after mitochondrial markers",
    "COL1A1;COL1A2;DCN;COL3A1", "ACTA2;TAGLN;MYH11;MCAM", "STMN1;CORO1A;LCP1;cycling genes", "ONECUT2;MKI67;CENPF;TOP2A after mitochondrial markers",
    "POU2F3;TRPM5;SH2D6;AVIL;IL17RB", "immunoglobulin and immune genes after mitochondrial markers", "CCL21;LYVE1;PROX1;MMRN1", "T-cell and myeloid genes mixed after mitochondrial markers",
    "COL1A1;COL1A2;DCN;TAGLN after mitochondrial markers", "LYZ;C1QA;HLA-DRA with immunoglobulin enrichment", "PIGR;ELF3;KRT8;EPCAM", "mitochondrial markers dominate; no supported lineage"
  ),
  expected_author_sheet = c(
    "", "B cell markers", "", "T_I_NK cell markers",
    "T_I_NK cell markers", "", "", "Endothelial cell markers",
    "", "Myeloid cell markers", "T_I_NK cell markers", "B cell markers",
    "", "", "", "",
    "Fibroblast markers", "Fibroblast markers", "", "",
    "", "", "Endothelial cell markers", "",
    "Fibroblast markers", "Myeloid cell markers", "", ""
  )
)

top10 <- fread(file.path(output_dir, "gse205506_cluster_top10_markers.tsv"), sep = "\t", header = TRUE, na.strings = NULL)
top50 <- fread(file.path(output_dir, "gse205506_cluster_top50_markers.tsv"), sep = "\t", header = TRUE, na.strings = NULL)
top10[, cluster := as.character(cluster)]
top50[, cluster := as.character(cluster)]
review[, top10_markers := vapply(cluster, function(value) {
  paste(top10[cluster == value][order(rank), gene], collapse = ";")
}, character(1))]
review[, top10_mitochondrial_markers := vapply(cluster, function(value) {
  sum(grepl("^MT-", top10[cluster == value, gene]))
}, integer(1))]
review[, quality_review := ifelse(
  top10_mitochondrial_markers >= 5L,
  "mitochondrial_marker_enriched",
  "retained_after_qc"
)]

object <- readRDS(input_rds)
if (ncol(object) != 238934L || length(unique(object$working_cluster)) != 28L) {
  stop("Annotated object does not match the verified 238934-cell, 28-cluster checkpoint", call. = FALSE)
}

metadata <- as.data.table(object@meta.data, keep.rownames = "cell_barcode")
metadata[, working_cluster_character := as.character(working_cluster)]
cluster_qc <- metadata[, .(
  cells = .N,
  samples = uniqueN(geo_accession),
  subjects = uniqueN(geo_subject),
  median_nCount_RNA = as.numeric(median(nCount_RNA_source_qc)),
  q05_nCount_RNA = as.numeric(quantile(nCount_RNA_source_qc, 0.05, type = 8)),
  q95_nCount_RNA = as.numeric(quantile(nCount_RNA_source_qc, 0.95, type = 8)),
  median_nFeature_RNA = as.numeric(median(nFeature_RNA_source_qc)),
  q05_nFeature_RNA = as.numeric(quantile(nFeature_RNA_source_qc, 0.05, type = 8)),
  q95_nFeature_RNA = as.numeric(quantile(nFeature_RNA_source_qc, 0.95, type = 8)),
  median_percent_mt = as.numeric(median(percent.mt)),
  q75_percent_mt = as.numeric(quantile(percent.mt, 0.75, type = 8)),
  q95_percent_mt = as.numeric(quantile(percent.mt, 0.95, type = 8))
), by = .(cluster = working_cluster_character)]

sample_counts <- metadata[, .N, by = .(cluster = working_cluster_character, geo_accession)]
setorder(sample_counts, cluster, -N, geo_accession)
top_sample <- sample_counts[, .SD[1L], by = cluster]
setnames(top_sample, c("geo_accession", "N"), c("largest_sample", "largest_sample_cells"))
cluster_qc <- merge(cluster_qc, top_sample, by = "cluster", all.x = TRUE, sort = FALSE)
cluster_qc[, largest_sample_fraction := largest_sample_cells / cells]
cluster_qc <- merge(cluster_qc, review, by = "cluster", all.x = TRUE, sort = FALSE)
cluster_qc[, cluster_numeric := as.integer(cluster)]
setorder(cluster_qc, cluster_numeric)
cluster_qc[, cluster_numeric := NULL]
if (anyNA(cluster_qc$detailed_annotation) || nrow(cluster_qc) != 28L) {
  stop("Reviewed annotation table did not join exactly to all 28 clusters", call. = FALSE)
}
write_tsv(cluster_qc, "gse205506_cluster_qc_and_reviewed_annotation.tsv")
write_tsv(cluster_qc, "gse205506_annotation_evidence.tsv")

detailed_map <- setNames(review$detailed_annotation, review$cluster)
broad_map <- setNames(review$broad_annotation, review$cluster)
quality_map <- setNames(review$quality_review, review$cluster)
cluster_values <- as.character(object$working_cluster)
object$cell_type_detailed_reviewed <- unname(detailed_map[cluster_values])
object$cell_type_broad_reviewed <- unname(broad_map[cluster_values])
object$cluster_quality_reviewed <- unname(quality_map[cluster_values])
if (anyNA(object$cell_type_detailed_reviewed) || anyNA(object$cell_type_broad_reviewed)) {
  stop("Reviewed annotation did not map to every cell", call. = FALSE)
}

metadata[, cell_type_detailed_reviewed := unname(detailed_map[working_cluster_character])]
metadata[, cell_type_broad_reviewed := unname(broad_map[working_cluster_character])]
metadata[, cluster_quality_reviewed := unname(quality_map[working_cluster_character])]
patient_counts <- metadata[, .N, by = .(
  geo_subject,
  geo_accession,
  geo_genotype,
  geo_treatment,
  derived_timepoint_from_exact_geo_treatment,
  table_s1_response,
  working_cluster = working_cluster_character,
  cell_type_detailed_reviewed,
  cell_type_broad_reviewed,
  cluster_quality_reviewed
)]
setnames(patient_counts, "N", "cells")
patient_counts[, cluster_numeric := as.integer(working_cluster)]
setorder(patient_counts, geo_subject, geo_accession, cluster_numeric)
patient_counts[, cluster_numeric := NULL]
write_tsv(patient_counts, "gse205506_patient_sample_cluster_cell_counts.tsv")

is_nuisance <- function(gene) {
  grepl("^(MT-|RPL|RPS|MTRNR)", gene) | gene %in% c("MALAT1")
}
author_top50 <- fread(
  file.path(output_dir, "gse205506_author_marker_reference_top50.tsv"),
  sep = "\t", header = TRUE, na.strings = NULL
)
our_sets <- top50[!is_nuisance(gene), .(our_genes = list(unique(gene))), by = cluster]
author_sets <- author_top50[!is_nuisance(gene), .(author_genes = list(unique(gene))), by = .(
  compartment_sheet,
  author_cluster = cluster
)]
clean_overlap <- rbindlist(lapply(seq_len(nrow(our_sets)), function(i) {
  current_cluster <- our_sets$cluster[[i]]
  current_genes <- our_sets$our_genes[[i]]
  rbindlist(lapply(seq_len(nrow(author_sets)), function(j) {
    author_genes <- author_sets$author_genes[[j]]
    shared <- intersect(current_genes, author_genes)
    data.table(
      cluster = current_cluster,
      compartment_sheet = author_sets$compartment_sheet[[j]],
      author_cluster = author_sets$author_cluster[[j]],
      our_non_nuisance_top50_genes = length(current_genes),
      author_non_nuisance_top50_genes = length(author_genes),
      overlap_count = length(shared),
      overlap_genes = paste(shared, collapse = ";"),
      jaccard = length(shared) / length(union(current_genes, author_genes))
    )
  }))
}))
clean_overlap <- merge(
  clean_overlap,
  review[, .(cluster, detailed_annotation, broad_annotation, expected_author_sheet)],
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)
clean_overlap[, expected_sheet_match := expected_author_sheet != "" & compartment_sheet == expected_author_sheet]
clean_overlap[, cluster_numeric := as.integer(cluster)]
setorder(clean_overlap, cluster_numeric, -expected_sheet_match, -overlap_count, -jaccard, compartment_sheet, author_cluster)
clean_overlap[, overlap_rank_all := seq_len(.N), by = cluster]
clean_overlap[, overlap_rank_expected_sheet := if (expected_sheet_match[[1]]) seq_len(.N) else NA_integer_, by = .(cluster, expected_sheet_match)]
clean_overlap[, cluster_numeric := NULL]
write_tsv(clean_overlap, "gse205506_author_marker_overlap_top50_non_nuisance_all.tsv")
write_tsv(clean_overlap[overlap_rank_all <= 3L], "gse205506_author_marker_overlap_top50_non_nuisance_top3.tsv")
write_tsv(
  clean_overlap[expected_sheet_match == TRUE & overlap_rank_expected_sheet <= 3L],
  "gse205506_author_marker_overlap_expected_compartment_top3.tsv"
)

broad_levels <- c(
  "Epithelial", "T/NK", "B cell", "Plasma", "Myeloid",
  "Endothelial", "Fibroblast/SMC", "Erythroid", "Immune unresolved", "Unresolved"
)
object$cell_type_broad_reviewed <- factor(object$cell_type_broad_reviewed, levels = broad_levels)
object$working_cluster <- factor(as.character(object$working_cluster), levels = as.character(0:27))
Idents(object) <- "working_cluster"

broad_palette <- c(
  "Epithelial" = "#9E2F00",
  "T/NK" = "#006B45",
  "B cell" = "#005A9C",
  "Plasma" = "#B56A00",
  "Myeloid" = "#9E3F75",
  "Endothelial" = "#0072A8",
  "Fibroblast/SMC" = "#5A3A00",
  "Erythroid" = "#8B0015",
  "Immune unresolved" = "#4A2E88",
  "Unresolved" = "#333333"
)
quality_palette <- c(
  "retained_after_qc" = "#333333",
  "mitochondrial_marker_enriched" = "#9B0000"
)
cluster_palette <- setNames(
  grDevices::hcl(h = seq(15, 375, length.out = 29L)[1:28], c = 85, l = 45),
  as.character(0:27)
)
timepoint_palette <- c("pre-treatment" = "#005A9C", "post-treatment" = "#9E2F00")
response_palette <- c("pCR" = "#006B45", "non-pCR" = "#9E3F75")

p_cluster <- DimPlot(
  object,
  reduction = "umap",
  group.by = "working_cluster",
  label = TRUE,
  repel = TRUE,
  label.size = 3.5,
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 0.65,
  cols = cluster_palette
) +
  labs(title = "GSE205506 working clusters (resolution 0.6)", color = "Cluster") +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 4))) +
  theme_classic(base_size = 11) +
  theme(aspect.ratio = 1)
save_figure(p_cluster, "gse205506_umap_working_clusters_reviewed", 9.5, 8)

p_broad <- DimPlot(
  object,
  reduction = "umap",
  group.by = "cell_type_broad_reviewed",
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  label.size = 3.6,
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 0.65,
  cols = broad_palette
) +
  labs(title = "GSE205506 reviewed broad annotation", color = "Cell type") +
  theme_classic(base_size = 11) +
  theme(aspect.ratio = 1)
save_figure(p_broad, "gse205506_umap_broad_cell_type_reviewed", 9.5, 8)

p_quality <- DimPlot(
  object,
  reduction = "umap",
  group.by = "cluster_quality_reviewed",
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 0.65,
  cols = quality_palette
) +
  labs(title = "Clusters with mitochondrial marker enrichment", color = "Cluster review") +
  theme_classic(base_size = 11) +
  theme(aspect.ratio = 1)
save_figure(p_quality, "gse205506_umap_cluster_quality_reviewed", 9.5, 8)

p_timepoint <- DimPlot(
  object,
  reduction = "umap",
  group.by = "derived_timepoint_from_exact_geo_treatment",
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 0.65,
  cols = timepoint_palette
) +
  labs(title = "Treatment timepoint", color = "Timepoint") +
  theme_classic(base_size = 11) +
  theme(aspect.ratio = 1)
save_figure(p_timepoint, "gse205506_umap_timepoint_reviewed", 9.5, 8)

p_response <- DimPlot(
  object,
  reduction = "umap",
  group.by = "table_s1_response",
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 0.65,
  cols = response_palette
) +
  labs(title = "Table S1 pathological response", color = "Response") +
  theme_classic(base_size = 11) +
  theme(aspect.ratio = 1)
save_figure(p_response, "gse205506_umap_response_reviewed", 9.5, 8)

dot_features <- c(
  "EPCAM", "KRT8", "KRT18", "KRT20", "CEACAM5", "MUC2", "BEST4", "OTOP2", "POU2F3", "TRPM5",
  "CD3D", "CD3E", "TRAC", "IL7R", "NKG7", "GNLY",
  "MS4A1", "CD79A", "CD74", "JCHAIN", "MZB1", "XBP1",
  "LST1", "TYROBP", "FCER1G", "C1QA", "S100A8",
  "PECAM1", "VWF", "CLDN5", "CCL21", "LYVE1",
  "COL1A1", "COL1A2", "DCN", "ACTA2", "TAGLN",
  "HBB", "HBA1", "ALAS2", "MKI67", "TOP2A", "STMN1"
)
dot_features <- intersect(dot_features, rownames(object))
p_dot <- DotPlot(
  object,
  features = dot_features,
  group.by = "working_cluster",
  assay = "RNA",
  dot.scale = 5
) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  labs(
    x = NULL,
    y = "Working cluster",
    title = "Cluster-level canonical marker evidence",
    color = "Scaled mean",
    size = "Percent expressed"
  ) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8))
save_figure(p_dot, "gse205506_marker_dotplot_working_clusters_reviewed", 17, 11)

object@misc$reviewed_annotation <- list(
  source = "Cluster Top10/Top50 presto markers, canonical marker dotplot, cluster QC, and author mmc3 marker reference",
  quality_rule = "top10_mitochondrial_markers >= 5 marks mitochondrial_marker_enriched; cells are retained",
  author_overlap_rule = "Top50 overlap after excluding MT-, RPL, RPS, MTRNR, and MALAT1; epithelial clusters have no matching mmc3 sheet"
)
saveRDS(object, reviewed_rds, compress = FALSE)

reviewed_summary <- data.table(
  object_path = normalizePath(reviewed_rds, winslash = "/", mustWork = TRUE),
  object_size_bytes = file.info(reviewed_rds)$size,
  cells = ncol(object),
  clusters = uniqueN(object$working_cluster),
  detailed_annotations = uniqueN(object$cell_type_detailed_reviewed),
  broad_annotations = uniqueN(object$cell_type_broad_reviewed),
  mitochondrial_marker_enriched_cells = sum(object$cluster_quality_reviewed == "mitochondrial_marker_enriched"),
  mitochondrial_marker_enriched_fraction = mean(object$cluster_quality_reviewed == "mitochondrial_marker_enriched")
)
write_tsv(reviewed_summary, "gse205506_reviewed_object_summary.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_review_annotation_and_author_markers.txt"))
message("Reviewed annotation and author marker comparison completed: ", reviewed_rds)
