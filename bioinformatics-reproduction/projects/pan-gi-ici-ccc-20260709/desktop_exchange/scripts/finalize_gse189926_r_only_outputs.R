suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(340)

project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
exchange_dir <- file.path(project_dir, "desktop_exchange")
out_dir <- file.path(exchange_dir, "uploads/20260710_gse189926_r_only_rerun")
object_dir <- file.path(project_dir, "04_objects/20260710_gse189926_r_only_rerun")
source_summary_dir <- file.path(exchange_dir, "uploads/20260709_response_object_construction")
prelabel_path <- file.path(object_dir, "gse189926_r_qc_pass_prelabel.rds")
unfiltered_plot_path <- file.path(object_dir, "gse189926_r_unfiltered_umap_plot_data.rds")
final_object_path <- file.path(object_dir, "gse189926_r_qc_pass_annotated.rds")
raw_object_path <- file.path(object_dir, "gse189926_r_unfiltered_raw.rds")

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "")
}

log_step <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
  flush.console()
}

if (!file.exists(prelabel_path)) stop(sprintf("Missing stage-1 object: %s", prelabel_path))
if (!file.exists(unfiltered_plot_path)) stop(sprintf("Missing unfiltered UMAP data: %s", unfiltered_plot_path))

label_audit <- data.table(
  cluster = as.character(0:18),
  resolution = 0.5,
  broad_label = c(
    "Granulocyte_like", "T_cells", "CD8_T_NK", "Plasma_cells", "Inflammatory_myeloid",
    "Epithelial_like", "Macrophages", "Monocytes", "Plasma_cells", "Mast_cells",
    "Stromal_perivascular_like", "Epithelial_like", "Endothelial_like", "Dendritic_cells",
    "Unresolved", "Epithelial_like", "Unresolved", "Plasma_cells", "Plasma_cells"
  ),
  refined_label = c(
    "Unresolved_low_signal", "CD4_T_naive_memory", "CD8_T_cytotoxic", "Plasma", "Monocyte_classical",
    "Unresolved_low_signal", "Macrophage_C1QC_APOE", "Monocyte_classical", "Plasma", "Mast",
    "Unresolved_low_signal", "Unresolved_low_signal", "Unresolved_low_signal", "DC_LAMP3",
    "Unresolved_low_signal", "Unresolved_low_signal", "Unresolved_low_signal", "Plasma", "Plasma"
  ),
  supporting_marker_genes = c(
    "CXCL8,S100A9,S100A8,FCGR3B,CSF3R,HCAR2",
    "TRBC2,TRAC,CD3D,CD3E,IL7R,LTB",
    "CCL5,NKG7,GZMB,GNLY,GZMA,TRBC2,TRAC,CD3D,CD8A,PRF1",
    "JCHAIN,MZB1,XBP1,DERL3,FKBP11,SDC1",
    "CXCL8,S100A9,S100A8,IL1B,FCER1G,IL1RN",
    "AGR2,KRT8,KRT18,KRT19,EPCAM,PGC,TFF2,PIGR",
    "C1QA,C1QB,C1QC,APOE,LYZ,FCER1G,CD74",
    "IL1B,LYZ,S100A9,S100A8,VCAN,FCER1G,LST1",
    "JCHAIN,MZB1,XBP1,DERL3,FKBP11",
    "TPSB2,TPSAB1,CPA3,MS4A2,KIT,GATA2",
    "COL3A1,COL1A1,COL1A2,RGS5,DCN,SPARC",
    "KRT13,KRT4,KRT5,KRT6A,TACSTD2,CLDN4",
    "PLVAP,PECAM1,VWF,AQP1,GNG11,SPARCL1",
    "IDO1,LAMP3,FSCN1,CD74,HLA-DRA,HLA-DPA1,GBP1",
    "MIR663AHG,LINC00632,RN7SK,NEAT1,JUN,HSPA1B",
    "EPCAM,KRT8,KRT18,CLDN4,GPX2,SERPINA1,CD24",
    "PGA4,PGA3,LIPF,CXCR4,TRBC1,BTG1",
    "JCHAIN,MZB1,DERL3,HSP90B1,IGLL5",
    "JCHAIN,MZB1,IGLL5,CD79A,XBP1"
  ),
  conflicting_marker_genes = c(
    "FCGR3B,CSF3R indicate granulocyte-like identity outside the requested refined label space",
    "CD74,JCHAIN ambient or mixed signal",
    "GNLY,NKG7 indicate an NK-like component but TRBC2,TRAC,CD3D support cytotoxic T identity",
    "S100A9,LYZ ambient myeloid signal",
    "FCGR3B indicates a granulocyte-like component",
    "LYZ and immune ambient signal",
    "S100A9 indicates an inflammatory myeloid component",
    "CD74 and HLA class II indicate activated antigen-presentation state",
    "S100A9,LYZ ambient myeloid signal",
    "STMN1 is present but canonical cycling genes do not dominate",
    "JCHAIN and CD74 ambient or mixed signal",
    "S100A9 and CD74 mixed inflammatory signal",
    "CD74 and JCHAIN mixed immune signal",
    "S100A9,S100A8 indicate an inflammatory component",
    "No stable canonical immune lineage program",
    "LYZ and SPP1 mixed signal",
    "T-cell genes and gastric epithelial genes coexist",
    "CD74 and myeloid ambient signal",
    "S100A9,S100A8 and LYZ mixed signal"
  ),
  marker_summary = c(
    "Strong granulocyte-like inflammatory program; not forced into a monocyte subtype.",
    "TCR genes with IL7R/naive-memory signal.",
    "Cytotoxic program with TCR genes; retained as cytotoxic CD8 T.",
    "Plasma-cell secretory program.",
    "Inflammatory classical-monocyte-like program.",
    "Gastric epithelial program outside the immune label space.",
    "C1Q/APOE macrophage program.",
    "Inflammatory monocyte program with VCAN and S100 genes.",
    "Plasma-cell secretory program.",
    "Canonical mast-cell program.",
    "Fibroblast/perivascular extracellular-matrix program outside the immune label space.",
    "Squamous-like epithelial keratin program outside the immune label space.",
    "Endothelial vascular program outside the immune label space.",
    "Mature LAMP3/FSCN1 dendritic-cell program.",
    "Stress and long-noncoding-RNA signal without a stable lineage program.",
    "Epithelial program outside the immune label space.",
    "Small mixed gastric epithelial/T-cell group without stable lineage separation.",
    "Small plasma-cell group.",
    "Small plasma-cell group with mixed inflammatory ambient signal."
  ),
  label_confidence_level = c(
    "high_broad_low_refined", "moderate", "high", "high", "moderate",
    "high_broad_low_refined", "high", "moderate", "high", "high",
    "high_broad_low_refined", "high_broad_low_refined", "high_broad_low_refined", "high",
    "low", "high_broad_low_refined", "low", "high", "moderate"
  ),
  label_decision_note = c(
    "Broad granulocyte-like identity retained; requested refined space has no granulocyte label, so refined label is Unresolved_low_signal.",
    "Assigned from TCR and IL7R evidence.",
    "Assigned from cytotoxic genes plus TCR evidence.",
    "Assigned from JCHAIN/MZB1/XBP1 secretory program.",
    "Assigned from inflammatory S100/IL1B myeloid program.",
    "Not forced into an immune subtype.",
    "Assigned from C1QA/C1QB/C1QC/APOE evidence.",
    "Assigned from LYZ/S100/VCAN inflammatory monocyte evidence.",
    "Assigned from JCHAIN/MZB1/XBP1 secretory program.",
    "Assigned from TPSAB1/TPSB2/CPA3/KIT evidence.",
    "Not forced into an immune subtype.",
    "Not forced into an immune subtype.",
    "Not forced into an immune subtype.",
    "Assigned from LAMP3/FSCN1/IDO1 antigen-presentation program.",
    "Retained unresolved because lineage evidence is weak.",
    "Not forced into an immune subtype.",
    "Retained unresolved because epithelial and T-cell evidence conflict.",
    "Assigned from plasma-cell secretory genes despite small cell count.",
    "Assigned from JCHAIN/MZB1/IGLL5 evidence with mixed inflammatory signal recorded."
  )
)

log_step("Loading QC-pass Seurat object")
object <- readRDS(prelabel_path)
cluster_values <- as.character(object[["cluster_harmony_r0_5", drop = TRUE]])
if (!identical(sort(unique(cluster_values)), sort(label_audit$cluster))) {
  stop("The exact Harmony resolution 0.5 cluster values do not match the reviewed label audit")
}
label_index <- match(cluster_values, label_audit$cluster)
object$broad_label <- label_audit$broad_label[label_index]
object$refined_label <- label_audit$refined_label[label_index]
object$annotation_resolution <- "Harmony_0.5"
object$cluster_label_display <- paste0(cluster_values, " | ", object$refined_label)
Idents(object) <- "cluster_harmony_r0_5"
write_tsv(label_audit, file.path(out_dir, "gse189926_r_refined_label_audit.tsv"))

resolution_summary_path <- file.path(out_dir, "gse189926_r_cluster_resolution_summary.tsv")
resolution_summary <- fread(resolution_summary_path)
resolution_summary[, selected_for_labeling := resolution == 0.5]
resolution_summary[, selection_reason := fifelse(
  resolution == 0.5,
  "Selected after marker review: stable biological groups without unsupported fine fragmentation.",
  fifelse(
    resolution == 0.3,
    "Not selected because major immune and nonimmune programs were less separated.",
    "Not selected because 25 groups added small mixed groups without additional supported immune label distinctions."
  )
)]
write_tsv(resolution_summary, resolution_summary_path)

metadata <- as.data.table(object[[]], keep.rownames = "cell_id")
metadata[, raw_outcome := get("characteristics_ch1::outcome")]
label_by_sample <- metadata[, .(
  cell_count = .N
), by = .(
  cluster = cluster_harmony_r0_5,
  resolution = annotation_resolution,
  broad_label,
  refined_label,
  sample_accession,
  patient_id,
  raw_outcome,
  binary_response_group,
  timepoint,
  timepoint_simple
)]
label_by_sample[, proportion_within_cluster_sample := cell_count / sum(cell_count), by = .(cluster, sample_accession)]
write_tsv(label_by_sample, file.path(out_dir, "gse189926_r_label_by_cluster_sample.tsv"))

celltype_counts <- metadata[, .(cell_count = .N), by = .(
  raw_outcome,
  binary_response_group,
  timepoint,
  timepoint_simple,
  broad_label,
  refined_label
)]
celltype_counts[, proportion_within_outcome_timepoint := cell_count / sum(cell_count), by = .(
  raw_outcome, binary_response_group, timepoint, timepoint_simple
)]
celltype_counts[, table_scope := "QC-pass cells only"]
celltype_counts[, metadata_note := "PR -> responder; SD/PD -> non_responder; NE and all other raw labels -> not_in_binary_response. Raw outcome is preserved."]
setorder(celltype_counts, raw_outcome, timepoint, broad_label, refined_label)
write_tsv(celltype_counts, file.path(out_dir, "gse189926_r_celltype_counts_by_response_timepoint.tsv"))

palette_for_levels <- function(levels, palette_name = "Dynamic") {
  colors <- grDevices::hcl.colors(max(3L, length(levels)), palette = palette_name)[seq_along(levels)]
  setNames(colors, levels)
}

unfiltered_data <- as.data.table(readRDS(unfiltered_plot_path))
qc_metadata_for_join <- metadata[, .(cell_id, broad_label, refined_label)]
unfiltered_data <- merge(unfiltered_data, qc_metadata_for_join, by = "cell_id", all.x = TRUE, sort = FALSE)
unfiltered_data[is.na(broad_label), broad_label := "QC_failed_not_labeled"]
unfiltered_data[is.na(refined_label), refined_label := "QC_failed_not_labeled"]
unfiltered_data[, raw_outcome := get("characteristics_ch1::outcome")]

qc_unintegrated_coords <- Embeddings(object, reduction = "umap.qc")
harmony_coords <- Embeddings(object, reduction = "umap.harmony")
qc_plot_data <- copy(metadata)
qc_plot_data[, `:=`(
  UMAP_1 = qc_unintegrated_coords[cell_id, 1L],
  UMAP_2 = qc_unintegrated_coords[cell_id, 2L]
)]
harmony_plot_data <- copy(metadata)
harmony_plot_data[, `:=`(
  UMAP_1 = harmony_coords[cell_id, 1L],
  UMAP_2 = harmony_coords[cell_id, 2L]
)]

all_sample_levels <- sort(unique(unfiltered_data$sample_accession))
all_patient_levels <- sort(unique(unfiltered_data$patient_id))
sample_palette <- palette_for_levels(all_sample_levels, "Dynamic")
patient_palette <- palette_for_levels(all_patient_levels, "Dark 3")
outcome_palette <- c(PD = "#D55E00", SD = "#0072B2", PR = "#009E73", NE = "#7A5195")
timepoint_palette <- c(
  "pre-treatment" = "#0072B2",
  "post-treatment" = "#D55E00",
  "post-treatment1" = "#D55E00",
  "post-treatment2" = "#009E73"
)
broad_palette <- c(
  Granulocyte_like = "#E69F00",
  T_cells = "#0072B2",
  CD8_T_NK = "#56B4E9",
  Plasma_cells = "#CC79A7",
  Inflammatory_myeloid = "#D55E00",
  Epithelial_like = "#B79F00",
  Macrophages = "#009E73",
  Monocytes = "#006D2C",
  Mast_cells = "#7A5195",
  Stromal_perivascular_like = "#8C564B",
  Endothelial_like = "#17BECF",
  Dendritic_cells = "#E15759",
  Unresolved = "#666666",
  QC_failed_not_labeled = "#BDBDBD"
)
refined_palette <- c(
  CD4_T_naive_memory = "#0072B2",
  CD8_T_cytotoxic = "#56B4E9",
  Plasma = "#CC79A7",
  Monocyte_classical = "#D55E00",
  Macrophage_C1QC_APOE = "#009E73",
  Mast = "#7A5195",
  DC_LAMP3 = "#E15759",
  Unresolved_low_signal = "#666666",
  QC_failed_not_labeled = "#BDBDBD"
)

plot_embedding_panel <- function(data, field, field_title, palette) {
  plot_data <- copy(data)
  plot_data[, plot_value := as.character(get(field))]
  set.seed(340)
  plot_data <- plot_data[sample.int(.N)]
  ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = plot_value)) +
    geom_point(size = 0.11, alpha = 0.90, stroke = 0) +
    scale_color_manual(values = palette, na.value = "#BDBDBD", drop = FALSE) +
    coord_equal() +
    labs(title = field_title, color = NULL) +
    guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1), ncol = 1)) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      legend.position = "right",
      legend.text = element_text(size = 6.5),
      legend.key.height = grid::unit(3.2, "mm"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

plot_embedding_grid <- function(data, path, title) {
  panels <- list(
    plot_embedding_panel(data, "sample_accession", "Sample accession", sample_palette),
    plot_embedding_panel(data, "patient_id", "Patient", patient_palette),
    plot_embedding_panel(data, "raw_outcome", "Raw outcome", outcome_palette),
    plot_embedding_panel(data, "timepoint", "Timepoint", timepoint_palette),
    plot_embedding_panel(data, "broad_label", "Broad label", broad_palette),
    plot_embedding_panel(data, "refined_label", "Refined label", refined_palette)
  )
  combined <- wrap_plots(panels, ncol = 3) + plot_annotation(
    title = title,
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )
  ggsave(path, combined, width = 24, height = 15, dpi = 240, bg = "white", limitsize = FALSE)
}

log_step("Generating three R UMAP comparison figures")
plot_embedding_grid(
  unfiltered_data,
  file.path(out_dir, "gse189926_r_umap_unintegrated_qc.png"),
  "GSE189926 R/Seurat unfiltered unintegrated UMAP"
)
plot_embedding_grid(
  qc_plot_data,
  file.path(out_dir, "gse189926_r_umap_qc_pass.png"),
  "GSE189926 R/Seurat QC-pass unintegrated UMAP"
)
plot_embedding_grid(
  harmony_plot_data,
  file.path(out_dir, "gse189926_r_umap_harmony.png"),
  "GSE189926 R/Seurat QC-pass Harmony UMAP"
)

marker_groups <- list(
  T_NK_ILC = c("CD3D", "CD3E", "TRAC", "CD4", "IL7R", "CCR7", "FOXP3", "IL2RA", "CTLA4", "CD8A", "CD8B", "GZMK", "CCL5", "NKG7", "GNLY", "GZMB", "PRF1", "KLRD1", "TRDC"),
  B_Plasma = c("MS4A1", "CD79A", "CD74", "CD37", "CD22", "CD27", "CD38", "MZB1", "JCHAIN", "XBP1", "SDC1"),
  Monocyte_Macrophage_DC = c("LYZ", "LST1", "S100A8", "S100A9", "FCGR3A", "VCAN", "C1QA", "C1QB", "C1QC", "APOE", "SPP1", "CD1C", "CLEC10A", "FCER1A", "CLEC9A", "XCR1", "IL3RA", "LAMP3", "FSCN1", "CCL19"),
  Mast_Cycling_Nonimmune_QC = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "MKI67", "TOP2A", "STMN1", "EPCAM", "KRT8", "KRT18", "KRT19", "COL1A1", "COL1A2", "COL3A1", "RGS5", "PECAM1", "VWF", "PLVAP")
)
marker_groups <- lapply(marker_groups, function(features) features[features %in% rownames(object)])
cluster_order <- label_audit[order(as.integer(cluster)), paste0(cluster, " | ", refined_label)]
object$cluster_label_display <- factor(object$cluster_label_display, levels = cluster_order)
dot_panels <- lapply(names(marker_groups), function(group_name) {
  DotPlot(object, features = marker_groups[[group_name]], group.by = "cluster_label_display", assay = "RNA", dot.scale = 4) +
    scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
    labs(title = gsub("_", " / ", group_name), x = NULL, y = NULL, color = "Scaled mean", size = "% expressing") +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 11, face = "bold"),
      panel.grid.major = element_line(linewidth = 0.2, color = "#E5E5E5"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.margin = margin(5, 8, 5, 5)
    )
})
marker_dotplot <- wrap_plots(dot_panels, ncol = 2) + plot_annotation(
  title = "GSE189926 R/Seurat marker review by Harmony cluster (resolution 0.5)",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
)
ggsave(
  file.path(out_dir, "gse189926_r_marker_dotplot.png"),
  marker_dotplot,
  width = 24,
  height = 22,
  dpi = 240,
  bg = "white",
  limitsize = FALSE
)

log_step("Building GSE235863 R patient-level summaries")
cd45_major <- fread(file.path(source_summary_dir, "gse235863_cd45_response_celltype_counts.tsv"))
cd45_sub <- fread(file.path(source_summary_dir, "gse235863_cd45_response_subcluster_counts.tsv"))
cd8_input <- fread(file.path(source_summary_dir, "gse235863_cd8_doublet_sensitivity_counts.tsv"))

cd45_patient <- cd45_major[, .(cell_count = sum(cell_count)), by = .(
  patient, raw_response, derived_response_group, response_ambiguity_flag, tissue, major_cluster
)]
cd45_patient[, proportion_within_patient_tissue := cell_count / sum(cell_count), by = .(patient, tissue)]
setorder(cd45_patient, patient, tissue, major_cluster)
write_tsv(cd45_patient, file.path(out_dir, "gse235863_cd45_r_patient_level_composition.tsv"))

cd45_sub_patient <- cd45_sub[, .(cell_count = sum(cell_count)), by = .(
  patient, raw_response, derived_response_group, response_ambiguity_flag, tissue, sub_cluster
)]
cd45_sub_patient[, proportion_within_patient_tissue := cell_count / sum(cell_count), by = .(patient, tissue)]
setorder(cd45_sub_patient, patient, tissue, sub_cluster)
write_tsv(cd45_sub_patient, file.path(out_dir, "gse235863_cd45_r_subcluster_patient_proportions.tsv"))

cd8_patient <- cd8_input[, .(cell_count = sum(cell_count)), by = .(
  patient, raw_response, derived_response_group, tissue, doublet_filter_state, sub_cluster
)]
cd8_patient[, proportion_within_patient_tissue := cell_count / sum(cell_count), by = .(patient, tissue, doublet_filter_state)]
setorder(cd8_patient, patient, tissue, doublet_filter_state, sub_cluster)
write_tsv(cd8_patient, file.path(out_dir, "gse235863_cd8_r_patient_level_doublet_sensitivity.tsv"))

major_palette <- palette_for_levels(sort(unique(cd45_patient$major_cluster)), "Dark 3")
cd45_plot_data <- copy(cd45_patient)
cd45_plot_data[, patient_tissue := paste(patient, tissue, sep = " | ")]
cd45_plot_data[, ambiguity_label := fifelse(as.logical(response_ambiguity_flag), "P18 response ambiguous", "Response metadata matched")]
cd45_plot_data[, ambiguity_short := fifelse(as.logical(response_ambiguity_flag), "P18 ambiguous", "Matched")]
cd45_plot_data[, response_short := fifelse(derived_response_group == "non_responder", "NR", "R")]
cd45_plot_data[, facet_label := paste(ambiguity_short, response_short, sep = " | ")]
cd45_plot <- ggplot(cd45_plot_data, aes(x = patient_tissue, y = proportion_within_patient_tissue, fill = major_cluster)) +
  geom_col(width = 0.82) +
  facet_grid(. ~ facet_label, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = major_palette) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "GSE235863 CD45 patient-level major-cluster composition (R summary)",
    x = "Patient | tissue", y = "Proportion", fill = "Major cluster",
    caption = "P18 ambiguity preserved: source response metadata contains CR/responder and PR/responder."
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )
ggsave(file.path(out_dir, "gse235863_cd45_r_response_composition_qc.png"), cd45_plot, width = 16, height = 8, dpi = 240, bg = "white")

cd8_palette <- palette_for_levels(sort(unique(cd8_patient$sub_cluster)), "Dynamic")
cd8_plot_data <- copy(cd8_patient)
cd8_plot_data[, patient_tissue := paste(patient, tissue, sep = " | ")]
cd8_plot_data[, filter_label := fifelse(
  doublet_filter_state == "before_doublet_exclusion",
  "Before doublet exclusion",
  "After predicted-doublet exclusion"
)]
cd8_plot_data[, facet_label := paste(filter_label, derived_response_group, sep = " | ")]
cd8_plot <- ggplot(cd8_plot_data, aes(x = patient_tissue, y = proportion_within_patient_tissue, fill = sub_cluster)) +
  geom_col(width = 0.82) +
  facet_wrap(~ facet_label, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = cd8_palette) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "GSE235863 CD8 patient-level composition before and after predicted-doublet exclusion (R summary)",
    x = "Patient | tissue", y = "Proportion", fill = "CD8 subcluster"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 7),
    plot.title = element_text(face = "bold")
  )
ggsave(file.path(out_dir, "gse235863_cd8_r_doublet_excluded_composition_qc.png"), cd8_plot, width = 16, height = 10, dpi = 240, bg = "white")

log_step(sprintf("Saving final annotated R object: %s", final_object_path))
saveRDS(object, final_object_path, compress = FALSE)

if (!requireNamespace("digest", quietly = TRUE)) stop("R package digest is required for SHA256 object inventory")
object_paths <- c(raw_object_path, prelabel_path, final_object_path, unfiltered_plot_path)
object_names <- c(
  "gse189926_r_unfiltered_raw",
  "gse189926_r_qc_pass_prelabel",
  "gse189926_r_qc_pass_annotated",
  "gse189926_r_unfiltered_umap_plot_data"
)
object_formats <- c("rds_seurat", "rds_seurat", "rds_seurat", "rds_data_frame")
inventory_rows <- rbindlist(lapply(seq_along(object_paths), function(idx) {
  data.frame(
    dataset_id = "GSE189926",
    object_name = object_names[[idx]],
    pc_local_path = normalizePath(object_paths[[idx]], winslash = "\\", mustWork = TRUE),
    object_format = object_formats[[idx]],
    file_size_bytes = file.info(object_paths[[idx]])$size,
    sha256 = digest::digest(file = object_paths[[idx]], algo = "sha256", serialize = FALSE),
    upload_status = "pc_local_only_large_file",
    notes = c(
      "Unfiltered sparse raw-count Seurat object with QC metadata; all 121452 cells retained.",
      "QC-pass Seurat object with Seurat PCA/UMAP, Harmony reduction, three clustering resolutions, and normalized data.",
      "QC-pass annotated Seurat object with reviewed broad and refined labels.",
      "Compact plotting data for the unfiltered R UMAP; kept PC-local with the R object set."
    )[[idx]],
    stringsAsFactors = FALSE
  )
}))
write_tsv(inventory_rows, file.path(out_dir, "object_inventory.tsv"))

run_status <- data.table(
  dataset_id = c("GSE189926", "GSE235863_CD45", "GSE235863_CD8"),
  step = c(
    "R_sparse_import_qc_seurat_harmony_clustering_marker_annotation",
    "R_patient_level_composition_summary",
    "R_patient_level_doublet_sensitivity_summary"
  ),
  status = "completed",
  language = "R",
  package_route = c(
    "Seurat 5.4.0; SeuratObject 5.3.0; Matrix 1.7.4; data.table 1.18.2.1; harmony 1.2.4; ggplot2 4.0.2; patchwork 1.3.2",
    "data.table 1.18.2.1; ggplot2 4.0.2",
    "data.table 1.18.2.1; ggplot2 4.0.2"
  ),
  input_files = c(
    "22 raw gzip text matrices; gse189926_sample_qc.tsv; gse189926_sample_parse_audit.tsv",
    "gse235863_cd45_response_celltype_counts.tsv; gse235863_cd45_response_subcluster_counts.tsv",
    "gse235863_cd8_doublet_sensitivity_counts.tsv"
  ),
  output_files = c(
    "All required GSE189926 R audit tables, marker tables, 3 UMAP figures, marker dotplot, and local RDS objects",
    "gse235863_cd45_r_patient_level_composition.tsv; gse235863_cd45_r_subcluster_patient_proportions.tsv; gse235863_cd45_r_response_composition_qc.png",
    "gse235863_cd8_r_patient_level_doublet_sensitivity.tsv; gse235863_cd8_r_doublet_excluded_composition_qc.png"
  ),
  cell_count_input = c(121452L, sum(cd45_patient$cell_count), sum(cd8_input[doublet_filter_state == "before_doublet_exclusion"]$cell_count)),
  cell_count_used = c(ncol(object), sum(cd45_patient$cell_count), sum(cd8_patient$cell_count)),
  notes = c(
    "All GSE189926 analysis, QC, UMAP, clustering, marker export, and labeling were performed in R. No downsampling was used. Harmony grouped by sample_accession only.",
    "Existing uploaded count tables summarized in R; raw response and P18 ambiguity flag preserved.",
    "Existing uploaded doublet-sensitivity counts summarized in R; before and after predicted-doublet exclusion retained."
  ),
  error_message = ""
)
write_tsv(run_status, file.path(out_dir, "run_status.tsv"))

warning_log <- data.table(
  dataset_id = c("GSE189926", "GSE189926", "GSE189926", "GSE189926", "GSE235863_CD45"),
  warning_type = c(
    "refined_label_space_limit",
    "nonimmune_signal_in_cd45_input",
    "small_or_mixed_clusters",
    "resolution_selection",
    "response_metadata_ambiguity"
  ),
  affected_field = c("refined_label", "refined_label", "refined_label", "annotation_resolution", "response_ambiguity_flag"),
  affected_count = c(
    sum(cluster_values == "0"),
    sum(cluster_values %in% c("5", "10", "11", "12", "15")),
    sum(cluster_values %in% c("14", "16", "18")),
    2L,
    sum(cd45_patient[patient == "P18"]$cell_count)
  ),
  explanation = c(
    "Cluster 0 has strong FCGR3B/CSF3R granulocyte-like evidence, but the requested refined label space has no granulocyte label.",
    "Several clusters show epithelial, stromal/perivascular, or endothelial programs despite CD45+ source metadata.",
    "Clusters 14, 16, and 18 are small or contain conflicting marker programs.",
    "Resolutions 0.3 and 0.8 were audited but not selected for final annotation.",
    "P18 retains the source CR/responder versus PR/responder ambiguity flag."
  ),
  action_taken = c(
    "Broad Granulocyte_like retained; refined label set to Unresolved_low_signal without forcing a monocyte subtype.",
    "Broad nonimmune identity retained; refined labels left Unresolved_low_signal.",
    "Labels assigned only where marker evidence was stable; conflicts are recorded in the label audit.",
    "Harmony resolution 0.5 selected after marker review; all three resolution summaries remain available.",
    "Raw response and ambiguity flag preserved in R composition outputs and figure facet."
  )
)
write_tsv(warning_log, file.path(out_dir, "warning_log.tsv"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

status_lines <- c(
  "# GSE189926 R-only rerun status",
  "",
  "Status: complete",
  "Language: R",
  "",
  "- Raw input: 22 gzip text matrices",
  "- Cells imported: 121452",
  sprintf("- QC-pass cells: %d", ncol(object)),
  "- Feature union: 51757",
  "- Unintegrated routes: unfiltered and QC-pass Seurat log-normalized PCA/UMAP",
  "- Sample-aware route: Harmony by sample_accession",
  "- Clustering resolutions audited: 0.3, 0.5, 0.8",
  "- Annotation resolution: 0.5",
  "- Marker export: raw and technical-gene-excluded top 50 per cluster",
  "- Large RDS objects: PC-local only and recorded in object_inventory.tsv",
  "- GSE235863 CD45/CD8 summaries: rebuilt in R from existing uploaded count tables",
  "- Downsampling: none",
  "",
  "Stop conditions respected: CellChat, LIANA, NicheNet, differential communication, final statistical testing, and manuscript figures were not run."
)
writeLines(status_lines, file.path(out_dir, "STATUS.md"))
log_step("Final R-only outputs completed")
