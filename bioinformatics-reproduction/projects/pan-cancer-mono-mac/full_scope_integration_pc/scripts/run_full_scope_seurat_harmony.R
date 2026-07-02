#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(data.table)
  library(ggplot2)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- normalizePath(sub(file_arg, "", script_args[grepl(file_arg, script_args)][1]), mustWork = FALSE)

find_project_root <- function(start_dir) {
  current <- normalizePath(start_dir, mustWork = FALSE)
  for (i in 1:8) {
    marker <- file.path(current, "inputs", "full_scope_standard_object_manifest_template.tsv")
    if (file.exists(marker)) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  stop("Could not locate full_scope_integration_pc root from: ", start_dir)
}

start_dir <- if (nzchar(script_path)) dirname(script_path) else getwd()
project_root <- Sys.getenv("FULL_SCOPE_PROJECT_ROOT", unset = find_project_root(start_dir))

manifest_path <- Sys.getenv(
  "FULL_SCOPE_INPUT_MANIFEST",
  unset = file.path(project_root, "inputs", "full_scope_standard_object_manifest.tsv")
)
output_root <- Sys.getenv("FULL_SCOPE_OUTPUT_ROOT", unset = file.path(project_root, "outputs"))
out_dir <- file.path(output_root, "04_data_processed", "single_cell", "full_scope_integration")
fig_dir <- file.path(output_root, "06_figures", "main", "full_scope_integration")
tab_dir <- file.path(output_root, "07_tables", "main", "full_scope_integration")
log_path <- file.path(output_root, "99_logs", "full_scope_integration_log.tsv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

write_log <- function(step, status, note = "") {
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    step = step,
    status = status,
    note = note,
    stringsAsFactors = FALSE
  )
  write.table(row, log_path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = !file.exists(log_path), append = file.exists(log_path))
}

save_plot <- function(plot, stem, width, height) {
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300, limitsize = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, limitsize = FALSE)
  write_log("save_plot", "created", paste0(stem, "; png=", png_path, "; pdf=", pdf_path))
  invisible(c(png_path, pdf_path))
}

resolve_path <- function(path_value, root) {
  if (grepl("^/", path_value)) {
    return(path_value)
  }
  file.path(root, path_value)
}

set.seed(20260702)
threads <- as.integer(Sys.getenv("FULL_THREADS", unset = "4"))
npcs <- as.integer(Sys.getenv("FULL_NPCS", unset = "30"))
resolution <- as.numeric(Sys.getenv("FULL_RESOLUTION", unset = "1.0"))
batch_column <- Sys.getenv("FULL_BATCH_COLUMN", unset = "batch_id")
run_markers <- identical(Sys.getenv("FULL_RUN_MARKERS", unset = "TRUE"), "TRUE")
marker_max_cells <- as.integer(Sys.getenv("FULL_MARKER_MAX_CELLS_PER_CLUSTER", unset = "0"))
variable_feature_method <- Sys.getenv("FULL_VARIABLE_FEATURE_METHOD", unset = "vst")

expected_dataset_ids <- c(
  "set1_ESCA", "set1_kidney", "set1_LYM", "set1_MYE", "set1_OV_FTC", "set1_PAAD", "set1_THCA", "set1_UCEC",
  "set2_TNBC", "set2_LUAD", "set2_COAD_READ", "set2_LIHC_CHOL", "set2_healthy_lung", "set2_healthy_colon",
  "set2_healthy_liver", "set2_GBM", "set2_healthy_PBMC"
)

required_manifest_columns <- c("dataset_id", "paper_label", "rds_path", "use_integration", "data_state", "notes")
required_metadata_columns <- c(
  "paper_dataset_id", "paper_label", "patient_or_subject", "origin_or_tissue",
  "source_cell_id", "manual_annotation", "manual_broad_lineage", "batch_id"
)
allowed_lineages <- c("Macrophages", "Monocytes", "Classical DCs", "Plasmacytoid DCs", "Mast cells")

if (!file.exists(manifest_path)) {
  stop("Missing input manifest: ", manifest_path)
}

manifest <- fread(manifest_path)
missing_manifest_columns <- setdiff(required_manifest_columns, names(manifest))
if (length(missing_manifest_columns) > 0L) {
  stop("Input manifest missing columns: ", paste(missing_manifest_columns, collapse = ", "))
}

unexpected_use_values <- setdiff(unique(manifest$use_integration), c("TRUE", "FALSE"))
if (length(unexpected_use_values) > 0L) {
  stop("use_integration must contain only TRUE or FALSE. Observed: ", paste(unexpected_use_values, collapse = ", "))
}

enabled <- manifest[use_integration == "TRUE"]
missing_dataset_ids <- setdiff(expected_dataset_ids, enabled$dataset_id)
extra_dataset_ids <- setdiff(enabled$dataset_id, expected_dataset_ids)
if (length(extra_dataset_ids) > 0L) {
  stop("Unexpected dataset_id values in enabled manifest rows: ", paste(extra_dataset_ids, collapse = ", "))
}
if (length(missing_dataset_ids) > 0L) {
  stop("Full 17-dataset integration requires all expected dataset_id values enabled. Missing: ", paste(missing_dataset_ids, collapse = ", "))
}
enabled[, dataset_id := factor(dataset_id, levels = expected_dataset_ids)]
setorder(enabled, dataset_id)
enabled[, dataset_id := as.character(dataset_id)]

enabled[, resolved_rds_path := vapply(rds_path, resolve_path, character(1), root = project_root)]
missing_rds <- enabled[!file.exists(resolved_rds_path)]
if (nrow(missing_rds) > 0L) {
  stop("Missing Seurat RDS files: ", paste(paste0(missing_rds$dataset_id, "=", missing_rds$resolved_rds_path), collapse = "; "))
}

write_log("manifest_loaded", "complete", paste0("datasets=", nrow(enabled), "; manifest=", manifest_path))

load_and_validate_object <- function(row) {
  dataset_id <- row[["dataset_id"]]
  paper_label <- row[["paper_label"]]
  rds_path <- row[["resolved_rds_path"]]
  write_log("read_rds_start", "running", paste0(dataset_id, "; path=", rds_path))
  obj <- readRDS(rds_path)
  if (!inherits(obj, "Seurat")) {
    stop(dataset_id, " RDS does not contain a Seurat object.")
  }
  if (!"RNA" %in% names(obj@assays)) {
    stop(dataset_id, " Seurat object lacks RNA assay.")
  }
  DefaultAssay(obj) <- "RNA"
  missing_metadata <- setdiff(required_metadata_columns, names(obj@meta.data))
  if (length(missing_metadata) > 0L) {
    stop(dataset_id, " object metadata missing columns: ", paste(missing_metadata, collapse = ", "))
  }
  if (!all(as.character(obj$paper_dataset_id) == dataset_id)) {
    stop(dataset_id, " object metadata paper_dataset_id does not exactly match manifest dataset_id.")
  }
  if (!all(as.character(obj$paper_label) == paper_label)) {
    stop(dataset_id, " object metadata paper_label does not exactly match manifest paper_label.")
  }
  observed_lineages <- sort(unique(as.character(obj$manual_broad_lineage)))
  outside_lineages <- setdiff(observed_lineages, allowed_lineages)
  if (length(outside_lineages) > 0L) {
    stop(dataset_id, " object has manual_broad_lineage values outside the allowed set: ", paste(outside_lineages, collapse = ", "))
  }
  if (anyDuplicated(as.character(obj$source_cell_id)) > 0L) {
    stop(dataset_id, " object has duplicated source_cell_id values.")
  }
  new_names <- make.unique(paste(dataset_id, as.character(obj$source_cell_id), sep = "__"))
  obj <- RenameCells(obj, new.names = new_names)
  obj$full_scope_cell_id <- colnames(obj)
  write_log("read_rds_done", "complete", paste0(dataset_id, "; cells=", ncol(obj), "; genes=", nrow(obj)))
  obj
}

object_list <- lapply(seq_len(nrow(enabled)), function(i) load_and_validate_object(enabled[i]))
names(object_list) <- enabled$dataset_id

common_features <- Reduce(intersect, lapply(object_list, rownames))
if (length(common_features) < 3000L) {
  stop("Common RNA feature set has fewer than 3000 genes: ", length(common_features))
}
write_log("common_features", "complete", paste0("genes=", length(common_features)))

object_list <- lapply(names(object_list), function(dataset_id) {
  obj <- subset(object_list[[dataset_id]], features = common_features)
  obj <- FindVariableFeatures(obj, selection.method = variable_feature_method, nfeatures = 3000, verbose = TRUE)
  obj
})
names(object_list) <- enabled$dataset_id

integration_features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
if (length(integration_features) == 0L) {
  stop("No integration features selected.")
}
write_log("integration_features", "complete", paste0("genes=", length(integration_features)))

message("Merging ", length(object_list), " standard Seurat objects")
combined <- Reduce(function(x, y) merge(x = x, y = y, merge.data = TRUE), object_list)
DefaultAssay(combined) <- "RNA"
VariableFeatures(combined) <- integration_features

if (!batch_column %in% names(combined@meta.data)) {
  stop("Batch column not present in merged metadata: ", batch_column)
}
if (anyNA(combined@meta.data[[batch_column]])) {
  stop("Batch column contains NA values: ", batch_column)
}

write_log("merged", "complete", paste0("cells=", ncol(combined), "; genes=", nrow(combined)))

message("Scaling, PCA, Harmony, UMAP, and Louvain clustering")
combined <- ScaleData(combined, features = integration_features, verbose = TRUE)
combined <- RunPCA(combined, features = integration_features, npcs = npcs, verbose = TRUE)
combined <- RunHarmony(
  combined,
  group.by.vars = batch_column,
  reduction.use = "pca",
  dims.use = seq_len(npcs),
  verbose = TRUE
)
combined <- RunUMAP(combined, reduction = "harmony", dims = seq_len(npcs), n.neighbors = 30, min.dist = 0.3, verbose = TRUE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = seq_len(npcs), verbose = TRUE)
combined <- FindClusters(combined, resolution = resolution, algorithm = 1, verbose = TRUE)

if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  combined <- JoinLayers(combined, assay = "RNA")
}

umap <- Embeddings(combined, "umap")
umap_dt <- data.table(
  full_scope_cell_id = colnames(combined),
  UMAP_1 = umap[, 1],
  UMAP_2 = umap[, 2],
  seurat_cluster = as.character(Idents(combined)),
  paper_dataset_id = as.character(combined$paper_dataset_id),
  paper_label = as.character(combined$paper_label),
  patient_or_subject = as.character(combined$patient_or_subject),
  origin_or_tissue = as.character(combined$origin_or_tissue),
  source_cell_id = as.character(combined$source_cell_id),
  manual_annotation = as.character(combined$manual_annotation),
  manual_broad_lineage = as.character(combined$manual_broad_lineage),
  batch_id = as.character(combined@meta.data[[batch_column]])
)
fwrite(umap_dt, file.path(tab_dir, "full_scope_umap_coordinates.tsv"), sep = "\t")

dataset_counts <- umap_dt[, .N, by = .(paper_dataset_id, paper_label)]
setnames(dataset_counts, "N", "cells")
fwrite(dataset_counts, file.path(tab_dir, "full_scope_cells_by_dataset.tsv"), sep = "\t")

origin_lineage_counts <- umap_dt[, .N, by = .(origin_or_tissue, manual_broad_lineage)]
setnames(origin_lineage_counts, "N", "cells")
origin_lineage_counts[, fraction := cells / sum(cells), by = origin_or_tissue]
fwrite(origin_lineage_counts, file.path(tab_dir, "full_scope_origin_lineage_counts.tsv"), sep = "\t")

lineage_colors <- c(
  "Macrophages" = "#E76F51",
  "Monocytes" = "#D9A21B",
  "Classical DCs" = "#C87F0A",
  "Plasmacytoid DCs" = "#8D5AA9",
  "Mast cells" = "#D8569C"
)

theme_full <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = grid::unit(3, "mm")
  )

plot_umap <- function(group_by, title, color_values = NULL, label = FALSE) {
  p <- DimPlot(combined, reduction = "umap", group.by = group_by, label = label, repel = TRUE, raster = TRUE) +
    ggtitle(title) +
    theme_full +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
  if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values, drop = FALSE)
  }
  p
}

save_plot(
  plot_umap("manual_broad_lineage", "Full-scope integrated myeloid cells by broad lineage", lineage_colors, label = TRUE),
  "full_scope_umap_by_manual_broad_lineage",
  width = 8,
  height = 6
)
save_plot(
  plot_umap("paper_label", "Full-scope integrated myeloid cells by dataset", NULL, label = FALSE),
  "full_scope_umap_by_dataset",
  width = 9,
  height = 6.4
)
save_plot(
  plot_umap("origin_or_tissue", "Full-scope integrated myeloid cells by origin", NULL, label = FALSE),
  "full_scope_umap_by_origin",
  width = 7.2,
  height = 5.8
)
save_plot(
  plot_umap("seurat_clusters", "Full-scope integrated myeloid cells by Seurat cluster", NULL, label = TRUE),
  "full_scope_umap_by_cluster",
  width = 8,
  height = 6
)

marker_genes <- c(
  "APOC1", "CD163", "C1QA", "LYVE1", "CTSB", "CTSS", "FCN1", "VCAN", "S100A8",
  "FTH1", "CD1C", "FCER1A", "CLEC10A", "LAMP3", "CLEC9A", "IL3RA", "GZMB", "CPA3", "KIT",
  "PTPRC", "LYZ", "LST1", "CST3", "CD68", "APOE", "C1QB", "TREM2", "SPP1", "THBS1",
  "FCGR3A", "MS4A7", "TPSAB1", "TPSB2", "MKI67", "TOP2A",
  "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "MS4A1", "CD79A", "EPCAM", "KRT19", "PECAM1", "VWF", "DCN", "COL1A1"
)
marker_genes <- intersect(marker_genes, rownames(combined))
if (length(marker_genes) > 0L) {
  p_dot <- DotPlot(combined, features = marker_genes, group.by = "manual_broad_lineage") +
    coord_flip() +
    ggtitle("Full-scope canonical marker audit by broad lineage") +
    theme_full +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  save_plot(p_dot, "full_scope_canonical_marker_dotplot", width = 8.5, height = 9.5)
}

if (run_markers) {
  message("Finding full-scope cluster markers")
  marker_args <- list(
    object = combined,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = TRUE
  )
  if (marker_max_cells > 0L) {
    marker_args$max.cells.per.ident <- marker_max_cells
  }
  markers <- do.call(FindAllMarkers, marker_args)
  fwrite(as.data.table(markers), file.path(tab_dir, "full_scope_cluster_markers_all.tsv"), sep = "\t")
  markers_dt <- as.data.table(markers)
  setorder(markers_dt, cluster, p_val_adj, -avg_log2FC)
  top10 <- markers_dt[, head(.SD, 10), by = cluster]
  fwrite(top10, file.path(tab_dir, "full_scope_cluster_top10_markers.tsv"), sep = "\t")
  write_log("markers", "complete", paste0("rows=", nrow(markers_dt), "; marker_max_cells_per_cluster=", marker_max_cells))
} else {
  write_log("markers", "skipped", "FULL_RUN_MARKERS was not TRUE")
}

run_summary <- data.table(
  metric = c(
    "paper_target_cells",
    "observed_integrated_cells",
    "difference_observed_minus_target",
    "datasets",
    "genes_after_common_feature_filter",
    "integration_features",
    "batch_column",
    "npcs",
    "resolution",
    "run_markers",
    "marker_max_cells_per_cluster"
  ),
  value = c(
    "131249",
    as.character(ncol(combined)),
    as.character(ncol(combined) - 131249L),
    as.character(uniqueN(combined$paper_dataset_id)),
    as.character(length(common_features)),
    as.character(length(integration_features)),
    batch_column,
    as.character(npcs),
    as.character(resolution),
    as.character(run_markers),
    as.character(marker_max_cells)
  )
)
fwrite(run_summary, file.path(tab_dir, "full_scope_integration_run_summary.tsv"), sep = "\t")

saveRDS(combined, file.path(out_dir, "full_scope_harmony_integrated_seurat.rds"))
write_log("script_done", "complete", paste0("cells=", ncol(combined), "; datasets=", uniqueN(combined$paper_dataset_id)))
message("Full-scope Harmony integration complete.")
