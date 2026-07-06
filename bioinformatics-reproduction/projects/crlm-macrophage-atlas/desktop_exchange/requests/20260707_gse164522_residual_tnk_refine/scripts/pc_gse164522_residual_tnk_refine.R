options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "Matrix", "Seurat", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
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

first_present <- function(values, universe, max_n = Inf) {
  out <- values[values %in% universe]
  unique(out)[seq_len(min(length(unique(out)), max_n))]
}

write_composition <- function(meta_dt, group_col, set_name, output_dir, output_key = group_col) {
  if (!group_col %in% names(meta_dt)) {
    return(NULL)
  }
  tab <- meta_dt[, .N, by = .(clean_cluster, value = get(group_col))]
  tab[, fraction := N / sum(N), by = clean_cluster]
  out_path <- file.path(output_dir, paste0("gse164522_", set_name, "_cluster_", output_key, "_composition.tsv"))
  fwrite(tab[order(clean_cluster, -fraction, value)], out_path, sep = "\t")
  invisible(tab)
}

script_path <- get_script_path()
request_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
output_dir <- file.path(request_dir, "outputs")
object_dir <- file.path(request_dir, "local_objects_not_for_github")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

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

cell_list_paths <- list(
  broad_refined_singlets = file.path(request_dir, "inputs", "gse164522_broad_refined_singlet_cells.tsv"),
  core_refined_singlets = file.path(request_dir, "inputs", "gse164522_core_refined_singlet_cells.tsv")
)
cleanup_labels_path <- file.path(request_dir, "inputs", "gse164522_cell_cleanup_labels.tsv.gz")
lineage_marker_path <- file.path(request_dir, "inputs", "lineage_marker_sets.tsv")
refinement_summary_path <- file.path(request_dir, "inputs", "gse164522_residual_tnk_refinement_input_summary.tsv")

missing_inputs <- c(
  unlist(cell_list_paths)[!file.exists(unlist(cell_list_paths))],
  cleanup_labels_path[!file.exists(cleanup_labels_path)],
  lineage_marker_path[!file.exists(lineage_marker_path)],
  refinement_summary_path[!file.exists(refinement_summary_path)]
)
if (length(missing_inputs) > 0) {
  stop("Missing input files: ", paste(missing_inputs, collapse = "; "))
}

message("Reading Seurat object: ", rds_path)
obj <- readRDS(rds_path)
old_cluster_col <- "RNA_snn_res.0.6"
if (!old_cluster_col %in% colnames(obj[[]])) {
  stop("Original cluster column not found in Seurat metadata: ", old_cluster_col)
}
obj[["original_cluster_res_0_6"]] <- as.character(obj[[old_cluster_col, drop = TRUE]])

cleanup_labels <- fread(cleanup_labels_path)
required_cleanup_cols <- c(
  "expression_cell_id", "cluster", "cleanup_label", "current_action",
  "review_label", "audit_status", "scDblFinder.class", "scDblFinder.score",
  "batch_sensitive_flag", "doublet_sensitive_cluster_flag"
)
missing_cleanup_cols <- setdiff(required_cleanup_cols, names(cleanup_labels))
if (length(missing_cleanup_cols) > 0) {
  stop("Cleanup label columns missing: ", paste(missing_cleanup_cols, collapse = ", "))
}
cleanup_labels[, expression_cell_id := as.character(expression_cell_id)]
cleanup_labels[, cluster := as.character(cluster)]

marker_sets <- fread(lineage_marker_path)
marker_sets[, gene := as.character(gene)]
refinement_summary <- fread(refinement_summary_path)
marker_features_requested <- c(
  "LYZ", "LST1", "FCN1", "S100A8", "S100A9", "C1QA", "C1QB", "C1QC",
  "APOE", "MRC1", "TREM2", "SPP1", "MS4A7", "CD3D", "TRAC", "NKG7",
  "CLEC10A", "FCER1A", "MS4A1", "MZB1", "TPSAB1", "PECAM1", "MKI67"
)
dotplot_features <- first_present(marker_features_requested, rownames(obj), max_n = 30)

run_one_set <- function(set_name, cell_path) {
  message("Running cleanup set: ", set_name)
  cells_dt <- fread(cell_path)
  required_cell_cols <- c("expression_cell_id", "cluster", "cleanup_label", "review_label", "current_action")
  missing_cell_cols <- setdiff(required_cell_cols, names(cells_dt))
  if (length(missing_cell_cols) > 0) {
    stop("Cell list columns missing for ", set_name, ": ", paste(missing_cell_cols, collapse = ", "))
  }
  cells_dt[, expression_cell_id := as.character(expression_cell_id)]
  requested_cells <- unique(cells_dt$expression_cell_id)
  missing_cells <- setdiff(requested_cells, colnames(obj))
  present_cells <- intersect(requested_cells, colnames(obj))
  if (length(missing_cells) > 0) {
    stop("Requested cells missing from Seurat object for ", set_name, ": ", length(missing_cells))
  }
  if (length(present_cells) < 100) {
    stop("Too few cells for ", set_name, ": ", length(present_cells))
  }

  input_audit <- data.table(
    set_name = set_name,
    requested_cells = length(requested_cells),
    present_cells = length(present_cells),
    missing_cells = length(missing_cells),
    input_clusters = uniqueN(cells_dt$cluster),
    input_cleanup_labels = paste(sort(unique(cells_dt$cleanup_label)), collapse = ";")
  )
  fwrite(input_audit, file.path(output_dir, paste0("gse164522_", set_name, "_input_cell_audit.tsv")), sep = "\t")

  clean_obj <- subset(obj, cells = present_cells)
  add_dt <- cleanup_labels[match(colnames(clean_obj), expression_cell_id)]
  if (!identical(add_dt$expression_cell_id, colnames(clean_obj))) {
    stop("Cleanup labels failed to align with Seurat columns for ", set_name)
  }
  add_meta <- as.data.frame(add_dt[, .(
    original_cluster_for_cleanup = cluster,
    cleanup_label,
    cleanup_action = current_action,
    cleanup_review_label = review_label,
    cleanup_audit_status = audit_status,
    scDblFinder_class = scDblFinder.class,
    scDblFinder_score = scDblFinder.score,
    batch_sensitive_flag,
    doublet_sensitive_cluster_flag
  )])
  rownames(add_meta) <- add_dt$expression_cell_id
  clean_obj <- AddMetaData(clean_obj, metadata = add_meta)

  clean_obj <- NormalizeData(clean_obj, verbose = TRUE)
  clean_obj <- FindVariableFeatures(clean_obj, selection.method = "vst", nfeatures = 3000, verbose = TRUE)
  clean_obj <- ScaleData(clean_obj, features = VariableFeatures(clean_obj), verbose = TRUE)
  clean_obj <- RunPCA(clean_obj, features = VariableFeatures(clean_obj), npcs = 50, verbose = TRUE)
  dims_use <- seq_len(min(30, ncol(Embeddings(clean_obj, "pca"))))
  clean_obj <- FindNeighbors(clean_obj, dims = dims_use, verbose = TRUE)
  clean_obj <- FindClusters(clean_obj, resolution = c(0.4, 0.6, 0.8), verbose = TRUE)
  clean_cluster_col <- "RNA_snn_res.0.6"
  if (!clean_cluster_col %in% colnames(clean_obj[[]])) {
    stop("Clean cluster column not found after clustering: ", clean_cluster_col)
  }
  clean_obj[["clean_cluster"]] <- as.character(clean_obj[[clean_cluster_col, drop = TRUE]])
  Idents(clean_obj) <- clean_cluster_col
  clean_obj <- RunUMAP(clean_obj, dims = dims_use, seed.use = 20260707, verbose = TRUE)

  umap <- as.data.table(Embeddings(clean_obj, "umap"), keep.rownames = "expression_cell_id")
  umap_cols <- setdiff(names(umap), "expression_cell_id")
  setnames(umap, umap_cols[1:2], c("umap_1", "umap_2"))
  meta_dt <- as.data.table(clean_obj[[]], keep.rownames = "expression_cell_id")
  meta_dt[, clean_cluster := as.character(get(clean_cluster_col))]
  umap_meta <- merge(umap, meta_dt, by = "expression_cell_id", all.x = TRUE, sort = FALSE)
  fwrite(umap_meta, file.path(output_dir, paste0("gse164522_", set_name, "_umap_coordinates.tsv.gz")), sep = "\t")

  if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
    clean_obj <- JoinLayers(clean_obj)
  }
  Idents(clean_obj) <- clean_cluster_col
  markers <- FindAllMarkers(
    clean_obj,
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  )
  markers_dt <- as.data.table(markers)
  fwrite(markers_dt, file.path(output_dir, paste0("gse164522_", set_name, "_cluster_markers.tsv.gz")), sep = "\t")
  top50 <- markers_dt[order(cluster, p_val_adj, -avg_log2FC), head(.SD, 50), by = cluster]
  fwrite(top50, file.path(output_dir, paste0("gse164522_", set_name, "_cluster_top50_markers.tsv")), sep = "\t")

  write_composition(meta_dt, "tissue", set_name, output_dir)
  write_composition(meta_dt, "sample", set_name, output_dir)
  write_composition(meta_dt, "patient", set_name, output_dir)
  write_composition(meta_dt, "original_cluster_for_cleanup", set_name, output_dir, output_key = "original_cluster")
  write_composition(meta_dt, "cleanup_label", set_name, output_dir)
  write_composition(meta_dt, "celltype_sub", set_name, output_dir)

  run_summary <- data.table(
    set_name = set_name,
    cells = ncol(clean_obj),
    genes = nrow(clean_obj),
    clean_clusters_res_0_6 = length(levels(Idents(clean_obj))),
    original_clusters = uniqueN(meta_dt$original_cluster_for_cleanup),
    cleanup_labels = paste(sort(unique(meta_dt$cleanup_label)), collapse = ";"),
    batch_sensitive_cells = sum(meta_dt$batch_sensitive_flag %in% TRUE),
    doublet_sensitive_cluster_cells = sum(meta_dt$doublet_sensitive_cluster_flag %in% TRUE)
  )
  fwrite(run_summary, file.path(output_dir, paste0("gse164522_", set_name, "_run_summary.tsv")), sep = "\t")

  clean_object_path <- file.path(object_dir, paste0("gse164522_", set_name, "_seurat.rds"))
  saveRDS(clean_obj, clean_object_path, compress = "xz")

  cluster_plot <- DimPlot(clean_obj, reduction = "umap", group.by = clean_cluster_col, label = TRUE, repel = TRUE, pt.size = 0.25) +
    ggtitle(paste0(set_name, ": clean clusters")) +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
  tissue_plot <- DimPlot(clean_obj, reduction = "umap", group.by = "tissue", pt.size = 0.25) +
    ggtitle("Tissue") +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
  cleanup_plot <- DimPlot(clean_obj, reduction = "umap", group.by = "cleanup_label", pt.size = 0.25) +
    ggtitle("Cleanup label") +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
  original_cluster_plot <- DimPlot(clean_obj, reduction = "umap", group.by = "original_cluster_for_cleanup", pt.size = 0.25) +
    ggtitle("Original cluster") +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
  layout_plot <- (cluster_plot | tissue_plot) / (cleanup_plot | original_cluster_plot)
  safe_write_png(layout_plot, file.path(output_dir, paste0("fig1_gse164522_", set_name, "_umap_layout.png")), 12, 9)
  safe_write_png(cluster_plot, file.path(output_dir, paste0("fig1_gse164522_", set_name, "_umap_clusters.png")), 6, 5)
  safe_write_png(tissue_plot, file.path(output_dir, paste0("fig1_gse164522_", set_name, "_umap_tissue.png")), 6.5, 5)
  safe_write_png(cleanup_plot, file.path(output_dir, paste0("fig1_gse164522_", set_name, "_umap_cleanup_label.png")), 7, 5)
  safe_write_png(original_cluster_plot, file.path(output_dir, paste0("fig1_gse164522_", set_name, "_umap_original_cluster.png")), 7, 5)

  if (length(dotplot_features) > 0) {
    dot <- DotPlot(clean_obj, features = dotplot_features, group.by = clean_cluster_col, cols = c("#edf2f7", "#b91c1c")) +
      RotatedAxis() +
      theme_bw(base_size = 8) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.title = element_blank()
      )
    safe_write_png(dot, file.path(output_dir, paste0("fig1_gse164522_", set_name, "_marker_dotplot.png")), 12, 6)
  }

  cluster_tissue <- meta_dt[, .N, by = .(clean_cluster, tissue)]
  cluster_tissue[, fraction := N / sum(N), by = clean_cluster]
  tissue_comp <- ggplot(cluster_tissue, aes(x = clean_cluster, y = fraction, fill = tissue)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.15) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()) +
    labs(y = "Cell fraction", fill = "Tissue")
  safe_write_png(tissue_comp, file.path(output_dir, paste0("fig1_gse164522_", set_name, "_cluster_tissue_composition.png")), 8, 4.5)

  run_summary
}

combined_summary <- rbindlist(
  lapply(names(cell_list_paths), function(set_name) run_one_set(set_name, cell_list_paths[[set_name]])),
  use.names = TRUE,
  fill = TRUE
)
fwrite(combined_summary, file.path(output_dir, "gse164522_residual_tnk_refined_combined_summary.tsv"), sep = "\t")

writeLines(
  c(
    "# GSE164522 residual T/NK refined macrophage/monocyte recluster",
    "",
    paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    paste0("Seurat object: ", rds_path),
    "",
    "Refinement input:",
    paste0(
      refinement_summary$set_name,
      ": source_cells=", refinement_summary$source_cells,
      "; removed_cluster=", refinement_summary$removed_postclean_cluster,
      "; removed_cells=", refinement_summary$removed_cells,
      "; retained_cells=", refinement_summary$retained_cells
    ),
    "",
    "Refined sets:",
    paste0(combined_summary$set_name, ": cells=", combined_summary$cells, "; clean_clusters_res_0_6=", combined_summary$clean_clusters_res_0_6),
    "",
    "No raw GEO files or RDS objects should be committed.",
    "Commit only the packaged upload directory."
  ),
  file.path(output_dir, "STATUS.md")
)

sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
