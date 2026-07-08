options(stringsAsFactors = FALSE)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Required R package missing: Matrix")
}

find_repo_root <- function(start) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, ".git"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Cannot find repository root")
    }
    current <- parent
  }
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  script_path <- normalizePath("bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_working_qc_sparse_objects/scripts/pc_working_qc_sparse_objects.R", winslash = "/", mustWork = TRUE)
}

request_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- find_repo_root(request_dir)
raw_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche", "desktop_exchange", "requests", "20260708_core_geo_download", "raw_core_geo")
object_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche", "desktop_exchange", "requests", "20260708_core_geo_download", "derived_working_qc_sparse_objects_20260708")
upload_dir <- file.path(repo_root, "bioinformatics-reproduction", "projects", "crlm-caf-ecm-spatial-niche", "desktop_exchange", "uploads", "20260708_working_qc_sparse_objects")
input_metadata <- file.path(request_dir, "inputs", "sample_metadata_from_geo_series.tsv")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

read_gz_table_base <- function(path, header = TRUE) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, header = header, check.names = FALSE, sep = "\t")
}

read_lines_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con)
}

safe_feature_symbols <- function(feature_df) {
  if (ncol(feature_df) >= 2) {
    return(as.character(feature_df[[2]]))
  }
  as.character(feature_df[[1]])
}

mt_percent_from_values <- function(mat, values) {
  mt <- grepl("^MT-", values)
  total <- Matrix::colSums(mat)
  mt_counts <- if (any(mt)) Matrix::colSums(mat[mt, , drop = FALSE]) else rep(0, ncol(mat))
  ifelse(total > 0, mt_counts / total * 100, NA_real_)
}

file_size_or_na <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  as.numeric(file.info(path)$size)
}

object_manifest_rows <- list()
status_lines <- c()
sample_meta <- read.delim(input_metadata, sep = "\t", header = TRUE, check.names = FALSE)

thresholds <- data.frame(
  dataset_id = c("GSE225857_nonimmune", "GSE178318", "GSE245552"),
  filter_scope = c("author_metadata_only", "working_qc", "working_qc"),
  min_features = c("", "500", "500"),
  min_counts = c("", "", "1000"),
  max_percent_mt = c("", "15", ""),
  notes = c(
    "Use exact author non-immune metadata and source clusters; no re-filtering in this request.",
    "Working checkpoint from sensitivity table; not final inclusion.",
    "Mitochondrial percentage excluded because feature files differ by tissue group."
  )
)
write_tsv(thresholds, file.path(upload_dir, "working_qc_thresholds.tsv"))

marker_sets <- data.frame(
  marker_group = c(
    "fibroblast_ecm", "fibroblast_ecm", "fibroblast_ecm", "fibroblast_ecm", "fibroblast_ecm",
    "caf_activation", "caf_activation", "caf_activation", "caf_activation",
    "pericyte_like", "pericyte_like", "pericyte_like",
    "endothelial", "endothelial", "endothelial",
    "epithelial_tumor", "epithelial_tumor", "epithelial_tumor",
    "immune_context", "immune_context", "immune_context"
  ),
  gene_symbol = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM",
    "FAP", "ACTA2", "POSTN", "THBS2",
    "MCAM", "RGS5", "PDGFRB",
    "PECAM1", "VWF", "KDR",
    "EPCAM", "KRT8", "KRT18",
    "PTPRC", "CD3D", "LYZ"
  )
)

gse225857_tar <- file.path(raw_dir, "GSE225857", "GSE225857_RAW.tar")
gse225857_temp <- tempfile("gse225857_working_qc_")
dir.create(gse225857_temp)
utils::untar(gse225857_tar, files = c("GSM7058755_non_immune_meta.txt.gz"), exdir = gse225857_temp)
nonimmune_meta <- read_gz_table_base(file.path(gse225857_temp, "GSM7058755_non_immune_meta.txt.gz"), header = TRUE)

cluster_values <- as.character(nonimmune_meta[["cluster"]])
cluster_prefix <- sub("[0-9_].*$", "", cluster_values)
cluster_summary <- aggregate(
  rep(1L, length(cluster_values)),
  by = list(cluster = cluster_values, source_cluster_prefix = cluster_prefix),
  FUN = sum
)
names(cluster_summary)[3] <- "cell_count"
cluster_summary <- cluster_summary[order(cluster_summary$source_cluster_prefix, -cluster_summary$cell_count, cluster_summary$cluster), ]
write_tsv(cluster_summary, file.path(upload_dir, "gse225857_nonimmune_source_cluster_summary.tsv"))

fibroblast_source <- grepl("fibroblast|fibrblast", cluster_values)
fibroblast_summary <- aggregate(
  rep(1L, sum(fibroblast_source)),
  by = list(
    cluster = cluster_values[fibroblast_source],
    organs = as.character(nonimmune_meta[["organs"]][fibroblast_source]),
    patients_organ = as.character(nonimmune_meta[["patients_organ"]][fibroblast_source])
  ),
  FUN = sum
)
names(fibroblast_summary)[4] <- "cell_count"
fibroblast_summary <- fibroblast_summary[order(fibroblast_summary$cluster, fibroblast_summary$organs, fibroblast_summary$patients_organ), ]
write_tsv(fibroblast_summary, file.path(upload_dir, "gse225857_nonimmune_fibroblast_source_clusters.tsv"))
status_lines <- c(status_lines, sprintf("GSE225857 nonimmune metadata cells reviewed: %s", nrow(nonimmune_meta)))

gse178318_dir <- file.path(raw_dir, "GSE178318")
gse178318_matrix <- Matrix::readMM(gzfile(file.path(gse178318_dir, "GSE178318_matrix.mtx.gz"), open = "rt"))
gse178318_features <- read_gz_table_base(file.path(gse178318_dir, "GSE178318_genes.tsv.gz"), header = FALSE)
gse178318_barcodes <- read_lines_gz(file.path(gse178318_dir, "GSE178318_barcodes.tsv.gz"))
gse178318_symbols <- safe_feature_symbols(gse178318_features)
gse178318_suffix <- sub("^[^_]+_", "", gse178318_barcodes)
gse178318_qc <- data.frame(
  barcode = gse178318_barcodes,
  barcode_suffix = gse178318_suffix,
  nCount_RNA = as.numeric(Matrix::colSums(gse178318_matrix)),
  nFeature_RNA = as.numeric(Matrix::colSums(gse178318_matrix > 0)),
  percent.mt = mt_percent_from_values(gse178318_matrix, gse178318_symbols)
)
gse178318_keep <- gse178318_qc$nFeature_RNA >= 500 & gse178318_qc$percent.mt <= 15
gse178318_retention <- aggregate(
  data.frame(start_cells = rep(1L, nrow(gse178318_qc)), retained_cells = as.integer(gse178318_keep)),
  by = list(barcode_suffix = gse178318_qc$barcode_suffix),
  FUN = sum
)
gse178318_retention$retained_fraction <- gse178318_retention$retained_cells / gse178318_retention$start_cells
gse178318_retention$min_features <- 500
gse178318_retention$max_percent_mt <- 15
gse178318_retention <- gse178318_retention[order(gse178318_retention$barcode_suffix), ]
write_tsv(gse178318_retention, file.path(upload_dir, "gse178318_working_qc_retention_by_suffix.tsv"))

gse178318_object_path <- file.path(object_dir, "gse178318_working_qc_sparse_matrix.rds")
saveRDS(
  list(
    dataset_id = "GSE178318",
    thresholds = thresholds[thresholds$dataset_id == "GSE178318", ],
    counts = gse178318_matrix[, gse178318_keep, drop = FALSE],
    features = gse178318_features,
    metadata = gse178318_qc[gse178318_keep, , drop = FALSE]
  ),
  gse178318_object_path
)
object_manifest_rows[[length(object_manifest_rows) + 1]] <- data.frame(
  dataset_id = "GSE178318",
  object_label = "gse178318_working_qc_sparse_matrix",
  object_path = gse178318_object_path,
  file_size_bytes = file_size_or_na(gse178318_object_path),
  features = nrow(gse178318_matrix),
  cells = sum(gse178318_keep)
)

png(file.path(upload_dir, "gse178318_working_qc_retention_by_suffix.png"), width = 1800, height = 900, res = 150)
par(mar = c(8, 5, 3, 1))
barplot(
  gse178318_retention$retained_fraction,
  names.arg = gse178318_retention$barcode_suffix,
  las = 2,
  ylim = c(0, 1),
  ylab = "Retained fraction",
  main = "GSE178318 working QC retention by barcode suffix",
  col = "#4C78A8",
  border = NA
)
abline(h = c(0.8, 0.9), lty = 2, col = "gray50")
dev.off()

marker_presence_rows <- list()
for (i in seq_len(nrow(marker_sets))) {
  marker_presence_rows[[length(marker_presence_rows) + 1]] <- data.frame(
    dataset_id = "GSE178318",
    sample_label = "all",
    tissue_site = "",
    marker_group = marker_sets$marker_group[i],
    gene_symbol = marker_sets$gene_symbol[i],
    present = marker_sets$gene_symbol[i] %in% gse178318_symbols
  )
}

status_lines <- c(
  status_lines,
  sprintf("GSE178318 retained cells: %s of %s", sum(gse178318_keep), ncol(gse178318_matrix))
)
rm(gse178318_matrix)
gc()

gse245552_tar <- file.path(raw_dir, "GSE245552", "GSE245552_RAW.tar")
tar_members <- utils::untar(gse245552_tar, list = TRUE)
matrix_members <- sort(tar_members[grepl("_matrix[.]mtx[.]gz$", tar_members)])
gse245552_temp <- tempfile("gse245552_working_qc_")
dir.create(gse245552_temp)
gse245552_retention_rows <- list()

for (matrix_member in matrix_members) {
  prefix <- sub("_matrix[.]mtx[.]gz$", "", matrix_member)
  barcode_member <- paste0(prefix, "_barcodes.tsv.gz")
  feature_member <- paste0(prefix, "_features.tsv.gz")
  utils::untar(gse245552_tar, files = c(matrix_member, barcode_member, feature_member), exdir = gse245552_temp)

  matrix_path <- file.path(gse245552_temp, matrix_member)
  barcode_path <- file.path(gse245552_temp, barcode_member)
  feature_path <- file.path(gse245552_temp, feature_member)
  mat <- Matrix::readMM(gzfile(matrix_path, open = "rt"))
  features <- read_gz_table_base(feature_path, header = FALSE)
  barcodes <- read_lines_gz(barcode_path)
  symbols <- safe_feature_symbols(features)

  sample_label <- sub("^GSM[0-9]+_", "", prefix)
  gsm <- sub("_.*$", "", prefix)
  meta_row <- sample_meta[sample_meta$sample_id == gsm, , drop = FALSE]
  tissue_site <- if (nrow(meta_row) > 0) meta_row$tissue_site[1] else ""
  lesion_type <- if (nrow(meta_row) > 0) meta_row$lesion_type[1] else ""

  sample_qc <- data.frame(
    barcode = barcodes,
    sample_label = sample_label,
    gsm = gsm,
    tissue_site = tissue_site,
    lesion_type = lesion_type,
    nCount_RNA = as.numeric(Matrix::colSums(mat)),
    nFeature_RNA = as.numeric(Matrix::colSums(mat > 0))
  )
  keep <- sample_qc$nFeature_RNA >= 500 & sample_qc$nCount_RNA >= 1000

  gse245552_retention_rows[[length(gse245552_retention_rows) + 1]] <- data.frame(
    dataset_id = "GSE245552",
    gsm = gsm,
    sample_label = sample_label,
    tissue_site = tissue_site,
    lesion_type = lesion_type,
    feature_rows = nrow(features),
    start_cells = ncol(mat),
    retained_cells = sum(keep),
    retained_fraction = sum(keep) / ncol(mat),
    min_features = 500,
    min_counts = 1000,
    mt_upper_prefix_count_col1 = sum(grepl("^MT-", as.character(features[[1]]))),
    mt_upper_prefix_count_col2 = if (ncol(features) >= 2) sum(grepl("^MT-", as.character(features[[2]]))) else NA_integer_
  )

  sample_object_path <- file.path(object_dir, paste0("gse245552_", sample_label, "_working_qc_sparse_matrix.rds"))
  saveRDS(
    list(
      dataset_id = "GSE245552",
      sample_label = sample_label,
      gsm = gsm,
      tissue_site = tissue_site,
      lesion_type = lesion_type,
      thresholds = thresholds[thresholds$dataset_id == "GSE245552", ],
      counts = mat[, keep, drop = FALSE],
      features = features,
      metadata = sample_qc[keep, , drop = FALSE]
    ),
    sample_object_path
  )
  object_manifest_rows[[length(object_manifest_rows) + 1]] <- data.frame(
    dataset_id = "GSE245552",
    object_label = paste0("gse245552_", sample_label, "_working_qc_sparse_matrix"),
    object_path = sample_object_path,
    file_size_bytes = file_size_or_na(sample_object_path),
    features = nrow(mat),
    cells = sum(keep)
  )

  for (i in seq_len(nrow(marker_sets))) {
    marker_presence_rows[[length(marker_presence_rows) + 1]] <- data.frame(
      dataset_id = "GSE245552",
      sample_label = sample_label,
      tissue_site = tissue_site,
      marker_group = marker_sets$marker_group[i],
      gene_symbol = marker_sets$gene_symbol[i],
      present = marker_sets$gene_symbol[i] %in% symbols
    )
  }

  rm(mat)
  unlink(c(matrix_path, barcode_path, feature_path))
  gc()
}

gse245552_retention <- do.call(rbind, gse245552_retention_rows)
gse245552_retention <- gse245552_retention[order(gse245552_retention$tissue_site, gse245552_retention$sample_label), ]
write_tsv(gse245552_retention, file.path(upload_dir, "gse245552_working_qc_retention_by_sample.tsv"))

png(file.path(upload_dir, "gse245552_working_qc_retention_by_sample.png"), width = 2200, height = 1000, res = 150)
par(mar = c(9, 5, 3, 1))
bar_cols <- ifelse(gse245552_retention$tissue_site == "liver metastasis", "#F58518",
  ifelse(gse245552_retention$tissue_site == "primary tumor", "#4C78A8",
    ifelse(gse245552_retention$tissue_site == "colon adjacent tissue", "#54A24B", "#B279A2")))
barplot(
  gse245552_retention$retained_fraction,
  names.arg = gse245552_retention$sample_label,
  las = 2,
  ylim = c(0, 1),
  ylab = "Retained fraction",
  main = "GSE245552 working QC retention by sample",
  col = bar_cols,
  border = NA
)
abline(h = c(0.7, 0.8, 0.9), lty = 2, col = "gray50")
legend("bottomleft", legend = c("primary tumor", "liver metastasis", "colon adjacent tissue", "liver adjacent tissue"),
  fill = c("#4C78A8", "#F58518", "#54A24B", "#B279A2"), bty = "n", cex = 0.8)
dev.off()

object_manifest <- do.call(rbind, object_manifest_rows)
write_tsv(object_manifest, file.path(upload_dir, "working_qc_object_manifest.tsv"))

marker_presence <- do.call(rbind, marker_presence_rows)
marker_presence <- marker_presence[order(marker_presence$dataset_id, marker_presence$sample_label, marker_presence$marker_group, marker_presence$gene_symbol), ]
write_tsv(marker_presence, file.path(upload_dir, "caf_ecm_marker_gene_presence.tsv"))

status_lines <- c(
  status_lines,
  sprintf("GSE245552 retained cells: %s of %s", sum(gse245552_retention$retained_cells), sum(gse245552_retention$start_cells)),
  sprintf("Working object manifest rows: %s", nrow(object_manifest)),
  sprintf("Marker presence rows: %s", nrow(marker_presence))
)

status_text <- c(
  "# CRLM CAF/ECM working QC sparse objects",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  paste0("Raw input directory: `", raw_dir, "`"),
  paste0("PC-local object directory: `", object_dir, "`"),
  "",
  status_lines,
  "",
  "No final filtering, integration, clustering, marker testing, or cell-type annotation was performed.",
  "RDS objects are PC-local and should not be committed to GitHub."
)
writeLines(status_text, file.path(upload_dir, "STATUS.md"))
