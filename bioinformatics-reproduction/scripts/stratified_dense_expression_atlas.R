#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(irlba)
  library(matrixStats)
  library(uwot)
  library(RANN)
  library(igraph)
  library(mclust)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript stratified_dense_expression_atlas.R <raw_dir> <output_dir> [file_prefix]")
}

raw_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
file_prefix <- if (length(args) >= 3) args[[3]] else ""

table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
object_dir <- file.path(output_dir, "objects")
log_path <- file.path(output_dir, "dense_expression_atlas_log.tsv")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260702)

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

dataset_from_path <- function(path, suffix) {
  name <- basename(path)
  name <- sub(paste0("^", file_prefix), "", name)
  sub(suffix, "", name)
}

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, limitsize = FALSE)
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, limitsize = FALSE)
}

metadata_files <- sort(list.files(raw_dir, pattern = "_metadata\\.csv\\.gz$", full.names = TRUE))
expression_files <- sort(list.files(raw_dir, pattern = "_normalized_expression\\.csv\\.gz$", full.names = TRUE))
required_metadata_columns <- c("index", "cancer", "tissue", "MajorCluster", "UMAP1", "UMAP2")

if (length(metadata_files) == 0 || length(expression_files) == 0) {
  stop("No matching metadata/expression csv.gz files found under ", raw_dir)
}

metadata_list <- lapply(metadata_files, function(path) {
  dataset <- dataset_from_path(path, "_metadata\\.csv\\.gz$")
  dt <- fread(path)
  missing <- setdiff(required_metadata_columns, names(dt))
  if (length(missing) > 0) {
    stop(basename(path), " missing columns: ", paste(missing, collapse = ", "))
  }
  dt[, dataset := dataset]
  dt[, cell_key := paste(dataset, index, sep = "||")]
  dt
})
metadata <- rbindlist(metadata_list, fill = TRUE)
metadata <- metadata[!is.na(UMAP1) & !is.na(UMAP2) & !is.na(MajorCluster) & !is.na(cancer)]

gene_sets <- list()
inspection <- list()
for (path in expression_files) {
  dataset <- dataset_from_path(path, "_normalized_expression\\.csv\\.gz$")
  header <- names(fread(cmd = paste("gzip -dc", shQuote(path), "| head -n 1"), nrows = 0))
  if (length(header) < 2 || header[[1]] != "index") {
    stop(basename(path), " must have first expression column named index")
  }
  genes <- header[-1]
  gene_sets[[dataset]] <- genes
  inspection[[dataset]] <- data.table(dataset = dataset, expression_file = path, expression_genes = length(genes))
}

common_genes <- Reduce(intersect, gene_sets)
fwrite(rbindlist(inspection), file.path(table_dir, "dense_expression_file_inspection.tsv"), sep = "\t")
writeLines(common_genes, file.path(table_dir, "common_genes.txt"))

metadata[, stratum := paste(cancer, MajorCluster, sep = "||")]
target_per_stratum <- 80L
max_cells <- 8000L
selected <- metadata[, .SD[sample.int(.N, min(.N, target_per_stratum))], by = stratum]
if (nrow(selected) > max_cells) {
  selected <- selected[sample.int(.N, max_cells)]
}
selected[, stratum := NULL]
setorder(selected, dataset, index)
fwrite(selected, file.path(table_dir, "selected_metadata.tsv"), sep = "\t")

read_selected_expression <- function(path) {
  dataset_label <- dataset_from_path(path, "_normalized_expression\\.csv\\.gz$")
  selected_file <- selected[dataset == dataset_label, .(index, cell_key)]
  if (nrow(selected_file) == 0) {
    return(NULL)
  }
  write_log("read_expression_start", "running", paste0(dataset_label, "; selected_cells=", nrow(selected_file)))
  dt <- fread(cmd = paste("gzip -dc", shQuote(path)), select = c("index", common_genes), nThread = 2, showProgress = FALSE)
  dt <- dt[index %in% selected_file$index]
  dt <- merge(selected_file, dt, by = "index", sort = FALSE)
  mat <- as.matrix(dt[, ..common_genes])
  storage.mode(mat) <- "double"
  rownames(mat) <- dt$cell_key
  rm(dt)
  gc()
  mat
}

expr_list <- lapply(expression_files, read_selected_expression)
expr_list <- expr_list[!vapply(expr_list, is.null, logical(1))]
expr <- do.call(rbind, expr_list)
selected <- selected[match(rownames(expr), cell_key)]

gene_vars <- matrixStats::colVars(expr)
names(gene_vars) <- colnames(expr)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(3000L, sum(gene_vars > 0)))]
scaled <- scale(expr[, top_genes, drop = FALSE])
scaled[!is.finite(scaled)] <- 0
rm(expr)
gc()

pca <- irlba::prcomp_irlba(scaled, n = 30, center = FALSE, scale. = FALSE)
umap_coords <- uwot::umap(pca$x, n_neighbors = 30, min_dist = 0.3, metric = "cosine", n_threads = 2, verbose = TRUE)
knn <- RANN::nn2(pca$x, k = 31)$nn.idx[, -1, drop = FALSE]
edges <- data.table(from = rep(seq_len(nrow(knn)), times = ncol(knn)), to = as.vector(knn))
graph <- igraph::make_empty_graph(n = nrow(pca$x), directed = FALSE)
graph <- igraph::add_edges(graph, as.vector(t(as.matrix(edges[from != to, .(from, to)]))))
graph <- igraph::simplify(graph)
clusters <- paste0("E", igraph::cluster_louvain(graph)$membership)

embedding <- data.table(
  cell_key = selected$cell_key,
  index = selected$index,
  dataset = selected$dataset,
  cancer = selected$cancer,
  tissue = selected$tissue,
  MajorCluster = selected$MajorCluster,
  AuthorUMAP1 = as.numeric(selected$UMAP1),
  AuthorUMAP2 = as.numeric(selected$UMAP2),
  ExprUMAP1 = umap_coords[, 1],
  ExprUMAP2 = umap_coords[, 2],
  expression_louvain = clusters
)
fwrite(embedding, file.path(table_dir, "expression_umap_coordinates.tsv"), sep = "\t")

metrics <- data.table(
  metric = c("selected_cells", "common_genes", "top_variable_genes", "expression_louvain_clusters",
             "ari_expression_louvain_vs_majorcluster", "ari_expression_louvain_vs_cancer", "ari_expression_louvain_vs_tissue"),
  value = c(nrow(embedding), length(common_genes), length(top_genes), uniqueN(embedding$expression_louvain),
            mclust::adjustedRandIndex(embedding$expression_louvain, embedding$MajorCluster),
            mclust::adjustedRandIndex(embedding$expression_louvain, embedding$cancer),
            mclust::adjustedRandIndex(embedding$expression_louvain, embedding$tissue))
)
fwrite(metrics, file.path(table_dir, "expression_first_pass_metrics.tsv"), sep = "\t")
saveRDS(list(embedding = embedding, metrics = metrics, top_genes = top_genes), file.path(object_dir, "dense_expression_atlas_embedding.rds"))

theme_atlas <- theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
p_major <- ggplot(embedding, aes(ExprUMAP1, ExprUMAP2, color = MajorCluster)) +
  geom_point(size = 0.35, alpha = 0.75) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 2, byrow = TRUE)) +
  labs(title = "Expression-derived UMAP by reference label", x = "Expression UMAP1", y = "Expression UMAP2", color = "Reference label") +
  theme_atlas + theme(legend.text = element_text(size = 7), legend.key.height = grid::unit(0.26, "cm"))
save_plot(p_major, "expression_umap_by_reference_label", 10.8, 7.2)

p_cluster <- ggplot(embedding, aes(ExprUMAP1, ExprUMAP2, color = expression_louvain)) +
  geom_point(size = 0.35, alpha = 0.8) +
  labs(title = "Expression-derived UMAP by Louvain cluster", x = "Expression UMAP1", y = "Expression UMAP2", color = "Expr cluster") +
  theme_atlas
save_plot(p_cluster, "expression_umap_by_louvain_cluster", 8.2, 6.4)

write_log("script_done", "complete", paste0("selected_cells=", nrow(embedding), "; clusters=", uniqueN(embedding$expression_louvain)))
message("Dense expression atlas first pass complete.")
