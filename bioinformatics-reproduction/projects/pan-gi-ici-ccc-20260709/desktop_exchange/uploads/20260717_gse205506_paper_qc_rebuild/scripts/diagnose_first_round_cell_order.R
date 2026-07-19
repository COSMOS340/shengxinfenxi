#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))

light_path <- "F:/pan-gi-ici-ccc-20260709/GSE205506/r_analysis_paper_qc_20260717/gse205506_first_round_rpca_light_checkpoint.rds"
basic_path <- "F:/pan-gi-ici-ccc-20260709/GSE205506/r_analysis_paper_qc_20260717/gse205506_article_basic_qc_primary.rds"

cat("loading_light\n")
light <- readRDS(light_path)$object
cat("loading_basic\n")
basic <- readRDS(basic_path)

light_cells <- colnames(light)
basic_cells <- colnames(basic)
cat(sprintf("light_cells=%d basic_cells=%d\n", length(light_cells), length(basic_cells)))
cat(sprintf(
  "identical_order=%s setequal=%s\n",
  identical(light_cells, basic_cells),
  setequal(light_cells, basic_cells)
))
cat(sprintf(
  "light_only=%d basic_only=%d same_positions=%d\n",
  length(setdiff(light_cells, basic_cells)),
  length(setdiff(basic_cells, light_cells)),
  sum(light_cells == basic_cells)
))
cat("light_first5=", paste(head(light_cells, 5L), collapse = ";"), "\n", sep = "")
cat("basic_first5=", paste(head(basic_cells, 5L), collapse = ";"), "\n", sep = "")

for (graph_name in c("rpca_first_nn", "rpca_first_snn")) {
  graph_cells <- colnames(light[[graph_name]])
  cat(sprintf(
    "graph=%s cells=%d identical_light=%s setequal_light=%s\n",
    graph_name,
    length(graph_cells),
    identical(graph_cells, light_cells),
    setequal(graph_cells, light_cells)
  ))
}

reordered <- basic[, light_cells]
cat(sprintf("bracket_reorder_identical=%s\n", identical(colnames(reordered), light_cells)))

for (graph_name in c("rpca_first_nn", "rpca_first_snn")) {
  graph_reordered <- light[[graph_name]][basic_cells, basic_cells, drop = FALSE]
  cat(sprintf(
    "graph_subset=%s class=%s identical_basic=%s\n",
    graph_name,
    paste(class(graph_reordered), collapse = ";"),
    identical(colnames(graph_reordered), basic_cells) && identical(rownames(graph_reordered), basic_cells)
  ))
  graph_reordered <- as(graph_reordered, "Graph")
  probe <- basic
  probe[[graph_name]] <- graph_reordered
  cat(sprintf(
    "graph_assignment=%s class=%s identical_basic=%s\n",
    graph_name,
    paste(class(probe[[graph_name]]), collapse = ";"),
    identical(colnames(probe[[graph_name]]), basic_cells) &&
      identical(rownames(probe[[graph_name]]), basic_cells)
  ))
  rm(graph_reordered, probe)
  invisible(gc())
}

probe <- basic
probe[["integrated"]] <- light[["integrated"]]
cat(sprintf(
  "integrated_assignment_identical_basic=%s\n",
  identical(colnames(probe[["integrated"]]), basic_cells)
))
for (reduction_name in c("pca.rpca.first", "umap.rpca.first")) {
  probe[[reduction_name]] <- light[[reduction_name]]
  cat(sprintf(
    "reduction_assignment=%s identical_basic=%s\n",
    reduction_name,
    identical(rownames(Embeddings(probe[[reduction_name]])), basic_cells)
  ))
}
