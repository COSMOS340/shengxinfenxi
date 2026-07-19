#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(harmony)
  library(presto)
})

options(stringsAsFactors = FALSE, future.globals.maxSize = 32 * 1024^3)
future::plan("sequential")
set.seed(340)

gse_root <- "F:/pan-gi-ici-ccc-20260709/GSE205506"
project_dir <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709"
output_dir <- file.path(project_dir, "desktop_exchange/uploads/20260717_gse205506_paper_qc_rebuild")
previous_output_dir <- file.path(project_dir, "desktop_exchange/uploads/20260716_stop_prjna932556_and_build_gse205506")
analysis_dir <- file.path(gse_root, "r_analysis_paper_qc_20260717")
input_path <- file.path(analysis_dir, "gse205506_formal_post_mt_rpca_atlas.rds")
author_path <- file.path(previous_output_dir, "gse205506_author_marker_reference_top50.tsv")

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(output_dir, filename), sep = "\t", quote = FALSE, na = "")
}

log_step <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
  flush.console()
}

if (!file.exists(input_path)) stop("Final RPCA atlas is missing", call. = FALSE)
if (!file.exists(author_path)) stop("Audited author Top50 marker reference is missing", call. = FALSE)

compartment_spec <- data.table(
  compartment = c("T/I/NK", "B", "Myeloid", "Endothelial", "Fibroblast", "Epithelial"),
  slug = c("T_I_NK", "B", "Myeloid", "Endothelial", "Fibroblast", "Epithelial"),
  author_sheet = c(
    "T_I_NK cell markers", "B cell markers", "Myeloid cell markers",
    "Endothelial cell markers", "Fibroblast markers", NA_character_
  )
)
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
working_resolution <- 0.6

broad_marker_sets <- list(
  `T/I/NK` = c("CD3D", "CD3E", "TRAC", "TRBC1"),
  B = c("CD79A", "CD79B", "MS4A1", "TNFRSF17", "MZB1"),
  Myeloid = c("CD14", "CD68"),
  Epithelial = c("EPCAM", "CD24"),
  Fibroblast = c("COL1A2", "COL3A1", "MYH11", "ACTA2"),
  Endothelial = c("VWF", "PECAM1")
)

is_mito_ribosomal <- function(gene) {
  grepl("^MT-", gene) | grepl("^RPS[0-9]", gene) | grepl("^RPL[0-9]", gene)
}

atlas <- readRDS(input_path)
required_compartments <- compartment_spec$compartment
if (!setequal(unique(atlas$paper_final_broad_compartment), required_compartments)) {
  stop("Final atlas does not contain exactly the six required broad compartments", call. = FALSE)
}
if (uniqueN(atlas$geo_accession) != 40L || uniqueN(atlas$geo_subject) != 19L) {
  stop("Final atlas lost sample or subject traceability", call. = FALSE)
}

author_top50 <- fread(author_path)
required_author_columns <- c("compartment_sheet", "cluster", "rank", "gene")
if (!all(required_author_columns %in% names(author_top50))) {
  stop("Audited author marker reference columns do not match the expected schema", call. = FALSE)
}
author_top50[, `:=`(
  compartment_sheet = as.character(compartment_sheet),
  cluster = as.character(cluster),
  gene = as.character(gene)
)]

inventory_rows <- list()
resolution_rows <- list()
for (spec_index in seq_len(nrow(compartment_spec))) {
  compartment_name <- compartment_spec$compartment[[spec_index]]
  slug <- compartment_spec$slug[[spec_index]]
  author_sheet <- compartment_spec$author_sheet[[spec_index]]
  log_step(paste("Reclustering compartment", compartment_name))

  cells <- colnames(atlas)[atlas$paper_final_broad_compartment == compartment_name]
  current <- subset(atlas, cells = cells)
  if ("integrated" %in% Assays(current)) current[["integrated"]] <- NULL
  current@reductions <- list()
  current@graphs <- list()
  current@neighbors <- list()
  DefaultAssay(current) <- "RNA"
  current <- NormalizeData(current, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
  current <- FindVariableFeatures(current, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
  current <- ScaleData(
    current,
    assay = "RNA",
    features = VariableFeatures(current),
    vars.to.regress = "nCount_RNA",
    verbose = TRUE
  )
  current <- RunPCA(
    current,
    assay = "RNA",
    features = VariableFeatures(current),
    npcs = 30,
    approx = TRUE,
    seed.use = 340,
    reduction.name = "pca.compartment",
    reduction.key = "COMPPCA_",
    verbose = TRUE
  )
  current <- RunHarmony(
    current,
    group.by.vars = "geo_accession",
    reduction.use = "pca.compartment",
    dims.use = 1:20,
    reduction.save = "harmony.compartment",
    project.dim = FALSE,
    max.iter.harmony = 20,
    verbose = TRUE
  )
  current <- RunUMAP(
    current,
    reduction = "harmony.compartment",
    dims = 1:20,
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    seed.use = 340,
    reduction.name = "umap.compartment",
    reduction.key = "COMPUMAP_",
    verbose = TRUE
  )
  graph_names <- c(paste0("compartment_", slug, "_nn"), paste0("compartment_", slug, "_snn"))
  current <- FindNeighbors(
    current,
    reduction = "harmony.compartment",
    dims = 1:20,
    k.param = 20,
    graph.name = graph_names,
    verbose = TRUE
  )

  for (resolution in resolutions) {
    resolution_text <- format(resolution, nsmall = 1)
    column_name <- paste0("compartment_cluster_res_", gsub("\\.", "_", resolution_text))
    current <- FindClusters(
      current,
      graph.name = graph_names[[2L]],
      resolution = resolution,
      algorithm = 1,
      n.start = 10,
      n.iter = 10,
      random.seed = 340,
      cluster.name = column_name,
      verbose = TRUE
    )
    resolution_rows[[length(resolution_rows) + 1L]] <- data.table(
      compartment = compartment_name,
      resolution = resolution,
      clusters = uniqueN(current[[column_name, drop = TRUE]]),
      cells = ncol(current),
      working_resolution = resolution == working_resolution
    )
  }

  working_column <- "compartment_cluster_res_0_6"
  current$compartment_working_cluster <- as.character(current[[working_column, drop = TRUE]])
  Idents(current) <- "compartment_working_cluster"

  log_step(paste("Calculating markers for", compartment_name))
  markers <- as.data.table(wilcoxauc(
    current,
    group_by = "compartment_working_cluster",
    assay = "data",
    seurat_assay = "RNA"
  ))
  setnames(markers, c("feature", "group"), c("gene", "cluster"))
  markers[, `:=`(
    cluster = as.character(cluster),
    positive_direction = auc > 0.5 & logFC > 0,
    passes_primary_filter = padj <= 0.05 & auc > 0.5 & logFC > 0 & pct_in >= 5
  )]
  markers[, cluster_numeric := as.integer(cluster)]
  setorder(markers, cluster_numeric, -passes_primary_filter, -positive_direction, -auc, -logFC, -pct_in, gene)
  markers[, rank := seq_len(.N), by = cluster]
  markers[, cluster_numeric := NULL]
  top10 <- markers[rank <= 10L]
  top50 <- markers[rank <= 50L]
  write_tsv(top10, paste0("gse205506_paper_qc_", slug, "_top10_markers.tsv"))
  write_tsv(top50, paste0("gse205506_paper_qc_", slug, "_top50_markers.tsv"))

  current_metadata <- as.data.table(current[[]], keep.rownames = "cell_id")
  distribution <- current_metadata[, .(
    cells = .N,
    samples = uniqueN(geo_accession),
    subjects = uniqueN(geo_subject),
    maximum_sample_fraction = max(table(geo_accession)) / .N,
    maximum_subject_fraction = max(table(geo_subject)) / .N
  ), by = .(cluster = compartment_working_cluster)]

  positive_genes <- broad_marker_sets[[compartment_name]]
  negative_genes <- unique(unlist(broad_marker_sets[names(broad_marker_sets) != compartment_name]))
  evidence <- rbindlist(lapply(sort(unique(current$compartment_working_cluster)), function(cluster_value) {
    top10_genes <- top10[cluster == cluster_value][order(rank), gene]
    top50_genes <- top50[cluster == cluster_value][order(rank), gene]
    positive_evidence <- intersect(top50_genes, positive_genes)
    conflicting_evidence <- intersect(top10_genes, negative_genes)
    mito_ribosomal_genes <- top10_genes[is_mito_ribosomal(top10_genes)]

    if (is.na(author_sheet)) {
      overlap_table <- data.table(
        author_subtype = NA_character_,
        overlap_count = NA_integer_,
        overlap_genes = NA_character_,
        jaccard = NA_real_
      )
    } else {
      reference <- author_top50[compartment_sheet == author_sheet]
      overlap_table <- reference[, .(
        author_genes = list(unique(gene))
      ), by = .(author_subtype = cluster)]
      overlap_table[, `:=`(
        overlap_genes = vapply(author_genes, function(genes) paste(intersect(top50_genes, genes), collapse = ";"), character(1)),
        overlap_count = vapply(author_genes, function(genes) length(intersect(top50_genes, genes)), integer(1)),
        jaccard = vapply(author_genes, function(genes) {
          union_size <- length(union(top50_genes, genes))
          if (union_size == 0L) NA_real_ else length(intersect(top50_genes, genes)) / union_size
        }, numeric(1))
      )]
      setorder(overlap_table, -overlap_count, -jaccard, author_subtype)
    }

    top_author <- overlap_table[1L]
    second_overlap <- if (nrow(overlap_table) >= 2L) overlap_table$overlap_count[[2L]] else NA_integer_
    author_label_pass <- !is.na(author_sheet) &&
      !is.na(top_author$overlap_count[[1L]]) && top_author$overlap_count[[1L]] >= 5L &&
      length(mito_ribosomal_genes) == 0L && length(conflicting_evidence) == 0L
    reviewed_label <- if (is.na(author_sheet)) {
      "Epithelial_no_author_sheet"
    } else if (author_label_pass) {
      top_author$author_subtype[[1L]]
    } else {
      "Unresolved_review"
    }

    data.table(
      compartment = compartment_name,
      cluster = cluster_value,
      top10_markers = paste(top10_genes, collapse = ";"),
      top50_markers = paste(top50_genes, collapse = ";"),
      canonical_positive_marker_evidence = paste(positive_evidence, collapse = ";"),
      canonical_negative_or_conflicting_evidence = paste(conflicting_evidence, collapse = ";"),
      top10_mito_ribosomal_markers = paste(mito_ribosomal_genes, collapse = ";"),
      top10_mito_ribosomal_count = length(mito_ribosomal_genes),
      top10_conflicting_broad_marker_count = length(conflicting_evidence),
      author_sheet = author_sheet,
      top_author_subtype = top_author$author_subtype[[1L]],
      top_author_overlap_count = top_author$overlap_count[[1L]],
      second_author_overlap_count = second_overlap,
      top_author_overlap_genes = top_author$overlap_genes[[1L]],
      top_author_jaccard = top_author$jaccard[[1L]],
      author_label_rule_pass = author_label_pass,
      final_reviewed_label = reviewed_label,
      review_flag = fifelse(
        is.na(author_sheet), "no_author_epithelial_sheet",
        fifelse(author_label_pass, "author_overlap_supported", "manual_review_required")
      ),
      action = "retain cluster; no whole-cluster deletion"
    )
  }), use.names = TRUE, fill = TRUE)
  evidence <- merge(evidence, distribution, by = "cluster", all.x = TRUE, sort = FALSE)
  evidence[, cluster_numeric := as.integer(cluster)]
  setorder(evidence, cluster_numeric)
  evidence[, cluster_numeric := NULL]

  if (!is.na(author_sheet)) {
    overlap_output <- rbindlist(lapply(sort(unique(current$compartment_working_cluster)), function(cluster_value) {
      genes <- top50[cluster == cluster_value][order(rank), gene]
      reference <- author_top50[compartment_sheet == author_sheet]
      rows <- reference[, .(author_genes = list(unique(gene))), by = .(author_subtype = cluster)]
      rows[, `:=`(
        compartment = compartment_name,
        cluster = cluster_value,
        overlap_count = vapply(author_genes, function(author_genes) length(intersect(genes, author_genes)), integer(1)),
        overlap_genes = vapply(author_genes, function(author_genes) paste(intersect(genes, author_genes), collapse = ";"), character(1)),
        jaccard = vapply(author_genes, function(author_genes) {
          union_size <- length(union(genes, author_genes))
          if (union_size == 0L) NA_real_ else length(intersect(genes, author_genes)) / union_size
        }, numeric(1))
      )]
      rows[, author_genes := NULL]
      setorder(rows, -overlap_count, -jaccard, author_subtype)
      rows[, author_overlap_rank := seq_len(.N)]
      rows
    }))
  } else {
    overlap_output <- evidence[, .(
      compartment,
      cluster,
      author_subtype = NA_character_,
      overlap_count = NA_integer_,
      overlap_genes = NA_character_,
      jaccard = NA_real_,
      author_overlap_rank = NA_integer_,
      note = "mmc3.xlsx has no epithelial marker sheet"
    )]
  }
  write_tsv(overlap_output, paste0("gse205506_paper_qc_", slug, "_author_overlap.tsv"))
  write_tsv(evidence, paste0("gse205506_paper_qc_", slug, "_annotation_evidence.tsv"))

  label_map <- setNames(evidence$final_reviewed_label, evidence$cluster)
  current$compartment_reviewed_label <- unname(label_map[current$compartment_working_cluster])
  if (anyNA(current$compartment_reviewed_label)) stop("Compartment label mapping failed for ", compartment_name, call. = FALSE)
  current@misc$compartment_reclustering <- list(
    compartment = compartment_name,
    integration = "Harmony project implementation by exact geo_accession after RNA PCA",
    pca_dims = 1:20,
    resolutions = resolutions,
    working_resolution = working_resolution,
    annotation_rule = "author Top50 overlap >=5, no Top10 mitochondrial/ribosomal marker, and no Top10 conflicting broad marker; otherwise Unresolved_review",
    cluster_deletion = "none"
  )

  object_path <- file.path(analysis_dir, paste0("gse205506_compartment_", slug, ".rds"))
  saveRDS(current, object_path, compress = FALSE)
  if (!requireNamespace("digest", quietly = TRUE)) stop("R package digest is required", call. = FALSE)
  inventory_rows[[length(inventory_rows) + 1L]] <- data.table(
    compartment = compartment_name,
    object_path = normalizePath(object_path, winslash = "/", mustWork = TRUE),
    size_bytes = as.numeric(file.info(object_path)$size),
    sha256 = digest::digest(file = object_path, algo = "sha256", serialize = FALSE),
    cells = ncol(current),
    clusters_at_resolution_0_6 = uniqueN(current$compartment_working_cluster),
    reviewed_labels = uniqueN(current$compartment_reviewed_label),
    delivery = "local_only"
  )
  rm(current, markers, top10, top50, evidence, overlap_output)
  invisible(gc())
}

write_tsv(rbindlist(resolution_rows), "gse205506_paper_qc_compartment_resolution_summary.tsv")
write_tsv(rbindlist(inventory_rows), "gse205506_paper_qc_compartment_object_inventory.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_gse205506_compartment_reclustering.txt"))
log_step("All six compartment objects complete; no whole cluster deleted")
