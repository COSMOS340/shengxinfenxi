suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

out <- "I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-cancer-mono-mac/fig4_snatac/outputs_fast_merge"
tab_dir <- file.path(out, "tables")
fig_dir <- file.path(out, "figures")
full_path <- file.path(tab_dir, "fig4c_chipseeker_peak_annotation_full.tsv.gz")
count_path <- file.path(tab_dir, "fig4c_chipseeker_annotation_counts.tsv")

anno_dt <- fread(full_path)
if (!"annotation" %in% names(anno_dt)) {
  stop("Missing annotation column in ", full_path)
}

anno_dt[, broad_annotation := fifelse(
  startsWith(annotation, "Promoter"),
  "Promoter",
  fifelse(
    startsWith(annotation, "Intron"),
    "Intron",
    fifelse(
      startsWith(annotation, "Exon"),
      "Exon",
      fifelse(
        startsWith(annotation, "Downstream"),
        "Downstream",
        annotation
      )
    )
  )
)]

detailed_counts_path <- file.path(tab_dir, "fig4c_chipseeker_annotation_counts_detailed.tsv")
if (file.exists(count_path) && !file.exists(detailed_counts_path)) {
  file.copy(count_path, detailed_counts_path, overwrite = FALSE)
}

anno_counts <- anno_dt[, .N, by = .(annotation = broad_annotation)][order(-N)]
anno_counts[, fraction := N / sum(N)]
anno_counts[, label := paste0(annotation, " (", percent(fraction, accuracy = 0.1), ")")]
anno_counts[, annotation := factor(annotation, levels = annotation)]

fwrite(anno_dt, full_path, sep = "\t", quote = FALSE)
fwrite(anno_counts[, .(annotation = as.character(annotation), N, fraction)], count_path, sep = "\t", quote = FALSE)

p_pie <- ggplot(anno_counts, aes(x = "", y = fraction, fill = annotation)) +
  geom_col(width = 1, color = "white", linewidth = 0.25) +
  coord_polar(theta = "y") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  guides(fill = guide_legend(ncol = 1, title = "Annotation")) +
  labs(title = "Fig. 4C ChIPSeeker annotation", x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

ggsave(file.path(fig_dir, "fig4c_chipseeker_annotation_pie.png"), p_pie, width = 6.2, height = 4.8, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(fig_dir, "fig4c_chipseeker_annotation_pie.pdf"), p_pie, width = 6.2, height = 4.8, bg = "white", limitsize = FALSE)

cat("fig4c_broad_annotation_regenerated\n")
print(anno_counts[, .(annotation = as.character(annotation), N, fraction)])
