# Formal Fig4A GSE306459 snATAC rerun status

Status: succeeded

Final run completed at 2026-07-05 05:18:19 +0800.

Summary:
- Formal fragments-backed ChromatinAssay workflow completed for 15 libraries.
- Cells before QC: 122218.
- Cells after QC: 102419.
- Peaks after assay creation: 1307563.
- Libraries after QC: 15.
- Samples after QC: 45.
- Clusters after QC: 63.
- Transferred labels: 8.
- Label transfer features: 2524.
- Blacklist ratio was not applied because a local hg38 blacklist resource was not available.
- Median prediction score across retained cells: 0.4132543.

Implementation note:
- The full GeneActivity matrix over all genes failed due to memory allocation in an earlier run.
- The final completed run restricted GeneActivity to reference-overlapping transfer genes, using 2781 requested genes and retaining 2524 genes for label transfer.
- This preserves the requested label-transfer workflow while avoiding an infeasible all-gene activity matrix on the local machine.

Primary outputs:
- `tables/fig4a_qc_metrics_by_cell.tsv.gz`
- `tables/fig4a_qc_thresholds.tsv`
- `tables/fig4a_cells_before_after_qc_by_library.tsv`
- `tables/fig4a_snatac_umap_labels.tsv.gz`
- `tables/fig4a_label_counts.tsv`
- `tables/fig4a_prediction_score_summary_by_label.tsv`
- `tables/fig4a_library_label_composition.tsv`
- `figures/fig4a_formal_qc_violin.png`
- `figures/fig4a_formal_umap_label_transfer.png`
- `figures/fig4a_formal_umap_by_library.png`
- `figures/fig4a_formal_umap_by_tissue_group.png`
- `figures/fig4a_formal_umap_prediction_score.png`

Quality-control and interpretation notes:
- The request checker equivalent passed after generating `STATUS.md`, `upload_manifest.tsv`, and `checksums.sha256`.
- Visual QA passed for the five required PNG figures: no title/axis/legend clipping was observed.
- The UMAP-by-library diagnostic shows strong library-associated structure, including several library-dominant regions.
- Several libraries are dominated by one transferred label; this is captured in `tables/fig4a_library_label_composition.tsv`.
- Label-transfer confidence is moderate overall, with median prediction score 0.4132543 and mean prediction score 0.4439141.
- The THBS1+ MDSCs label has only 29 cells and low prediction scores, so it should be interpreted cautiously.

The full processed Seurat object is not included in the GitHub handoff package.
