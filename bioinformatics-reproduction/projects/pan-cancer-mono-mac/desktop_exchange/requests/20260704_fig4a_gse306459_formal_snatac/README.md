# PC Request: Formal Fig. 4A GSE306459 snATAC Rerun

Date: 2026-07-04

This request replaces the previous fast-merge Fig. 4A diagnostic output. Do not polish or reuse the fast-merge UMAP as a final Fig. 4A reproduction.

## Why this rerun is needed

The uploaded fast-merge Fig. 4A is not equivalent to the paper panel.

Evidence from Mac-side review:

- The UMAP is strongly structured by `Library` and `tissue_group`, not by transferred myeloid cell type.
- Label transfer confidence is low for key paper labels.
- `THBS1+ MDSCs` has only 14 transferred cells in the fast-merge output.
- Several original Fig. 4A labels are absent from the transferred labels.
- The fast-merge script used H5 peak matrices only and skipped fragments-based snATAC QC.

## Required rerun scope

Run a formal GSE306459 snATAC workflow for Fig. 4A using fragments and paper-style QC before LSI/UMAP/label transfer.

This is a Fig. 4A-focused rerun. Fig. 4C can be regenerated from the same filtered object if feasible, but the priority is a correct Fig. 4A UMAP and label-transfer table.

## Required inputs

Use the existing inputs already provided in this exchange:

- Request file list and original script:
  - `desktop_exchange/requests/20260703_fig4_gse306459_snatac/inputs/gse306459_geo_filelist.tsv`
  - `desktop_exchange/requests/20260703_fig4_gse306459_snatac/scripts/run_fig4_gse306459_snatac.R`
- Mac-provided reference and supplementary files:
  - `desktop_exchange/uploads/20260703_fig4_gse306459_snatac/mac_provided_inputs/files/fig4_coad_read_reference_with_fig4_display_label.rds`
  - `desktop_exchange/uploads/20260703_fig4_gse306459_snatac/mac_provided_inputs/supplementary/cir-24-1255_table_s1_suppst1.xlsx`
  - `desktop_exchange/uploads/20260703_fig4_gse306459_snatac/mac_provided_inputs/supplementary/cir-24-1255_table_s4_suppst4.xlsx`

If the RDS is still split into `.part000` and `.part001`, reconstruct it using the command in:

- `desktop_exchange/uploads/20260703_fig4_gse306459_snatac/MAC_RESPONSE.md`

R packages do not need to be uploaded by Mac. Install missing R/Bioconductor packages locally on PC.

## Required method changes versus fast-merge

Do not use the fast-merge workflow as final Fig. 4A.

Required changes:

1. Extract both `raw_peak_bc_matrix.h5` and `fragments.tsv.gz` members from `GSE306459_RAW.tar`.
2. Build a fragments-backed `ChromatinAssay`.
3. Compute and export snATAC QC metrics before filtering:
   - `nCount_peaks`
   - `nFeature_peaks`
   - `TSS.enrichment`
   - `nucleosome_signal`
   - `pct_reads_in_peaks`
   - `blacklist_ratio`, if the hg38 blacklist resource is available locally
4. Apply a documented QC filter before LSI/UMAP. If exact thresholds are not in the paper, choose thresholds from the metric distributions, but export the thresholds and before/after counts.
5. Run `RunTFIDF`, `FindTopFeatures`, `RunSVD`, `RunUMAP`, `FindNeighbors`, and `FindClusters` only after filtering.
6. Perform label transfer using `FIG4_REFERENCE_LABEL_COLUMN=fig4_display_label`.
7. Export diagnostics that show whether the final UMAP is still library-driven:
   - UMAP colored by transferred label
   - UMAP colored by `Library`
   - UMAP colored by `tissue_group`
   - UMAP colored by `prediction.score.max`
   - label counts and prediction-score summary
   - Library-by-label composition table

## Required environment variables

Use paths adapted to the PC checkout, for example:

```powershell
$env:FIG4_GEO_DIR = "fig4_snatac/raw/GSE306459"
$env:FIG4_SUPP_DIR = "desktop_exchange/uploads/20260703_fig4_gse306459_snatac/mac_provided_inputs/supplementary"
$env:FIG4_OUTPUT_DIR = "desktop_exchange/uploads/20260704_fig4a_gse306459_formal_snatac/pc_formal_fig4a_results"
$env:FIG4_REFERENCE_RDS = "desktop_exchange/uploads/20260703_fig4_gse306459_snatac/mac_provided_inputs/files/fig4_coad_read_reference_with_fig4_display_label.rds"
$env:FIG4_REFERENCE_LABEL_COLUMN = "fig4_display_label"
$env:FIG4_RUN_SCOPE = "fig4a_formal"
```

## Required upload directory

Upload results to:

`desktop_exchange/uploads/20260704_fig4a_gse306459_formal_snatac/pc_formal_fig4a_results/`

## Required output files

Create these files:

- `STATUS.md`
- `upload_manifest.tsv`
- `checksums.sha256`
- `logs/fig4a_formal_snatac.log`
- `tables/fig4a_formal_run_summary.tsv`
- `tables/fig4a_qc_metrics_by_cell.tsv.gz`
- `tables/fig4a_qc_thresholds.tsv`
- `tables/fig4a_cells_before_after_qc_by_library.tsv`
- `tables/fig4a_snatac_umap_labels.tsv.gz`
- `tables/fig4a_label_counts.tsv`
- `tables/fig4a_prediction_score_summary_by_label.tsv`
- `tables/fig4a_library_label_composition.tsv`
- `figures/fig4a_formal_umap_label_transfer.png`
- `figures/fig4a_formal_umap_by_library.png`
- `figures/fig4a_formal_umap_by_tissue_group.png`
- `figures/fig4a_formal_umap_prediction_score.png`
- `figures/fig4a_formal_qc_violin.png`

Do not upload the full processed Seurat object unless it is small enough for GitHub. If the object is large, save only:

- `files/fig4a_formal_processed_object_not_uploaded.txt`

## Acceptance criteria

The rerun is useful only if:

- The filtered UMAP is not dominated by single-library islands.
- Key original Fig. 4A labels are recovered at plausible counts.
- `THBS1+ MDSCs`, `SPP1+ TAMs`, and C1QC+ TAM labels are not near-zero artifacts.
- The prediction-score summaries support the transferred labels.

If these criteria fail, upload the diagnostics anyway and state the failure explicitly in `STATUS.md`.
