# Figure 4 GSE306459 snATAC Request

## Purpose

Run the heavy GSE306459 CRC snATAC-seq part of Fig. 4 on the PC.

The Mac has already generated current-scope Fig. 4C-F from Supplementary Tables S2/S4 and COAD-READ scRNA expression. The remaining panels need the real snATAC object and fragment-level outputs:

- Fig. 4A: snATAC UMAP with transferred labels.
- Fig. 4B: adjacent normal versus tumor proportions by transferred label.
- Fig. 4C: full ChIPSeeker annotation pie from differential peaks.
- Fig. 4G: coverage tracks for `C1QB`, `CHID1`, `IL1B`, and `SLC25A37`.
- Fig. 4H: chromVAR motif accessibility UMAPs for `TFEC`, `TFE3`, `BHLHE41`, `BACH1`, `Nfe2l2`, and `FOSL2`.
- Fig. 4I: motif logos for the same TF set.

## Exact Inputs

Required GEO file list:

```text
desktop_exchange/requests/20260703_fig4_gse306459_snatac/inputs/gse306459_geo_filelist.tsv
```

NCBI GEO exposes only this archive in the public supplementary directory:

```text
GSE306459_RAW.tar
```

The individual `raw_peak_bc_matrix.h5` and `fragments.tsv.gz` entries in `filelist.txt` are archive members, not standalone FTP files. The script downloads or reuses `GSE306459_RAW.tar`, then extracts only the files needed for the selected run scope.

Required supplementary tables:

```text
cir-24-1255_table_s1_suppst1.xlsx
cir-24-1255_table_s4_suppst4.xlsx
```

Table S1 is required because the script uses the exact `barcode` sheet columns:

```text
Cell barcode
Library
Tissue
Sample
```

Required scRNA reference for label transfer:

```text
FIG4_REFERENCE_RDS
FIG4_REFERENCE_LABEL_COLUMN
```

`FIG4_REFERENCE_RDS` must be a Seurat object. `FIG4_REFERENCE_LABEL_COLUMN` must be the exact metadata column containing the labels to transfer. The script will stop and print observed metadata columns if this column is absent.

## Run

From this repository root:

```bash
cd bioinformatics-reproduction/projects/pan-cancer-mono-mac
```

Run:

```bash
FIG4_RUN_SCOPE="fig4c" \
FIG4_GEO_DIR="fig4_snatac/raw/GSE306459" \
FIG4_SUPP_DIR="fig4_snatac/inputs/supplementary" \
FIG4_OUTPUT_DIR="fig4_snatac/outputs" \
FIG4_REFERENCE_RDS="/path/to/coad_read_reference.rds" \
FIG4_REFERENCE_LABEL_COLUMN="fig4_display_label" \
Rscript desktop_exchange/requests/20260703_fig4_gse306459_snatac/scripts/run_fig4_gse306459_snatac.R
```

`FIG4_RUN_SCOPE=fig4c` is the current requested run. It extracts the H5 peak matrices, runs label transfer, finds differential peaks, and generates the complete Fig. 4C ChIPSeeker annotation pie. It does not extract fragment files.

Use `FIG4_RUN_SCOPE=all` only when the PC should also generate coverage tracks and motif accessibility panels. That mode extracts fragment files and is much heavier.

The output directory can be changed, but keep it outside `desktop_exchange/uploads` while the job runs.

## Upload Back

Create:

```text
desktop_exchange/uploads/20260703_fig4_gse306459_snatac/
```

Upload these files:

- `STATUS.md`
- `upload_manifest.tsv`
- `checksums.sha256`
- `logs/fig4_gse306459_snatac.log`
- `tables/fig4a_snatac_umap_labels.tsv.gz`
- `tables/fig4b_label_tissue_proportions.tsv`
- `tables/fig4b_wilcoxon_by_label.tsv`
- `tables/fig4c_chipseeker_annotation_counts.tsv`
- `tables/fig4c_chipseeker_peak_annotation_full.tsv.gz`
- `tables/fig4_differential_peaks_by_label.tsv.gz`
- `tables/fig4_gse306459_run_summary.tsv`
- `tables/fig4e_chromvar_motif_activity_by_label.tsv`
- `tables/fig4g_track_region_summary.tsv`
- `figures/fig4a_snatac_umap_label_transfer.png`
- `figures/fig4b_tissue_proportion_boxplots.png`
- `figures/fig4c_chipseeker_annotation_pie.png`
- `figures/fig4g_gene_coverage_tracks.png`
- `figures/fig4h_motif_accessibility_umaps.png`
- `figures/fig4i_motif_logos.png`
- `files/fig4_gse306459_snatac_processed.rds.sha256`
- split RDS chunks if the processed object is needed for Mac-side review

For `FIG4_RUN_SCOPE=fig4c`, the motif/coverage outputs are not expected. Upload the Fig. 4A/B/C figures and C-related tables listed above.

Example split command:

```bash
python3 desktop_exchange/requests/20260703_fig4_gse306459_snatac/scripts/split_file_for_github.py \
  --input fig4_snatac/outputs/files/fig4_gse306459_snatac_processed.rds \
  --out-dir desktop_exchange/uploads/20260703_fig4_gse306459_snatac/files \
  --chunk-size-mb 90
```

Do not upload the 87.55 GB GEO raw archive or the fragment files.
