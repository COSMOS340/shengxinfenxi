# GSE225857 Nonimmune and Spatial Mapping Request

Date: 2026-07-08

## Purpose

Generate review-grade evidence from `GSE225857`, the CRLM dataset with author-annotated nonimmune single-cell data and spatial transcriptomics.

This step should connect the author-defined nonimmune fibroblast clusters with spatial spot-level CAF/ECM programs. It should not recluster the nonimmune single-cell object, should not create final cell labels, and should not upload large matrices or RDS objects.

## Inputs

Use the PC-local raw tar from the core GEO download step:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_core_geo_download/raw_core_geo/GSE225857/GSE225857_RAW.tar`

Expected raw members include:

- `GSM7058755_non_immune_counts.txt.gz`
- `GSM7058755_non_immune_meta.txt.gz`
- `GSM7058756_C1.matrix.mtx.gz`
- `GSM7058756_C1.features.tsv.gz`
- `GSM7058756_C1.barcodes.tsv.gz`
- `GSM7058756_C1_tissue_positions_list.csv.gz`
- `GSM7058757_C2.matrix.mtx.gz`
- `GSM7058757_C2.features.tsv.gz`
- `GSM7058757_C2.barcodes.tsv.gz`
- `GSM7058757_C2_tissue_positions_list.csv.gz`
- `GSM7058758_C3.matrix.mtx.gz`
- `GSM7058758_C3.features.tsv.gz`
- `GSM7058758_C3.barcodes.tsv.gz`
- `GSM7058758_C3_tissue_positions_list.csv.gz`
- `GSM7058759_C4.matrix.mtx.gz`
- `GSM7058759_C4.features.tsv.gz`
- `GSM7058759_C4.barcodes.tsv.gz`
- `GSM7058759_C4_tissue_positions_list.csv.gz`
- `GSM7058760_L1.matrix.mtx.gz`
- `GSM7058760_L1.features.tsv.gz`
- `GSM7058760_L1.barcodes.tsv.gz`
- `GSM7058760_L1_tissue_positions_list.csv.gz`
- `GSM7058761_L2.matrix.mtx.gz`
- `GSM7058761_L2.features.tsv.gz`
- `GSM7058761_L2.barcodes.tsv.gz`
- `GSM7058761_L2_tissue_positions_list.csv.gz`

If any required member is missing, stop with an explicit error.

## Command

Run from the repository root:

```bash
Rscript bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_gse225857_nonimmune_spatial_mapping/scripts/pc_gse225857_nonimmune_spatial_mapping.R
```

The script uses:

- `Matrix`
- `data.table`
- `ggplot2`

The script can install missing CRAN packages on the PC if needed.

## Analysis Scope

1. Read `GSM7058755_non_immune_meta.txt.gz`.
2. Preserve exact author metadata fields including `cluster`, `organs`, `patients`, `patients_organ`, `predicted.doublet`, and `doublet`.
3. Summarize author nonimmune clusters, with explicit focus on these fibroblast clusters:
   - `F01_fibroblast_PRELP`
   - `F02_fibrblast_MCAM`
   - `F03_fibroblast_CXCL14`
   - `F04_fibroblast_C3`
   - `F05_fibroblast_COCH`
   - `F06_cycling_MKI67`
4. Stream only the requested marker genes from `GSM7058755_non_immune_counts.txt.gz`; do not read the full count matrix into memory as a dense object.
5. Compute nonimmune marker expression summaries by author cluster.
6. Parse six spatial samples: `C1`, `C2`, `C3`, `C4`, `L1`, and `L2`.
7. Normalize spatial counts per spot with library-size scaling and log transform.
8. Score predefined CAF/ECM, fibroblast subtype, endothelial, pericyte, epithelial, and cycling programs per spot.
9. Produce review plots for spatial score distributions and spatial coordinates.

## Required Marker Programs

Use the marker-program table embedded in the script. At minimum it should include:

- `fibroblast_ecm`: `COL1A1`, `COL1A2`, `COL3A1`, `DCN`, `LUM`, `PRELP`
- `caf_activation`: `FAP`, `ACTA2`, `POSTN`, `THBS2`, `MMP2`, `MMP11`
- `f01_prelp`: `PRELP`, `DCN`, `LUM`, `COL1A1`, `COL1A2`
- `f02_mcam`: `MCAM`, `RGS5`, `PDGFRB`, `CSPG4`, `ACTA2`, `TAGLN`
- `f03_cxcl14`: `CXCL14`, `CFD`, `C3`, `CXCL12`, `SFRP1`
- `f04_c3`: `C3`, `CFD`, `CXCL12`, `CCL2`
- `f05_coch`: `COCH`, `COL14A1`, `PI16`
- `f06_cycling`: `MKI67`, `TOP2A`, `CENPF`
- `endothelial`: `PECAM1`, `VWF`, `KDR`
- `epithelial_tumor`: `EPCAM`, `KRT8`, `KRT18`, `KRT19`

Do not invent additional sample labels. Use the exact sample IDs and file prefixes from the raw tar.

## Expected GitHub Upload Directory

Upload only small tables, PNG figures, and status files to:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/uploads/20260708_gse225857_nonimmune_spatial_mapping/`

Do not commit raw files, dense count matrices, spatial images, or RDS objects.

## Expected Uploaded Files

- `STATUS.md`
- `gse225857_nonimmune_cluster_counts.tsv`
- `gse225857_nonimmune_fibroblast_cluster_counts_by_organ.tsv`
- `gse225857_nonimmune_marker_expression_by_cluster.tsv`
- `gse225857_nonimmune_marker_gene_coverage.tsv`
- `gse225857_spatial_sample_summary.tsv`
- `gse225857_spatial_marker_gene_coverage.tsv`
- `gse225857_spatial_program_score_summary.tsv`
- `gse225857_spatial_program_scores_wide.tsv`
- `gse225857_nonimmune_marker_dotplot.png`
- `gse225857_spatial_program_score_heatmap.png`
- `gse225857_spatial_fibroblast_ecm_score.png`
- `gse225857_spatial_f02_mcam_score.png`
- `gse225857_spatial_f03_cxcl14_score.png`
- `gse225857_spatial_epithelial_tumor_score.png`

## Interpretation Boundary

This request creates evidence for CAF/ECM program review and spatial localization. It does not finalize cell labels, does not infer histology regions, does not perform deconvolution, and does not claim patient-level clinical associations.
