# 20260704 Supplementary Fig S4-S5 Spatial Feature Plots Request

## Purpose

Run the complete spatial transcriptomics feature-plot reproduction for Supplementary Fig. S4 and Supplementary Fig. S5.

The Mac has only the 10x Genomics colorectal cancer Visium data locally. The original S4 and S5 PDFs are visually and textually duplicated except for the figure number, and both show spatial feature plots for:

- 10x Genomics GBM: `EGFR`, `C1QC`, `FOLR2`
- GSE226997 colorectal cancer: `EPCAM`, `C1QC`, `FOLR2`
- GSE274103 PDAC: `EPCAM`, `C1QC`, `PLTP`
- 10x Genomics colorectal cancer: `EPCAM`, `C1QC`, `FOLR2`

## Run

From this request directory:

```bash
Rscript scripts/run_supp_fig_s4_s5_spatial_featureplots.R
```

Do not upload R packages into the repository. Install packages locally on the PC if needed.

## Data Sources

### 10x Genomics GBM

Use the 10x Genomics sample:

`CytAssist_11mm_FFPE_Human_Glioblastoma`

Download URLs verified from the 10x dataset page on 2026-07-04:

- `https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_filtered_feature_bc_matrix.h5`
- `https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_spatial.tar.gz`
- `https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_image.tif`

### 10x Genomics Colorectal Cancer

Use the same 10x sample already used by the Mac Fig. 3 spatial reproduction:

`CytAssist_11mm_FFPE_Human_Colorectal_Cancer`

The script downloads it if it is absent on the PC.

### GSE274103 PDAC

Official GEO download:

- `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE274103&format=file`
- Expected HTTP attachment: `GSE274103_RAW.tar`
- Expected size from GEO header on 2026-07-04: `186603520` bytes

The script extracts all nested tar files, finds Visium-like sample folders, and uses the first sorted sample unless `S4S5_GSE274103_SAMPLE` is set.

### GSE226997 Colorectal Cancer

Official GEO download:

- `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE226997&format=file`
- Expected attachment: `GSE226997_RAW.tar`
- Expected size from GEO page on 2026-07-04: `44193679360` bytes, displayed as `41.2 Gb`

This is too large for the 16 GB Mac local workflow. The PC should run it if space allows. If not, set:

```bash
export S4S5_SKIP_GSE226997=1
```

and upload the partial output plus the recorded skip reason.

## Expected Upload

Upload the generated `outputs/` directory under:

`desktop_exchange/uploads/20260704_supp_fig_s4_s5_spatial_featureplots/`

Important files:

- `outputs/figures/supp_fig_s4_s5_spatial_featureplots.png`
- `outputs/figures/supp_fig_s4_s5_spatial_featureplots.pdf`
- `outputs/tables/supp_fig_s4_s5_dataset_manifest.tsv`
- `outputs/tables/supp_fig_s4_s5_gene_presence.tsv`
- `outputs/tables/supp_fig_s4_s5_selected_spatial_sources.tsv`
- `outputs/logs/run_supp_fig_s4_s5_spatial_featureplots.log`

If a dataset fails, do not silently replace it with another dataset. Upload the logs, the manifest, and any completed panels.
