# 20260704 Supplementary Fig S6 spatial signatures PC request

## Purpose

Run the complete Supplementary Fig. S6 spatial signature reproduction on the PC.

Original Supplementary Fig. S6 contains:

- A-D: spatial feature plots for 10x Genomics GBM, GSE226997 colorectal cancer, GSE274103 PDAC, and 10x Genomics colorectal cancer.
- E-H: Spearman correlation plots between SPP1+ TAM signature and MDSC signature for the four datasets above.

The Mac currently has only the 10x Genomics colorectal cancer object, so local output is CRC-only.

## Signature source

This request includes:

```text
inputs/fig5_signature_gene_sets_unique_current_scope.tsv
```

The script uses exact `components` values:

```text
SPP1+ TAMs top10
MDSC 19 genes
```

SPP1+ TAM signature is scored from the 10 genes tagged `SPP1+ TAMs top10`. MDSC signature is scored from the 19 genes tagged `MDSC 19 genes`.

## Run

From the repository root on PC:

```powershell
cd I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac
powershell -ExecutionPolicy Bypass -File desktop_exchange\requests\20260704_supp_fig_s6_spatial_signatures\scripts\run_supp_fig_s6_spatial_signatures.ps1
```

The runner defaults to reusing the data folder from the S4/S5 request:

```text
desktop_exchange\requests\20260704_supp_fig_s4_s5_spatial_featureplots\data
```

Override with:

```powershell
$env:S6_DATA_DIR = "I:\path\to\existing\spatial\data"
```

If `GSE226997` is too large for the available disk, set:

```powershell
$env:S6_SKIP_GSE226997 = "1"
```

If `GSE274103` should be skipped:

```powershell
$env:S6_SKIP_GSE274103 = "1"
```

If skipped, upload partial outputs plus the log and selected-source table.

## Data sources

The script uses the same source definitions as the S4/S5 request:

- 10x Genomics GBM: `CytAssist_11mm_FFPE_Human_Glioblastoma`
- 10x Genomics colorectal cancer: `CytAssist_11mm_FFPE_Human_Colorectal_Cancer`
- GSE274103 PDAC: `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE274103&format=file`
- GSE226997 colorectal cancer: `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE226997&format=file`

Do not upload R packages, raw archives, H5 files, image files, or extracted spatial folders into GitHub.

## Expected upload

Upload generated outputs under:

```text
desktop_exchange/uploads/20260704_supp_fig_s6_spatial_signatures/
```

Important files:

```text
outputs/figures/supp_fig_s6_spatial_signatures.png
outputs/figures/supp_fig_s6_spatial_signatures.pdf
outputs/tables/supp_fig_s6_dataset_manifest.tsv
outputs/tables/supp_fig_s6_selected_spatial_sources.tsv
outputs/tables/supp_fig_s6_gene_presence.tsv
outputs/tables/supp_fig_s6_signature_gene_coverage.tsv
outputs/tables/supp_fig_s6_signature_scores.tsv.gz
outputs/tables/supp_fig_s6_spearman_correlations.tsv
outputs/logs/run_supp_fig_s6_spatial_signatures.log
```

If a dataset fails, do not replace it silently. Upload the log, manifest, selected-source table, and completed panels.
