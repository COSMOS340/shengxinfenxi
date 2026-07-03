# PC Request for Mac: Fig4 GSE306459 snATAC

## Status

PC found the request on branch `fig4-snatac-request-local-20260703` at:

```text
bioinformatics-reproduction/projects/pan-cancer-mono-mac/desktop_exchange/requests/20260703_fig4_gse306459_snatac/
```

The workflow cannot be run to completion yet because required inputs are missing on the PC side.

## Required From Mac

Please provide the exact scRNA reference object and label column required by the request README:

```text
FIG4_REFERENCE_RDS=<exact RDS path or uploaded split chunks>
FIG4_REFERENCE_LABEL_COLUMN=<exact metadata column name present in that RDS>
```

The script explicitly requires `FIG4_REFERENCE_RDS` to be a Seurat object and will stop unless `FIG4_REFERENCE_LABEL_COLUMN` is an exact metadata column in that object. The README example uses `fig4_display_label`, but PC cannot assume that column exists without the RDS or a metadata-column listing.

Please also provide the supplementary xlsx files or confirm a working download URL:

```text
cir-24-1255_table_s1_suppst1.xlsx
cir-24-1255_table_s4_suppst4.xlsx
```

PC searched local `I:\` and `H:\` and did not find these files. The request TSV lists PMC paths under `/articles/instance/12865363/bin/`, but direct checks returned HTML rather than an xlsx payload, so PC did not treat those responses as valid Excel files.

## PC Environment Notes

The following R packages are currently missing from the PC R library and will need installation before the full run:

```text
Signac
ensembldb
EnsDb.Hsapiens.v86
BSgenome.Hsapiens.UCSC.hg38
JASPAR2020
TFBSTools
ChIPseeker
TxDb.Hsapiens.UCSC.hg38.knownGene
ggseqlogo
```

PC can install these to `H:\R\library` after the reference/supplementary inputs are available.

## Files Already Synced To PC

The request files were synced locally under:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\requests\20260703_fig4_gse306459_snatac
```

No GEO raw archive or fragment files have been uploaded back to GitHub. PC has not started the heavy GEO download/run because the reference RDS and exact label column are blocking prerequisites.
