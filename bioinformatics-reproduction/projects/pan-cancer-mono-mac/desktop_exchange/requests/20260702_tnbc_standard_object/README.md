# TNBC Standard Object Request

## Purpose

The uploaded TNBC full run completed successfully, but the uploaded GitHub files do not include the large Seurat RDS. The PC still has the processed object locally:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\tnbc_full_run\outputs\04_data_processed\single_cell\full_scope_preprocess\tnbc_bpcells_full\tnbc_full_marker_clustering_seurat.rds
```

This request asks the PC to subset that object to the reviewed myeloid clusters and write the standard object needed by the full-scope integration workflow.

## Run

From:

```bash
cd shengxinfenxi/bioinformatics-reproduction/projects/pan-cancer-mono-mac/tnbc_full_run
```

Run:

```bash
STANDARDIZE_INPUT_RDS="outputs/04_data_processed/single_cell/full_scope_preprocess/tnbc_bpcells_full/tnbc_full_marker_clustering_seurat.rds" \
STANDARDIZE_CLUSTER_ANNOTATION="../desktop_exchange/requests/20260702_tnbc_standard_object/tnbc_cluster_annotation_for_standard_object.tsv" \
STANDARDIZE_OUTPUT_RDS="outputs/04_data_processed/single_cell/full_scope_standard_objects/set2_TNBC_standard_myeloid.rds" \
STANDARDIZE_DATASET_ID="set2_TNBC" \
STANDARDIZE_PAPER_LABEL="TNBC" \
STANDARDIZE_CLUSTER_COLUMN="seurat_clusters" \
STANDARDIZE_SOURCE_CELL_COLUMN="colnames" \
STANDARDIZE_PATIENT_COLUMN="patient_from_title" \
STANDARDIZE_ORIGIN_COLUMN="tissue" \
STANDARDIZE_BATCH_COLUMN="sample_id_from_title" \
Rscript ../full_scope_integration_pc/scripts/make_standard_seurat_object_from_review.R
```

Expected retained clusters:

```text
4, 5, 8, 10, 11
```

Expected retained cells from uploaded PC cluster tables:

```text
75151
```

## Upload Back

Create:

```text
desktop_exchange/uploads/20260702_tnbc_standard_object/
```

Upload:

- `STATUS.md`
- `upload_manifest.tsv`
- `checksums.sha256`
- `tables/set2_TNBC_standard_object_summary.tsv`
- the standard RDS if it is below GitHub's file-size limit
- otherwise split the RDS using `split_file_for_github.py` in this request folder and upload all generated `.partNNN` chunks plus the generated chunk manifest

Example split command:

```bash
python3 ../desktop_exchange/requests/20260702_tnbc_standard_object/split_file_for_github.py \
  --input outputs/04_data_processed/single_cell/full_scope_standard_objects/set2_TNBC_standard_myeloid.rds \
  --out-dir ../desktop_exchange/uploads/20260702_tnbc_standard_object/files \
  --chunk-size-mb 90
```

Do not upload raw GEO inputs or the BPCells matrix directory.
