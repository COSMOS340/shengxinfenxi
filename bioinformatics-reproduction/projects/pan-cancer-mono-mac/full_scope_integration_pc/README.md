# Full-Scope Integration PC Run

This folder is the PC-side workflow for the complete 17-dataset monocyte/macrophage reproduction.

Raw data are not stored in this GitHub repository. This workflow expects one reviewed, standardized Seurat RDS input per paper dataset, then runs the full Harmony integration and exports the Fig. 1/Fig. 2 upstream objects.

## Scope

Paper target:

- 17 datasets
- 131,249 integrated myeloid cells in Fig. 1
- 5 broad Fig. 1 cell types: macrophages, monocytes, classical DCs, plasmacytoid DCs, mast cells
- downstream monocyte/macrophage subclustering for Fig. 2

Current unresolved evidence blockers:

- `set2_TNBC`: run `../tnbc_full_run/` first and upload marker-review outputs.
- `set2_GBM`: exact 3 GBM sample IDs are not identified in the public UCSC metadata inspected so far.
- `set2_healthy_PBMC`: exact 15 healthy PBMC donor IDs are not identified in the AIDA/CELLxGENE metadata inspected so far.

Do not choose GBM or PBMC IDs by matching cell-count bars. Use only explicit source evidence or user-provided evidence.

## Required Standard Input Objects

Each dataset must be supplied as a Seurat RDS object with:

- RNA assay present.
- normalized expression in the RNA `data` layer or slot.
- metadata columns with exact names:
  - `paper_dataset_id`
  - `paper_label`
  - `patient_or_subject`
  - `origin_or_tissue`
  - `source_cell_id`
  - `manual_annotation`
  - `manual_broad_lineage`
  - `batch_id`
- cells already filtered by dataset-level marker review and contamination removal.

For exact Fig. 1 reproduction, `manual_broad_lineage` should use these exact values:

```text
Macrophages
Monocytes
Classical DCs
Plasmacytoid DCs
Mast cells
```

Use `inputs/full_scope_standard_object_manifest_template.tsv` as the manifest format. Do not rename manifest columns.

`batch_id` is the Harmony batch field. Use the exact sample/library field used for dataset-level correction when it is available; otherwise document the evidence and limitation in the upload `STATUS.md`.

## Make a Standard Object From a Reviewed Cluster Run

Use `inputs/cluster_annotation_template.tsv` as the review table format. After marker review, fill one row per cluster and run:

```bash
STANDARDIZE_INPUT_RDS=/path/to/dataset_marker_clustering_seurat.rds \
STANDARDIZE_CLUSTER_ANNOTATION=/path/to/cluster_annotation.tsv \
STANDARDIZE_OUTPUT_RDS=/path/to/set2_TNBC_standard_myeloid.rds \
STANDARDIZE_DATASET_ID=set2_TNBC \
STANDARDIZE_PAPER_LABEL=TNBC \
STANDARDIZE_CLUSTER_COLUMN=seurat_clusters \
STANDARDIZE_SOURCE_CELL_COLUMN=cell_id \
STANDARDIZE_PATIENT_COLUMN=patient_from_title \
STANDARDIZE_ORIGIN_COLUMN=tissue \
STANDARDIZE_BATCH_COLUMN=sample_id_from_title \
Rscript scripts/make_standard_seurat_object_from_review.R
```

For objects whose cell names are already the source cell IDs, set `STANDARDIZE_SOURCE_CELL_COLUMN=colnames`.

## Run

From this folder:

```bash
Rscript scripts/run_full_scope_seurat_harmony.R 2>&1 | tee full_scope_integration.log
```

Recommended PC settings:

```bash
FULL_THREADS=8 FULL_NPCS=30 FULL_RESOLUTION=1.0 FULL_RUN_MARKERS=TRUE \
Rscript scripts/run_full_scope_seurat_harmony.R 2>&1 | tee full_scope_integration.log
```

If marker testing is too slow:

```bash
FULL_MARKER_MAX_CELLS_PER_CLUSTER=30000 \
Rscript scripts/run_full_scope_seurat_harmony.R 2>&1 | tee full_scope_integration_marker_limited.log
```

Use the unrestricted marker run first if hardware allows it.

## Outputs

Outputs are written under `outputs/` by default:

- `outputs/04_data_processed/single_cell/full_scope_integration/full_scope_harmony_integrated_seurat.rds`
- `outputs/07_tables/main/full_scope_integration/full_scope_integration_run_summary.tsv`
- `outputs/07_tables/main/full_scope_integration/full_scope_umap_coordinates.tsv`
- `outputs/07_tables/main/full_scope_integration/full_scope_cluster_markers_all.tsv`
- `outputs/07_tables/main/full_scope_integration/full_scope_cluster_top10_markers.tsv`
- `outputs/06_figures/main/full_scope_integration/full_scope_umap_by_manual_broad_lineage.png`
- `outputs/06_figures/main/full_scope_integration/full_scope_umap_by_dataset.png`
- `outputs/06_figures/main/full_scope_integration/full_scope_umap_by_origin.png`
- `outputs/06_figures/main/full_scope_integration/full_scope_canonical_marker_dotplot.png`

## Upload Results for Codex Review

Upload review outputs to:

```text
bioinformatics-reproduction/projects/pan-cancer-mono-mac/desktop_exchange/uploads/
```

Create a dated folder such as:

```text
desktop_exchange/uploads/20260702_full_scope_integration/
```

Include:

- `STATUS.md`
- `upload_manifest.tsv`
- `checksums.sha256`
- `full_scope_integration.log`
- run summary table
- UMAP coordinate table, if file size allows
- marker tables, if file size allows
- generated PNG/PDF figures

Do not upload raw GEO inputs, CELLxGENE H5AD files, BPCells caches, or large temporary files.
