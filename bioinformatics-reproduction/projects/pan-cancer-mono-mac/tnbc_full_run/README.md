# TNBC Full Run for Desktop

This folder runs the full GSE169246 TNBC RNA analysis for the pan-cancer monocyte/macrophage lineage reproduction.

The raw data are not stored in this GitHub repository. The desktop machine downloads the four public GEO files directly from NCBI GEO, verifies SHA256 checksums, parses sample metadata, and then runs the BPCells/Seurat/Harmony workflow.

## Scope

- GEO accession: `GSE169246`
- RNA matrix: `GSE169246_TNBC_RNA.counts.mtx.gz`
- Matrix dimensions recorded from the source audit: 27,085 features, 489,490 cells, 646,867,807 nonzero entries
- GEO series metadata rows parsed locally: 83 total sample rows, including 78 RNA samples
- Patient scope note: the RNA barcode suffixes expose 22 patient IDs, while the target paper reports TNBC `n = 21`. No exact public file in the current evidence set identifies which one patient should be removed, so this workflow keeps all 22 patient IDs and records the mismatch for later evidence-based filtering.

## Hardware

Recommended:

- RAM: 64 GB or more
- Disk: at least 80 GB free space for downloaded inputs, BPCells disk matrix, Seurat object, marker tables, figures, and temporary files
- CPU: 8 or more cores

## R Dependencies

Install R packages before running:

```r
install.packages(c("Seurat", "harmony", "data.table", "ggplot2"))
install.packages("BPCells")
```

If `BPCells` is not available from the active CRAN mirror, install it using the official BPCells instructions for your R version.

## Run From a Desktop

Clone the repository:

```bash
git clone https://github.com/COSMOS340/shengxinfenxi.git
cd shengxinfenxi/bioinformatics-reproduction/projects/pan-cancer-mono-mac/tnbc_full_run
```

Download and verify GEO inputs:

```bash
python3 scripts/download_tnbc_inputs.py --out inputs
```

Parse GEO sample metadata:

```bash
python3 scripts/parse_tnbc_series_matrix.py
```

Run the full analysis:

```bash
TNBC_THREADS=8 TNBC_NPCS=30 TNBC_RESOLUTION=0.1 \
Rscript scripts/run_tnbc_bpcells_full.R 2>&1 | tee tnbc_full_run.log
```

By default, marker testing uses all cells. If marker testing is too slow, rerun with a per-cluster marker cap:

```bash
TNBC_MARKER_MAX_CELLS_PER_CLUSTER=20000 \
Rscript scripts/run_tnbc_bpcells_full.R 2>&1 | tee tnbc_full_run_marker_limited.log
```

Use the unrestricted marker run first if the desktop hardware allows it.

## Outputs

Outputs are written under `outputs/`:

- `outputs/04_data_processed/single_cell/full_scope_preprocess/tnbc_bpcells_full/GSE169246_TNBC_RNA_counts_bpcells/`
- `outputs/04_data_processed/single_cell/full_scope_preprocess/tnbc_bpcells_full/tnbc_full_marker_clustering_seurat.rds`
- `outputs/07_tables/main/full_scope_preprocess/tnbc_bpcells_full/tnbc_full_run_summary.tsv`
- `outputs/07_tables/main/full_scope_preprocess/tnbc_bpcells_full/tnbc_full_cluster_top10_markers.tsv`
- `outputs/07_tables/main/full_scope_preprocess/tnbc_bpcells_full/tnbc_full_cluster_markers_all.tsv`
- `outputs/07_tables/main/full_scope_preprocess/tnbc_bpcells_full/tnbc_full_umap_coordinates.tsv`
- `outputs/06_figures/main/full_scope_preprocess/tnbc_bpcells_full/tnbc_full_umap_by_cluster.png`
- `outputs/06_figures/main/full_scope_preprocess/tnbc_bpcells_full/tnbc_full_canonical_marker_dotplot.png`

After the run, inspect `tnbc_full_cluster_top10_markers.tsv` and the marker dotplot manually before assigning cell-type names. Do not use automatic labels as final cell-type names.

## Downloaded Files and Checksums

The downloader verifies these exact SHA256 values:

```text
ad6c784ec3f4d965e2ad15e22467c66aa185eceead8808b956e75b33cfb7c76d  GSE169246_TNBC_RNA.counts.mtx.gz
c61d02185e54905a1f7ed658d69fd588db2526a9130951272240338c083a9141  GSE169246_TNBC_RNA.barcode.tsv.gz
83f6d996d1a539ec02161884a3ee29920ab3c56207a4242dd93ff62be2702b60  GSE169246_TNBC_RNA.feature.tsv.gz
80d720ca21d28916d26fc31ec7cfb93fcbe4f9445dc16ea2f0382cf36c8f3c88  GSE169246_series_matrix.txt.gz
```

To recheck after download:

```bash
cd inputs
shasum -a 256 -c checksums.sha256
```

## Failure Handling

If the R run fails, keep these files for debugging:

- `tnbc_full_run.log`
- the last printed step
- the generated `outputs/` directory
- `inputs/checksums.sha256`

Do not delete the BPCells directory unless the failure is clearly during import and a clean import rerun is needed.
