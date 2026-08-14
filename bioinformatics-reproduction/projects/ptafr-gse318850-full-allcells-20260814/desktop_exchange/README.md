# PTAFR GSE318850 all-cell PC handoff

## Scope

This handoff rebuilds the GSE318850 singlet atlas with all 649,163 retained singlets. It does not use a sketch, leverage-score sampling, uniform sampling, or SingleR. Cell-type annotation is performed manually after review of the marker table generated from all cells.

The source H5 file is an author-provided processed count matrix. The supplied object is already cell-filtered; it is not FASTQ or a Cell Ranger raw-feature matrix. Sample-wise `scDblFinder` results were produced previously for all 718,118 input cells, leaving 649,163 singlets.

## Git branch

Repository: `COSMOS340/shengxinfenxi`

Branch: `agent/ptafr-full-allcells-20260814`

On the PC:

```bash
cd /mnt/i
git clone git@github.com:COSMOS340/shengxinfenxi.git
cd shengxinfenxi
git fetch origin agent/ptafr-full-allcells-20260814
git switch --track -c agent/ptafr-full-allcells-20260814 origin/agent/ptafr-full-allcells-20260814
```

If the repository already exists, use its existing checkout and run only `git fetch` and `git switch`.

## Required project files

The default project root is:

`/mnt/i/ptafr/05_公共数据库复现_图ABCD`

The three required inputs and their exact hashes are recorded in `input_manifest.tsv`. The 1.76 GB H5 count matrix is intentionally not stored in GitHub.

The launcher verifies the byte size and SHA-256 hash of every input before R starts.

## Launch

From the repository checkout:

```bash
cd /mnt/i/shengxinfenxi/bioinformatics-reproduction/projects/ptafr-gse318850-full-allcells-20260814/desktop_exchange
chmod +x run_on_pc.sh
nohup bash run_on_pc.sh > /mnt/i/ptafr/05_公共数据库复现_图ABCD/99_logs/GSE318850_PC_launcher.log 2>&1 &
```

If the project root differs, set it explicitly without changing the script:

```bash
PTAFR_PROJECT_ROOT=/mnt/i/exact/project/path nohup bash run_on_pc.sh > launcher.log 2>&1 &
```

Monitor progress:

```bash
tail -f /mnt/i/ptafr/05_公共数据库复现_图ABCD/99_logs/GSE318850_PC_launcher.log
```

## Analysis sequence

1. Validate the exact cell order across the H5 counts and two RDS metadata files.
2. Retain exactly 649,163 sample-wise `scDblFinder` singlets.
3. Normalize all singlets and select 2,500 variable genes using all cells.
4. Run BPCells disk-backed centering, scaling, and 40-component truncated SVD on all cells.
5. Save an all-cell PCA checkpoint.
6. Run Harmony across the exact `orig.ident` library field using all cells.
7. Save an all-cell Harmony checkpoint.
8. Build a 30-neighbor index and SNN graph from all cells.
9. Scan clustering resolutions 0.2 to 0.6 and use resolution 0.4 for the first marker audit.
10. Run UMAP for all 649,163 cells.
11. Calculate all-cell cluster markers with BPCells Wilcoxon tests and expression proportions.
12. Manually annotate clusters from the current marker evidence. Neutrophils remain one broad `Neutrophils` class.

## Expected checkpoints and outputs

Output directory:

`/mnt/i/ptafr/05_公共数据库复现_图ABCD/04_data_processed/single_cell/GSE318850/reanalysis_full`

Key files:

- `GSE318850_full_pca_checkpoint.rds`
- `GSE318850_full_harmony_checkpoint.rds`
- `GSE318850_full_harmony_stage1.rds`
- `GSE318850_full_harmony_umap_coordinates.tsv.gz`
- `GSE318850_full_harmony_cluster_markers.tsv.gz`
- `GSE318850_full_harmony_cluster_top50_markers.tsv`
- `GSE318850_full_resolution_summary.tsv`
- `GSE318850_full_library_centroid_separation_summary.tsv`
- `PC_resource_snapshot_*.json`
- `PC_R_sessionInfo_*.txt`

Do not interpret cluster numbers as cell types. Manual annotation begins only after the all-cell marker tables are complete.
