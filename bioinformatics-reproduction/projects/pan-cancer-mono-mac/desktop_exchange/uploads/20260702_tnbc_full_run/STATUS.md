# Desktop Upload Status

## Run Identity

- Upload folder: uploads/20260702_tnbc_full_run
- Dataset: GSE169246 TNBC RNA
- Workflow: 	nbc_full_run BPCells/Seurat/Harmony full run
- Desktop machine: Windows desktop
- Operator: Codex desktop worker
- Start time: 2026-07-02, Asia/Shanghai
- End time: 2026-07-02, Asia/Shanghai
- Git commit used: $commit

## Commands

``bash
python3 scripts/download_tnbc_inputs.py --out inputs
python3 scripts/parse_tnbc_series_matrix.py
TNBC_THREADS=8 TNBC_NPCS=30 TNBC_RESOLUTION=0.1 Rscript scripts/run_tnbc_bpcells_full.R 2>&1 | tee tnbc_full_run.log
``

## Environment

- OS: Microsoft Windows 11 Professional, version 10.0.22631
- CPU: AMD Ryzen 5 5600 6-Core Processor, 6 cores / 12 logical processors
- RAM: 31.9 GiB visible memory
- Free disk before run: not recorded before download/run; H: had insufficient space for Rtools initially, then was cleared before Rtools install
- R version: R version 4.4.2 (2024-10-31 ucrt)
- Python version: Python 3.12.7
- Key R packages and versions: Seurat 5.4.0; harmony 1.2.4; data.table 1.18.2.1; ggplot2 4.0.2; BPCells 0.3.1

## Inputs

- Input source: public GEO files downloaded by scripts/download_tnbc_inputs.py
- Checksum verification result: passed; see original inputs/checksums.sha256 in desktop run folder
- Any missing input: none detected

## Outputs

- Run status: completed; Rscript exit code 0; log ended with Done
- Output root: I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\tnbc_full_run\outputs
- Marker test mode: unrestricted marker test; TNBC_MARKER_MAX_CELLS_PER_CLUSTER not set, summary value 0
- Number of cells: 489,490 pre-filter; 489,490 post-filter
- Number of genes: 27,085 source features; 3,000 variable features
- Number of clusters: 13 final clusters
- Known warnings: R packages built under R 4.4.3 warnings; Harmony Quick-TRANSfer maximum-step warnings; UMAP method change message; dense matrix warning during marker calculation
- Failed steps: initial dependency setup failed until harmony, SingleCellExperiment, Rtools44, and BPCells were installed; final workflow run completed

## Manual Notes

- Visual issues noticed on desktop: key PNG figures visually inspected by Codex before upload; no obvious axis clipping, legend overflow, or severe text overlap detected. Group/tissue UMAPs are visually dense because of 489,490 cells and should be interpreted with tables.
- Biological concerns: cluster identities are not manually assigned; Mac should review top markers and canonical marker dotplot before annotation
- Files too large for GitHub: tnbc_full_umap_coordinates.tsv (~61.1 MB), tnbc_full_qc_before_filter.tsv (~43.7 MB), tnbc_full_marker_clustering_seurat.rds (~400.7 MB), BPCells matrix directory, raw GEO inputs
- External paths or storage locations: large outputs remain under I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\tnbc_full_run\outputs
- Mac request check before upload: checked desktop_exchange on GitHub before upload; no Mac-side request/inbox/status file was present
