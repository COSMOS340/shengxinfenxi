# Pan-GI ICI CCC matrix download and lightweight inspection

Completed: 2026-07-09

GitHub handoff directory: bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_matrix_download
Local matrix directory: I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\03_matrices

Summary:
- Manifest rows: 32
- Downloaded: 32
- Failed: 0
- Local matrix bytes: 8538523347
- All matrix files are local_only_large_file in matrix_file_inventory.tsv; no matrix payloads were uploaded to GitHub.

Lightweight inspection:
- GSE189926 gzip text matrices: rows appear to be genes and columns appear to be cells based on first-line structure.
- GSE235863 h5ad.gz files: HDF5 structure inspected after temporary decompression; obs/var dimensions and obs columns recorded.
- GSE236581 metadata/count matrix set: metadata has 975275 rows and 9 columns; no response-related metadata column was detected; MatrixMarket dimensions are 36027 x 975275 with 1310816895 nonzero entries.
- ICB_Zenodo_Liver_Ma RDS: readRDS inspected a Seurat object with 21324 genes x 9451 cells and Cancer_type labels HCC/iCCA.
- OMIX001073 zip files: archive entries and uncompressed byte totals recorded only; matrices were not extracted.

No normalization, clustering, integration, CellChat, LIANA, plotting, or full expression export was run.

