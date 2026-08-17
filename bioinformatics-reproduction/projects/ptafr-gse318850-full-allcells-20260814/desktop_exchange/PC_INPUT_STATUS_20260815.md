# PC input status - 2026-08-15

The PTAFR all-cell runner and manifest were found and synchronized on the PC.

The following exact manifest inputs are absent from all mounted local data drives checked by the PC:

- `04_data_processed/single_cell/GSE318850/GSE318850_RNA_counts_10x.h5`
- `04_data_processed/single_cell/GSE318850/reannotation/GSE318850_full_independent_cell_metadata.rds`
- `04_data_processed/single_cell/GSE318850/qc_doublet/GSE318850_qc_doublet_metadata.rds`

The PC has started a resumable, range-based download of the official GEO source object:

- Source: `GSE318850_seurat_all_cells.qs.gz`
- Exact remote size: `8695855276` bytes
- Local path: `I:\ptafr\05_公共数据库复现_图ABCD\00_data_raw\GSE318850\GSE318850_seurat_all_cells.qs.gz`

Please add one of the following to this branch:

1. Exact external download locations for the three manifest inputs, or
2. The exact preparation scripts and package versions used to generate the three files from the GEO source object.

Do not add the 1.76 GB H5 file directly to GitHub. The PC will continue the official source download and will verify all final inputs against `input_manifest.tsv` before launching the R analysis.
