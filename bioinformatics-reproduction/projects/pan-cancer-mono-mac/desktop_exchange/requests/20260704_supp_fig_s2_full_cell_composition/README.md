# 20260704 Supplementary Fig S2 Full Cell Composition Request

## Purpose

Run the complete Supplementary Fig. S2 reproduction on the PC.

Supplementary Fig. S2 is not a macrophage-only plot. It shows global all-cell UMAPs and paired pie charts for:

- `LUAD`
- `LIHC-CHOL`
- `Healthy lung`
- `Healthy liver`
- `NSCLC-GSE148071`
- `HCC-GSE156625`

The Mac completed a current-scope version for `LIHC-CHOL`, `Healthy lung`, and `Healthy liver`, plus a LUAD composition-only pie. The Mac could not compute the LUAD all-cell UMAP because reading the full GSE131907 all-cell matrix exceeded the 16 GB R memory limit. `NSCLC-GSE148071` and `HCC-GSE156625` were not present locally.

## Run

From this request directory:

```bash
Rscript scripts/run_supp_fig_s2_full_cell_composition_pc.R
```

Do not upload R packages into the repository. Install packages locally on the PC if needed.

The script is a structured runner and audit scaffold. It writes data-source manifests, logs, and the known-file-structure dataset outputs. For `NSCLC-GSE148071` and `HCC-GSE156625`, inspect the downloaded files and use exact observed metadata or object columns. Do not guess column names.

## Official Data Sources

All URLs below are official NCBI GEO supplementary file URLs checked on 2026-07-04.

### LUAD: GSE131907

Base:

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/`

Files:

- `GSE131907_Lung_Cancer_cell_annotation.txt.gz`
- `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`
- `GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz`
- `GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz`
- `GSE131907_Lung_Cancer_Feature_Summary.xlsx`

Mac current-scope composition used `Sample_Origin %in% c("tLung", "nLung")` and mapped `T lymphocytes` plus `NK cells` to `T cells`, because the original LUAD panel does not display an NK label and the resulting major fractions match the original pie closely.

### LIHC-CHOL: GSE125449

Base:

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/`

Files:

- `GSE125449_Set1_matrix.mtx.gz`
- `GSE125449_Set1_genes.tsv.gz`
- `GSE125449_Set1_barcodes.tsv.gz`
- `GSE125449_Set1_samples.txt.gz`
- `GSE125449_Set2_matrix.mtx.gz`
- `GSE125449_Set2_genes.tsv.gz`
- `GSE125449_Set2_barcodes.tsv.gz`
- `GSE125449_Set2_samples.txt.gz`

Known annotation column in the samples files: `Type`.

### Healthy lung: GSE122960

Base:

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122960/suppl/`

Archive:

- `GSE122960_RAW.tar`

Use the eight donor filtered H5 files listed by GEO:

- `GSM3489182_Donor_01_filtered_gene_bc_matrices_h5.h5`
- `GSM3489185_Donor_02_filtered_gene_bc_matrices_h5.h5`
- `GSM3489187_Donor_03_filtered_gene_bc_matrices_h5.h5`
- `GSM3489189_Donor_04_filtered_gene_bc_matrices_h5.h5`
- `GSM3489191_Donor_05_filtered_gene_bc_matrices_h5.h5`
- `GSM3489193_Donor_06_filtered_gene_bc_matrices_h5.h5`
- `GSM3489195_Donor_07_filtered_gene_bc_matrices_h5.h5`
- `GSM3489197_Donor_08_filtered_gene_bc_matrices_h5.h5`

If no original all-cell annotation is available in the downloaded files, cluster the donor cells and upload the cluster top-marker table and final global cell type mapping.

### Healthy liver: GSE115469

Base:

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/suppl/`

Files:

- `GSE115469_CellClusterType.txt.gz`
- `GSE115469_Data.csv.gz`

Known annotation columns:

- `CellName`
- `Sample`
- `Cell#`
- `Cluster#`
- `CellType`

### NSCLC-GSE148071

Base:

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/`

Archive:

- `GSE148071_RAW.tar`

The official file list shows per-patient processed expression text files such as `GSM4453576_P1_exp.txt.gz` through `GSM4453617_P42_exp.txt.gz`. GEO states that processed data are provided as supplementary files and raw data are not provided for this record.

PC task:

1. Extract `GSE148071_RAW.tar`.
2. Inspect the exact column and row structure of every `*_exp.txt.gz` file.
3. Build the all-cell expression object.
4. If a cell type annotation is present in the files, use the exact observed field name and record it.
5. If no cell type annotation is present, cluster and assign global cell types from marker evidence. Upload the cluster marker table and the exact final mapping table.

### HCC-GSE156625

Base:

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156625/suppl/`

Human HCC files:

- `GSE156625_HCCmatrix.mtx.gz`
- `GSE156625_HCCbarcodes.tsv.gz`
- `GSE156625_HCCgenes.tsv.gz`
- `GSE156625_HCCscanpyobj.h5ad.gz`

Additional human combined files:

- `GSE156625_HCCF_integrated.RData.gz`
- `GSE156625_HCCFmatrix.mtx.gz`
- `GSE156625_HCCFbarcodes.tsv.gz`
- `GSE156625_HCCFgenes.tsv.gz`
- `GSE156625_HCCFscanpyobj.h5ad.gz`

PC task:

1. Prefer `GSE156625_HCCscanpyobj.h5ad.gz` if it contains UMAP and cell type annotations.
2. Inspect and record exact `obs` column names from the h5ad object.
3. If the HCC h5ad object lacks a global cell type annotation, use the matrix files to cluster and annotate from marker evidence.
4. Do not use mouse files for the human HCC panel.

## Expected Output

Upload under:

`desktop_exchange/uploads/20260704_supp_fig_s2_full_cell_composition/`

Required files:

- `outputs/figures/supp_fig_s2_full_cell_composition.png`
- `outputs/figures/supp_fig_s2_full_cell_composition.pdf`
- `outputs/figures/comparison_original_supp_fig_s2_vs_pc_full.png`
- `outputs/tables/supp_fig_s2_dataset_manifest.tsv`
- `outputs/tables/supp_fig_s2_cell_type_composition.tsv`
- `outputs/tables/supp_fig_s2_umap_coordinates.tsv.gz`
- `outputs/tables/supp_fig_s2_annotation_columns_audit.tsv`
- `outputs/tables/supp_fig_s2_cluster_markers.tsv.gz`
- `outputs/tables/supp_fig_s2_final_cluster_to_global_cell_type.tsv`
- `outputs/logs/run_supp_fig_s2_full_cell_composition_pc.log`

If a dataset fails, do not replace it silently with another dataset. Upload the completed panels, the log, and a dataset-level failure table with the exact error text.

## Visual Target

Match the original page structure:

- Three rows.
- Two datasets per row.
- Each dataset has a UMAP on the left and a pie chart on the right.
- UMAP labels use global cell type names.
- Myeloid cells should be visibly labelled.
- Pie chart percentages should be readable and should not overlap.

Original row order:

1. `LUAD`, `LIHC-CHOL`
2. `Healthy lung`, `Healthy liver`
3. `NSCLC-GSE148071`, `HCC-GSE156625`

