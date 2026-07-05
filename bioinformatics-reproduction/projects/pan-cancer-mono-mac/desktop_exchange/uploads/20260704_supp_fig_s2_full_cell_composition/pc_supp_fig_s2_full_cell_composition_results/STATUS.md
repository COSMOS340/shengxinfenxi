# Supplementary Fig. S2 full cell composition PC run

Status: complete

Run time: 2026-07-05 16:34:34 to 2026-07-05 16:46:45 Asia/Shanghai

Local result directory:

`I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_supp_fig_s2_full_cell_composition\pc_supp_fig_s2_full_cell_composition_results`

## Completion

- LUAD: 88,144 cells, loaded
- LIHC-CHOL: 9,946 cells, loaded
- Healthy lung: 43,632 cells, loaded
- Healthy liver: 8,444 cells, loaded
- NSCLC-GSE148071: 89,887 cells, loaded
- HCC-GSE156625: 73,589 cells, loaded
- Dataset failures: 0

## Annotation sources

- LUAD used `GSE131907_Lung_Cancer_cell_annotation.txt.gz` field `Cell_type.refined`, filtered by `Sample_Origin` values `tLung` and `nLung`.
- LIHC-CHOL used `GSE125449_Set1_samples.txt.gz` and `GSE125449_Set2_samples.txt.gz` field `Type`.
- Healthy liver used `GSE115469_CellClusterType.txt.gz` field `CellType`.
- Healthy lung used the eight GSE122960 donor filtered h5 files and marker-score assignment because no direct global cell type annotation was observed in the downloaded files.
- NSCLC-GSE148071 used 42 `*_exp.txt.gz` files and marker-score assignment because no direct global cell type annotation was observed in the downloaded files.
- HCC-GSE156625 used `GSE156625_HCCscanpyobj.h5ad` with `louvain` clusters mapped from `rank_genes_groups` markers.

## Visual QA

Visual inspection passed for `outputs/figures/supp_fig_s2_full_cell_composition.png`.

The first rendered version had overlapping UMAP inline labels. The final version uses major global cell type labels with white label boxes and connector lines. `Myeloid cells` is visibly labelled in all six UMAP panels. The final figure has no observed text overlap, clipped axis labels, clipped legend, or canvas-size problem.

## Included outputs

- `outputs/figures/supp_fig_s2_full_cell_composition.png`
- `outputs/figures/supp_fig_s2_full_cell_composition.pdf`
- `outputs/figures/comparison_original_supp_fig_s2_vs_pc_full.png`
- `outputs/tables/supp_fig_s2_dataset_manifest.tsv`
- `outputs/tables/supp_fig_s2_cell_type_composition.tsv`
- `outputs/tables/supp_fig_s2_umap_coordinates.tsv.gz`
- `outputs/tables/supp_fig_s2_annotation_columns_audit.tsv`
- `outputs/tables/supp_fig_s2_cluster_markers.tsv.gz`
- `outputs/tables/supp_fig_s2_final_cluster_to_global_cell_type.tsv`
- `outputs/tables/supp_fig_s2_dataset_failures.tsv`
- `outputs/logs/run_supp_fig_s2_full_cell_composition_pc.log`
- `scripts/run_supp_fig_s2_full_cell_composition_pc.py`

Raw GEO inputs, extracted input directories, h5ad files, tar files, and temporary caches are not included.

## Note

`comparison_original_supp_fig_s2_vs_pc_full.png` currently contains the PC-generated full S2 panel because no original Supp Fig. S2 image was available in the local request package.
