# GSE164522 Macrophage/Monocyte doublet and batch audit

Completed: 2026-07-06T22:18:41+0800
Seurat object: I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-macrophage-atlas\desktop_exchange\requests\20260706_gse164522_full_mac_mono_analysis\local_objects_not_for_github\gse164522_full_mac_mono_seurat.rds
Cells: 15853
Genes: 24662
Cluster column: RNA_snn_res.0.6
Clusters: 23

No clusters were removed in this audit.

Package availability:
SingleCellExperiment: TRUE
scDblFinder: TRUE
BiocParallel: TRUE
FNN: TRUE
harmony: TRUE

Doublet status: scDblFinder completed
Nearest-neighbor batch mixing status: FNN nearest-neighbor batch mixing completed
Harmony status: Harmony by sample completed

Review thresholds written as boolean flags:
sample_dominated_70: max_sample_fraction >= 0.70
patient_dominated_70: max_patient_fraction >= 0.70
tissue_dominated_90: max_tissue_fraction >= 0.90
author_subtype_dominated_80: max_celltype_sub_fraction >= 0.80
doublet_enriched_15: scDblFinder_doublet_fraction >= 0.15
non_myeloid_score_enriched: top_non_myeloid_fraction >= 0.50

Main tables:
gse164522_audit_cluster_integrated_review.tsv
gse164522_audit_cluster_review_status.tsv

PC local fixes applied before upload:
- Cast QC medians and means to numeric to avoid mixed integer/double group output.
- Used Seurat v5 `layer = "data"` access with legacy fallback.
- Fixed data.table multi-column assignment syntax for non-myeloid score reset.
- Regenerated UMAP and Harmony UMAP PNGs with larger opaque points and higher-contrast colors after visual QC.

PC visual QC:
- All PNG files were opened locally through a contact sheet.
- Plots were not blank.
- Axes, legends, and main labels were not clipped.
- Initial UMAP colors were too pale; the UMAP PNGs were regenerated with higher-contrast colors before packaging.
