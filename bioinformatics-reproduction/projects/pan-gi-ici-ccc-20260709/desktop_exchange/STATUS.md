# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-17

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

- Request: `20260717_gse205506_paper_qc_rebuild_request.md`
- Manifest: `20260717_gse205506_paper_qc_rebuild_manifest.tsv`
- Mac review: `20260717_gse205506_pc_run_review.md`
- PC status: `REQUESTED_GSE205506_PAPER_QC_REBUILD`
- Required upload: `uploads/20260717_gse205506_paper_qc_rebuild`
- Current message: Rebuild GSE205506 from the 40 author matrices using the Cancer Cell source-paper QC. The formal object must retain cells with 500-5000 genes and 400-25000 UMIs, use RPCA with 2000 variable genes, 20 PCs and first-round resolution 1.2, assign six broad compartments, and then apply compartment-aware mitochondrial filtering. Do not delete the nine high-mitochondrial clusters as whole units. scDblFinder is sensitivity-only because the source paper did not report an explicit computational doublet caller. Redraw all UMAPs with the approved bright palette copied from the first R/Seurat dataset.

## Most Recent Completed Request

- Completed request: `20260716_stop_prjna932556_and_build_gse205506_request.md`
- Completed upload: `uploads/20260716_stop_prjna932556_and_build_gse205506`
- Completion commit: `9a093b32ff4425a6b8ce14652b8c3ca7dcbb5803`
- Accepted: PRJNA932556 deletion audit; GSE205506 archive and matrix audit; Table S1 response mapping; R scripts, logs, session information, and author marker references.
- Not accepted for formal Figure 1: the 238,934-cell object, its UMAP, clusters, markers, and annotations because the preceding mitochondrial QC did not reproduce the source article and retained large populations with extreme mitochondrial fractions.

## Completed Upload Directories

- `uploads/20260716_stop_prjna932556_and_build_gse205506`
- `uploads/20260714_cellranger72_archive_downloaded_to_f`
- `uploads/20260714_prjna932556_cellranger72_single_sample_hardrun`
- `uploads/20260714_prjna932556_cellranger72_preflight`
- `uploads/20260714_prjna932556_srr23490337_full_bamtofastq`
- `uploads/20260713_prjna932556_bamtofastq_feasibility`
- `uploads/20260713_prjna932556_ena_bam_pair_download`
- `uploads/20260711_prjna932556_ena_bam_acquisition`
- `uploads/20260711_gse235863_timepoint_prjna932556`
- `uploads/20260710_gse189926_r_only_rerun`
- `uploads/20260710_gse189926_annotation_qc_refinement`
- `uploads/20260709_metadata_download`
- `uploads/20260709_small_file_download`
- `uploads/20260709_matrix_download`
- `uploads/20260709_preprocess_priority_response`
- `uploads/20260709_response_object_construction`
