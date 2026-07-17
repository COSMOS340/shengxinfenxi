# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-17

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Requests

- Priority 1 quick sweep: `20260717_gse189926_resolution_umap_topology_sweep_request.md`
- Priority 1 manifest: `20260717_gse189926_resolution_umap_topology_sweep_manifest.tsv`
- Priority 2 long rebuild: `20260717_gse205506_paper_qc_rebuild_request.md`
- Priority 2 manifest: `20260717_gse205506_paper_qc_rebuild_manifest.tsv`
- Mac review: `20260717_gse205506_pc_run_review.md`
- PC status: `REQUESTED_GSE189926_SWEEP_THEN_GSE205506_REBUILD`
- Required upload 1: `uploads/20260717_gse189926_resolution_umap_topology_sweep`
- Required upload 2: `uploads/20260717_gse205506_paper_qc_rebuild`
- Current message: First run the lightweight GSE189926 fixed-coordinate resolution sweep and fixed-Harmony UMAP parameter sweep. The existing 0.3, 0.5 and 0.8 resolutions share one UMAP, so resolution and UMAP topology must be evaluated separately. Upload all seven resolution panels, seven UMAP routes, local-linearity metrics, sample and patient mixing, QC overlays, and long-arm cell audits. Then run the longer GSE205506 source-paper QC rebuild. For GSE205506, do not delete the nine high-mitochondrial clusters as whole units; scDblFinder remains sensitivity-only; redraw UMAPs with the approved bright palette.

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
