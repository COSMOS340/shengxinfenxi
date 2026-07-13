# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-13

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

Please run:

- `20260713_prjna932556_ena_bam_pair_download_request.md`
- `20260713_prjna932556_ena_bam_pair_download_manifest.tsv`

Upload lightweight outputs to:

- `uploads/20260713_prjna932556_ena_bam_pair_download`

The task is to download and verify the exact ENA-submitted BAM/BAI pair for `SRR23490337`. The BAM and BAI files must remain PC-local and must not be uploaded to GitHub. Do not run `bamtofastq`, Cell Ranger, Seurat, or quantification in this task.

## Most Recent Completed Request

The previous request, `20260711_prjna932556_ena_bam_acquisition_request.md`, was completed as a gate-failed audit. The ENA read-run report for `SRR23490337` was retrieved exactly, but Phase 2 failed because ENA returned `submitted_format=BAM;BAI` with two submitted file, checksum, and byte-size entries. Per that request, no submitted file was manually selected and no BAM was downloaded.

## Completed Request Files

- `20260711_prjna932556_ena_bam_acquisition_request.md`
- `20260711_prjna932556_ena_bam_acquisition_manifest.tsv`
- `20260711_prjna932556_windows_sra_pilot_request.md`
- `20260711_prjna932556_windows_sra_pilot_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_request.md`
- `20260710_gse189926_r_only_rerun_manifest.tsv`
- `20260710_gse189926_r_only_rerun_request.md`
- `20260710_gse189926_annotation_qc_refinement_manifest.tsv`
- `20260710_gse189926_annotation_qc_refinement_request.md`
- `20260710_gse189926_umap_repair_policy.md`

## Completed Upload Directories

- `uploads/20260711_prjna932556_ena_bam_acquisition`
- `uploads/20260711_prjna932556_windows_sra_pilot`
- `uploads/20260711_gse235863_timepoint_prjna932556`
- `uploads/20260710_gse189926_r_only_rerun`
- `uploads/20260710_gse189926_annotation_qc_refinement`
- `uploads/20260709_metadata_download`
- `uploads/20260709_small_file_download`
- `uploads/20260709_matrix_download`
- `uploads/20260709_preprocess_priority_response`
- `uploads/20260709_response_object_construction`

Large RDS/h5ad objects and any full FASTQ/BAM outputs remain PC-local and are recorded in the relevant upload inventories when generated.
