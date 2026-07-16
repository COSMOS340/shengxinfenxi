# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-16

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

- Request: `20260716_stop_prjna932556_and_build_gse205506_request.md`
- Manifest: `20260716_stop_prjna932556_and_build_gse205506_manifest.tsv`
- PC status: `REQUESTED_STOP_DELETE_PRJNA932556_SWITCH_GSE205506`
- Required upload: `uploads/20260716_stop_prjna932556_and_build_gse205506`
- Current message: Stop all PRJNA932556 compute. Inventory exact PRJNA932556 paths, delete its verified BAM, BAI, FASTQ, Cell Ranger work directories and expression outputs, and report bytes released. Keep Cell Ranger software, reference genomes, other datasets and lightweight audit files. Then download the author-processed GSE205506 archive, obtain the actual article Table S1 response mapping and build the formal CRC object entirely in R.

## Most Recent Completed Request

The Cell Ranger 7.2.0 archive download and prior PRJNA932556 BAM-to-FASTQ work are complete as technical records. The project has now removed PRJNA932556 from formal analysis because the six-sample reconstruction is not acceptable for current PC storage and compute. The new request explicitly requires deletion of PRJNA932556 large data after an exact-path inventory.

## Completed Upload Directories

- `uploads/20260714_cellranger72_archive_downloaded_to_f` (CELLRANGER72_ARCHIVE_DOWNLOADED_TO_F)
- `uploads/20260714_prjna932556_cellranger72_single_sample_hardrun` (BLOCKED_CELLRANGER72_MANUAL_DOWNLOAD)
- `uploads/20260714_prjna932556_cellranger72_preflight` (Cell Ranger 7.2.0 preflight audit)
- `uploads/20260714_prjna932556_srr23490337_full_bamtofastq` (full bamtofastq success audit; full FASTQs PC-local only)
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
