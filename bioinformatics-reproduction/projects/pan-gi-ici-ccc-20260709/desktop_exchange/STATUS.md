# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-14

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

- Request: `20260714_prjna932556_cellranger72_preflight_request.md`
- Manifest: `20260714_prjna932556_cellranger72_preflight_manifest.tsv`
- PC result: `LOCAL_CELLRANGER_COUNT_NOT_READY`
- Audit upload: `uploads/20260714_prjna932556_cellranger72_preflight`
- Stop boundary honored: no large downloads, no installation, no Cell Ranger count, no other PRJNA932556 runs

## Most Recent Completed Request

The current request, `20260714_prjna932556_cellranger72_preflight_request.md`, is complete as an audit. The PC-local `SRR23490337` FASTQs were reverified. The current QEMU audit environment does not satisfy the requested Cell Ranger count gate because it exposes 4 guest CPU cores and 4104351744 guest RAM bytes, and the largest fixed-volume free space observed is 417535246336 bytes. The decision is `LOCAL_CELLRANGER_COUNT_NOT_READY`.

## Completed Upload Directories

- `uploads/20260714_prjna932556_cellranger72_preflight` (Cell Ranger 7.2.0 preflight audit)
- `uploads/20260714_prjna932556_srr23490337_full_bamtofastq` (full bamtofastq success audit; full FASTQs PC-local only)
- `uploads/20260713_prjna932556_bamtofastq_feasibility`
- `uploads/20260713_prjna932556_ena_bam_pair_download`
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
