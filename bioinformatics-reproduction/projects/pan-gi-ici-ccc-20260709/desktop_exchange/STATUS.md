# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-14

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

- Request: `20260714_prjna932556_cellranger72_preflight_request.md`
- Manifest: `20260714_prjna932556_cellranger72_preflight_manifest.tsv`
- Scope: audit exact host, Linux runtime, Cell Ranger 7.2.0, GRCh38-2020-A, and storage readiness for one later `SRR23490337` count run
- Required decision: `LOCAL_CELLRANGER_COUNT_NOT_READY` or `READY_FOR_CELLRANGER72_ACQUISITION_AND_COUNT_REQUEST`
- Stop boundary: no large downloads, no installation, no count, no other PRJNA932556 runs

## Most Recent Completed Request

The previous request, `20260714_prjna932556_srr23490337_full_bamtofastq_request.md`, is complete. The existing `SRR23490337` BAM/BAI pair was reverified, the official 10x `bamtofastq` v1.4.1 binary was used without `--locus`, and full conversion exited with code 0. The run wrote 401831188 read pairs to 804 gzip FASTQ files totaling 29268262419 bytes. Full FASTQ files remain PC-local and were not uploaded to GitHub.

## Completed Upload Directories

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
