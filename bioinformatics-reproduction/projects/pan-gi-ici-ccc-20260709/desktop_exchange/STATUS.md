# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-14

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

Please run:

- `20260714_prjna932556_srr23490337_full_bamtofastq_request.md`
- `20260714_prjna932556_srr23490337_full_bamtofastq_manifest.tsv`

Upload lightweight outputs to:

- `uploads/20260714_prjna932556_srr23490337_full_bamtofastq`

The task is full official 10x `bamtofastq` conversion for `SRR23490337` only. Reverify the already downloaded BAM/BAI pair, enforce the 60 GB free-space gate, run full conversion without `--locus`, and upload only lightweight audits plus first-100-record summaries. Keep full FASTQ files PC-local. Do not process other PRJNA932556 runs and do not run Cell Ranger, Seurat, CellChat, LIANA, NicheNet, or quantification.

## Most Recent Completed Request

The previous request, `20260713_prjna932556_bamtofastq_feasibility_request.md`, is complete. The existing `SRR23490337` BAM/BAI pair was reverified, the exact contig rule selected `MT:1-16569`, and official 10x `bamtofastq` v1.4.1 completed the bounded pilot with exit code 0. The run wrote 76,559,661 read pairs to 154 gzip FASTQ files totaling 5,386,519,028 bytes. Full FASTQ files remain PC-local and were not uploaded to GitHub.

## Completed Request Files

- `20260713_prjna932556_bamtofastq_feasibility_request.md`
- `20260713_prjna932556_bamtofastq_feasibility_manifest.tsv`
- `20260713_prjna932556_ena_bam_pair_download_request.md`
- `20260713_prjna932556_ena_bam_pair_download_manifest.tsv`
- `20260711_prjna932556_ena_bam_acquisition_request.md`
- `20260711_prjna932556_ena_bam_acquisition_manifest.tsv`
- `20260711_prjna932556_windows_sra_pilot_request.md`
- `20260711_prjna932556_windows_sra_pilot_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_request.md`

## Completed Upload Directories

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
