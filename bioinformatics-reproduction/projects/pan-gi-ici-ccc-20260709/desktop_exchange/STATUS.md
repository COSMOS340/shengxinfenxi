# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-13

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

Please run:

- `20260713_prjna932556_bamtofastq_feasibility_request.md`
- `20260713_prjna932556_bamtofastq_feasibility_manifest.tsv`

Upload lightweight outputs to:

- `uploads/20260713_prjna932556_bamtofastq_feasibility`

The task is a bounded feasibility pilot for official 10x `bamtofastq` using the already verified `SRR23490337` BAM/BAI pair. Reverify the BAM/BAI first, audit runtime and tool versions, then run only the locus-restricted pilot if the exact contig rule passes. Keep all FASTQ outputs PC-local and do not run Cell Ranger, Seurat, CellChat, LIANA, NicheNet, or quantification.

## Most Recent Completed Request

The previous request, `20260713_prjna932556_ena_bam_pair_download_request.md`, is complete. The ENA gate passed, and both exact submitted files for `SRR23490337` were downloaded to PC-local storage and verified by byte size and MD5. The BAM and BAI remain PC-local and were not uploaded to GitHub.

## Completed Request Files

- `20260713_prjna932556_ena_bam_pair_download_request.md`
- `20260713_prjna932556_ena_bam_pair_download_manifest.tsv`
- `20260711_prjna932556_ena_bam_acquisition_request.md`
- `20260711_prjna932556_ena_bam_acquisition_manifest.tsv`
- `20260711_prjna932556_windows_sra_pilot_request.md`
- `20260711_prjna932556_windows_sra_pilot_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_request.md`

## Completed Upload Directories

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
