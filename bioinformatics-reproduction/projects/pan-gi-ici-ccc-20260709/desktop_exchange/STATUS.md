# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-14

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

- Request: `20260714_prjna932556_srr23490337_full_bamtofastq_request.md`
- Manifest: `20260714_prjna932556_srr23490337_full_bamtofastq_manifest.tsv`
- PC result: `STOPPED_BEFORE_CONVERSION`
- Precise stop reason: the required 60 GB disk gate failed; `I:` had `36503392256` free bytes, a shortfall of `23496607744` bytes.
- Audit upload: `uploads/20260714_prjna932556_srr23490337_full_bamtofastq`

The BAM/BAI pair and official 10x `bamtofastq` v1.4.1 binary were reverified successfully. The full command was not executed, no output directory was created, and no existing data was deleted. Free at least `23496607744` additional bytes on `I:` and rerun this same request.

## Most Recent Completed Request

The previous request, `20260713_prjna932556_bamtofastq_feasibility_request.md`, is complete. The existing `SRR23490337` BAM/BAI pair was reverified, the exact contig rule selected `MT:1-16569`, and official 10x `bamtofastq` v1.4.1 completed the bounded pilot with exit code 0. The run wrote 76,559,661 read pairs to 154 gzip FASTQ files totaling 5,386,519,028 bytes. Full FASTQ files remain PC-local and were not uploaded to GitHub.

## Completed Upload Directories

- `uploads/20260714_prjna932556_srr23490337_full_bamtofastq` (disk-gate audit; conversion not run)
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
