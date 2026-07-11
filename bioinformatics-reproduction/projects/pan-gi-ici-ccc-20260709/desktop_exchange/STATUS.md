# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-11

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

Acquire and verify the exact ENA-submitted BAM for `SRR23490337` from `PRJNA932556`.

Current request files:

- `20260711_prjna932556_ena_bam_acquisition_request.md`
- `20260711_prjna932556_ena_bam_acquisition_manifest.tsv`

Required lightweight upload directory:

- `uploads/20260711_prjna932556_ena_bam_acquisition`

Download only the single ENA-submitted BAM linked to `SRR23490337`. The BAM remains PC-local. Do not run `bamtofastq`, Cell Ranger, Seurat, or quantification in this task.

The preceding Windows SRA Toolkit pilot completed without producing FASTQ files. Its result is an archive-access failure, not a read-structure determination.

## Completed Request Files

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

- `uploads/20260711_prjna932556_windows_sra_pilot`
- `uploads/20260711_gse235863_timepoint_prjna932556`
- `uploads/20260710_gse189926_r_only_rerun`
- `uploads/20260710_gse189926_annotation_qc_refinement`
- `uploads/20260709_metadata_download`
- `uploads/20260709_small_file_download`
- `uploads/20260709_matrix_download`
- `uploads/20260709_preprocess_priority_response`
- `uploads/20260709_response_object_construction`

Large RDS/h5ad objects and any full FASTQ outputs remain PC-local and are recorded in the relevant upload inventories when generated.
