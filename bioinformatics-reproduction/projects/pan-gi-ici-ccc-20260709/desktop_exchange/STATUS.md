# Pan-GI ICI CCC Metadata Download

Created: 2026-07-09

This folder contains the PC-side request for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

Download metadata and file lists only.

Do not download large expression matrices, FASTQ files, h5ad objects, RDS objects, or spatial image files in this task.

## Files

- `20260709_metadata_download_manifest.tsv`: exact resource list, URLs, and output paths.
- `20260709_metadata_download_request.md`: instructions, required output tables, and checksum requirements.

## Required PC Outputs

Write outputs to:

`/Volumes/research/pan_gi_ici_ccc_20260709/00_metadata/source_metadata`

Required completion files:

- `download_status.tsv`
- `SHA256SUMS.txt`

After completion, upload the status table, checksum file, and downloaded metadata files back through this repository or the existing desktop exchange route.
