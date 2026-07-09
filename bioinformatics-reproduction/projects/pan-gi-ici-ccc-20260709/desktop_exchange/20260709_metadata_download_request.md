# PC Task: Metadata Download for Pan-GI ICI CCC

Created: 2026-07-09

## Why PC Is Needed

This Mac failed at DNS resolution for Zenodo, CELLxGENE, and NCBI FTP/API during metadata download. The failure happened before file transfer, so it does not indicate that the resources are unavailable.

## Goal

Download only metadata and file lists for the current project. Do not download large count matrices, FASTQ files, RDS objects, h5ad objects, or spatial image files in this task.

## Project Paths

Project root:

`/Volumes/research/pan_gi_ici_ccc_20260709`

Metadata output directory:

`/Volumes/research/pan_gi_ici_ccc_20260709/00_metadata/source_metadata`

Manifest:

`/Volumes/research/pan_gi_ici_ccc_20260709/99_logs/pc_tasks/20260709_metadata_download_manifest.tsv`

## Required Outputs

1. Download every row in the manifest where possible.
2. Preserve the exact output paths in the manifest.
3. Generate checksum file:

`/Volumes/research/pan_gi_ici_ccc_20260709/00_metadata/source_metadata/SHA256SUMS.txt`

4. Generate download status table:

`/Volumes/research/pan_gi_ici_ccc_20260709/00_metadata/source_metadata/download_status.tsv`

Required columns for `download_status.tsv`:

- resource_id
- status
- output_path
- file_size_bytes
- sha256
- error_message

## Suggested Shell Workflow

Use the manifest as the source of truth. For rows with `url` not equal to `manual`, run:

```bash
mkdir -p /Volumes/research/pan_gi_ici_ccc_20260709/00_metadata/source_metadata
while IFS=$'\t' read -r resource_id resource_type url output_path priority notes; do
  if [ "$resource_id" = "resource_id" ]; then
    continue
  fi
  if [ "$url" = "manual" ]; then
    printf '%s\tmanual_required\t%s\t0\t\tmanual export needed\n' "$resource_id" "$output_path"
    continue
  fi
  mkdir -p "$(dirname "$output_path")"
  if curl -L --retry 3 --connect-timeout 30 --max-time 300 -o "$output_path" "$url"; then
    bytes=$(wc -c < "$output_path" | tr -d ' ')
    sum=$(shasum -a 256 "$output_path" | awk '{print $1}')
    printf '%s\tdownloaded\t%s\t%s\t%s\t\n' "$resource_id" "$output_path" "$bytes" "$sum"
  else
    printf '%s\tfailed\t%s\t0\t\tcurl_failed\n' "$resource_id" "$output_path"
  fi
done < /Volumes/research/pan_gi_ici_ccc_20260709/99_logs/pc_tasks/20260709_metadata_download_manifest.tsv \
  > /Volumes/research/pan_gi_ici_ccc_20260709/00_metadata/source_metadata/download_status.tsv

cd /Volumes/research/pan_gi_ici_ccc_20260709/00_metadata/source_metadata
find . -maxdepth 1 -type f ! -name SHA256SUMS.txt -print0 | xargs -0 shasum -a 256 > SHA256SUMS.txt
```

## Manual Items

`STAD-PRJEB25780` needs CIDE metadata. If no direct API is obvious, use the CIDE web UI and export the sample table with response labels. Save it to:

`/Volumes/research/pan_gi_ici_ccc_20260709/00_metadata/source_metadata/STAD_PRJEB25780_CIDE_metadata.tsv`

For `OMIX001073`, save the HTML page and any sample or file table visible on the OMIX page.

## Important

Do not start heavy analysis in this task. The only purpose is to get metadata that lets us decide which datasets are scientifically usable.
