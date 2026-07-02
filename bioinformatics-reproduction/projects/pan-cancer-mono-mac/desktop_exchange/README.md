# Desktop Exchange for Pan-Cancer Mono/Mac Reproduction

This folder is the GitHub handoff point between the desktop machine and Codex for the pan-cancer monocyte/macrophage reproduction project.

Use it for desktop-run outputs, status notes, checksums, and short machine-readable manifests. Do not upload raw GEO inputs or temporary cache folders here.

## Upload Layout

After a desktop run completes, create one dated upload folder under `uploads/`:

```text
uploads/
  20260702_tnbc_full_run/
    STATUS.md
    upload_manifest.tsv
    checksums.sha256
    logs/
    tables/
    figures/
```

Use a new upload folder for every rerun. Do not overwrite a previous upload folder.

## Required Files

Each upload folder should contain:

- `STATUS.md`: copy from `STATUS_TEMPLATE.md` and fill every section.
- `upload_manifest.tsv`: copy from `upload_manifest_template.tsv` and list every uploaded file.
- `checksums.sha256`: SHA256 checksums for all uploaded files except `checksums.sha256` itself.
- `logs/`: run logs, especially the `Rscript ... | tee ...` output.
- `tables/`: summary tables, cluster markers, manual-review inputs, UMAP coordinates.
- `figures/`: generated PNG/PDF figures for visual QA.

Large raw matrices, downloaded GEO files, BPCells directories, and temporary files should stay outside this repository. If a required output is larger than GitHub can store safely, write its external path and checksum in `upload_manifest.tsv` instead of committing the file.

## Manifest Columns

Do not rename these columns:

```text
relative_path	file_type	dataset	description	bytes	sha256	created_by	source_command	notes
```

Codex will parse these exact field names.

## Minimum TNBC Upload

For the `tnbc_full_run` workflow, upload at least:

- `STATUS.md`
- `upload_manifest.tsv`
- `checksums.sha256`
- the full run log
- `tnbc_full_run_summary.tsv`
- `tnbc_full_cluster_top10_markers.tsv`
- `tnbc_full_cluster_markers_all.tsv`, if file size allows
- `tnbc_full_umap_coordinates.tsv`, if file size allows
- `tnbc_full_umap_by_cluster.png`
- `tnbc_full_canonical_marker_dotplot.png`

If marker testing was run with `TNBC_MARKER_MAX_CELLS_PER_CLUSTER`, record the exact value in `STATUS.md`.

## Codex Review Steps

After an upload appears in this folder, Codex should:

1. Verify `checksums.sha256`.
2. Read `STATUS.md` and `upload_manifest.tsv`.
3. Inspect the log for failed steps or marker-test limits.
4. Visually inspect uploaded figures.
5. Run marker review from `tnbc_full_cluster_top10_markers.tsv` and canonical marker plots.
6. Add an incremental original-vs-generated comparison entry before continuing to the next dataset.
