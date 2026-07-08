# PC Request: CRLM CAF/ECM Raw Structure Audit

Created: 2026-07-08

## Context

The previous PC task completed all required raw GEO downloads for the CRLM CAF/ECM project. The raw files are stored locally on PC and were intentionally not committed to GitHub.

Previous PC raw directory:

`I:\shengxinfenxi\bioinformatics-reproduction\projects\crlm-caf-ecm-spatial-niche\desktop_exchange\requests\20260708_core_geo_download\raw_core_geo`

Completed raw files from PC status:

| Dataset | File | Bytes | SHA256 | Integrity |
|---|---|---:|---|---|
| `GSE225857` | `GSE225857_RAW.tar` | 636508160 | `f697d211c0c9239c79f0044f4777adb7091e4639e2a92a80122874e148629d42` | `tar_pass` |
| `GSE178318` | `GSE178318_matrix.mtx.gz` | 545942023 | `304943daca7795ab5fef630d6e127d385a8469ccbfb3d0c196788d9102f8cd74` | `gzip_pass` |
| `GSE245552` | `GSE245552_RAW.tar` | 1065748480 | `012c4c3e44e763bd63a14808e5f6ba17829239ee39ba7eb19bac1ff0671b4c1a` | `tar_pass` |

Mac has only the two small `GSE178318` files and cannot run raw-dependent analysis locally yet.

## Task

Run the raw structure audit script:

```bash
python bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_raw_structure_audit/scripts/pc_audit_raw_structures.py
```

The script should read the existing PC raw directory from the previous request and write small TSV outputs under:

`bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/uploads/20260708_raw_structure_audit/`

## Required Outputs

Commit only these small outputs:

- `STATUS.md`
- `raw_file_manifest.tsv`
- `gse225857_tar_members.tsv`
- `gse225857_inner_file_summaries.tsv`
- `gse178318_matrix_dimensions.tsv`
- `gse245552_tar_members.tsv`
- `gse245552_inner_file_summaries.tsv`

Do not commit raw `.tar`, `.mtx.gz`, large extracted matrices, Seurat objects, or image files.

## Interpretation Boundary

This request is only file and matrix structure audit.

Do not assign cell types.
Do not infer patient IDs or tissue labels from partial filename matching.
Do not decide final inclusion/exclusion.

The output should preserve exact raw filenames and exact GEO accessions where present.

