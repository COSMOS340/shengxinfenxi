# TNBC Pre-Supplementary-Fig-S1 Portable Standard Object Request

## Purpose

The previous TNBC standard object upload reconstructed successfully on Mac, but its RNA assay layers still point to the PC-only BPCells directory:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\tnbc_full_run\outputs\04_data_processed\single_cell\full_scope_preprocess\tnbc_bpcells_full\GSE169246_TNBC_RNA_counts_bpcells
```

That object can be read for metadata on Mac, but any expression calculation fails when it tries to materialize the BPCells matrix. This request asks the PC to export a portable TNBC standard object whose `counts` and `data` layers are stored inside the RDS as sparse matrices.

This export also applies the Supplementary Fig. S1 scope correction before upload:

- keep only the 21 `Pre_P*` TNBC labels from the paper legend
- remove `P003`
- remove Post/Prog cells from the retained Fig. S1 patients
- expected retained cells: `30209`
- expected retained patients: `21`
- expected retained batches: `36`

## Run

From:

```bash
cd shengxinfenxi/bioinformatics-reproduction/projects/pan-cancer-mono-mac/tnbc_full_run
```

Run:

```bash
Rscript ../desktop_exchange/requests/20260702_tnbc_pre_supp_fig_s1_portable/export_tnbc_pre_supp_fig_s1_portable.R
```

Default input:

```text
outputs/04_data_processed/single_cell/full_scope_standard_objects/set2_TNBC_standard_myeloid.rds
```

Default output:

```text
outputs/04_data_processed/single_cell/full_scope_standard_objects/set2_TNBC_standard_myeloid_pre_supp_fig_s1_portable.rds
```

## Upload Back

Create:

```text
desktop_exchange/uploads/20260702_tnbc_pre_supp_fig_s1_portable/
```

Include:

- `STATUS.md`
- `upload_manifest.tsv`
- `checksums.sha256`
- `logs/export_tnbc_pre_supp_fig_s1_portable.log`
- `tables/set2_TNBC_pre_supp_fig_s1_portable_summary.tsv`
- `tables/set2_TNBC_pre_supp_fig_s1_portable_layer_check.tsv`
- split RDS chunks under `files/`

Use the existing split helper:

```bash
python3 ../desktop_exchange/requests/20260702_tnbc_standard_object/split_file_for_github.py \
  --input outputs/04_data_processed/single_cell/full_scope_standard_objects/set2_TNBC_standard_myeloid_pre_supp_fig_s1_portable.rds \
  --out-dir ../desktop_exchange/uploads/20260702_tnbc_pre_supp_fig_s1_portable/files \
  --chunk-size-mb 90
```

Do not upload raw GEO files or BPCells directories.

## Pass Criteria

- R script exits with code 0.
- Summary table reports `retained_cells=30209`, `retained_patients=21`, and `retained_batches=36`.
- Layer check table reports no BPCells classes for `counts` or `data` after reloading the saved RDS.
- Reconstructed RDS checksum matches the uploaded `.sha256` file.
