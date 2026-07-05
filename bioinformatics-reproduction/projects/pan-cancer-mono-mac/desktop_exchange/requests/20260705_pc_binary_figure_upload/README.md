# 20260705 PC binary figure upload request

## Purpose

Do not rerun the analyses. Package and upload the PC-generated PNG figures that were reported as complete but were not included in the text-only GitHub handoffs.

The Mac needs these figures for local visual QC and original-versus-generated comparison.

## Run

From the repository root on the PC:

```powershell
cd I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac
powershell -ExecutionPolicy Bypass -File desktop_exchange\requests\20260705_pc_binary_figure_upload\scripts\run_package_pc_binary_figures.ps1
```

The script writes a text-only handoff under:

```text
desktop_exchange\uploads\20260705_pc_binary_figure_upload\pc_binary_figure_text_handoff
```

Commit and push that whole output directory. It contains base64 text chunks, a manifest, SHA256 hashes, and a reassembly script. This avoids the GitHub connector binary-upload limitation.

If direct binary upload works in your environment, also upload the original PNG files under:

```text
desktop_exchange\uploads\20260705_pc_binary_figure_upload\pc_binary_figure_text_handoff\figures_direct
```

The base64 text chunks are still the required upload.

## Required PNG files

The packaging script uses these exact PC paths:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_supp_fig_s2_full_cell_composition\pc_supp_fig_s2_full_cell_composition_results\outputs\figures\supp_fig_s2_full_cell_composition.png
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_supp_fig_s2_full_cell_composition\pc_supp_fig_s2_full_cell_composition_results\outputs\figures\comparison_original_supp_fig_s2_vs_pc_full.png
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_supp_fig_s4_s5_spatial_featureplots\pc_supp_fig_s4_s5_spatial_featureplots_results\figures\supp_fig_s4_s5_spatial_featureplots.png
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_supp_fig_s6_spatial_signatures\pc_supp_fig_s6_spatial_signatures_results\figures\supp_fig_s6_spatial_signatures.png
```

## Optional audit PNG

Upload this only if it exists. It is for record-keeping because the PC status reported this Fig. 5E route as an all-zero scientific nonmatch.

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_fig5e_gse120575_full_seurat\pc_fig5e_full_seurat_results\outputs\figures\fig5e_full_seurat_proportions.png
```

## Required output checks

After running the script, verify:

- `STATUS.md` says `required_found` is `4`.
- `MISSING_FILES.tsv` has no row with `priority` equal to `required`.
- `upload_manifest.tsv` records bytes and SHA256 for each packaged PNG.
- `files/base64_chunks/` contains one or more `.b64.txt` files per packaged PNG.

Do not upload raw GEO data, R packages, Seurat objects, extracted Visium folders, H5 files, or temporary caches for this request.
