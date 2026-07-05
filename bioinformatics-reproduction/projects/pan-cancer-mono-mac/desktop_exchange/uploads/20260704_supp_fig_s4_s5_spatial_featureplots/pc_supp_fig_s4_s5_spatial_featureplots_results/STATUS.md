# Supplementary Fig. S4/S5 spatial feature plots PC run status

Run completed on 2026-07-05 Asia/Shanghai.

## Scope

Request: `20260704_supp_fig_s4_s5_spatial_featureplots`

Workflow executed on the PC:

1. Loaded 10x Genomics GBM Visium data.
2. Loaded 10x Genomics colorectal cancer Visium data.
3. Downloaded GSE274103 raw tar from NCBI GEO official FTP.
4. Verified `GSE274103_RAW.tar` by byte size and `tar -tf` listing.
5. Extracted five GSE274103 PDAC spatial folders and selected the first sorted folder for plotting under the request script's existing selection logic.
6. Generated spatial feature plots and summary tables.

## Dataset loading summary

Loaded datasets: 3

- `10x_GBM`: loaded, 10,878 spots, 18,085 features. Genes present: `EGFR`, `C1QC`, `FOLR2`.
- `10x_CRC`: loaded, 9,080 spots, 18,085 features. Genes present: `EPCAM`, `C1QC`, `FOLR2`.
- `GSE274103_PDAC`: loaded from `52906-10_spatial`, 4,134 spots, 17,943 features. Genes present: `EPCAM`, `C1QC`, `PLTP`.
- `GSE226997_CRC`: skipped by `S4S5_SKIP_GSE226997=1` because the 41.2 GB GEO RAW tar was not feasible for this run after observed single-connection download speed was only tens of KB/s.

## Generated outputs

Local handoff package:

`I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_supp_fig_s4_s5_spatial_featureplots\pc_supp_fig_s4_s5_spatial_featureplots_results`

Key files:

- `figures/supp_fig_s4_s5_spatial_featureplots.png`
- `figures/supp_fig_s4_s5_spatial_featureplots.pdf`
- `tables/supp_fig_s4_s5_dataset_manifest.tsv`
- `tables/supp_fig_s4_s5_gene_presence.tsv`
- `tables/supp_fig_s4_s5_gse274103_discovered_visium_dirs.tsv`
- `tables/supp_fig_s4_s5_selected_spatial_sources.tsv`
- `logs/run_supp_fig_s4_s5_spatial_featureplots.log`
- `logs/supp_fig_s4_s5_spatial_featureplots.console.log`
- `scripts/run_supp_fig_s4_s5_spatial_featureplots.R`
- `scripts/run_supp_fig_s4_s5_spatial_featureplots_pc.ps1`
- `files/binary_and_raw_files_not_uploaded.txt`

## Visual QA

The PNG was opened locally and visually inspected. The figure is 3150 x 3750 pixels. Titles, gene labels, colorbars, and panels are visible without obvious clipping or overlap. The blank `GSE226997_CRC not loaded` row is intentional and reflects the skipped dataset.

## Local script adjustments

- Switched GEO downloads from `www.ncbi.nlm.nih.gov/geo/download/?acc=...` to official GEO FTP URLs.
- Added reliable `curl.exe` fallback with retries and resume support.
- Disabled `aria2c` in the PC runner after `aria2c` stalled on GSE274103 at 19,031,408 bytes.
- Used verified byte-range `curl` chunks to complete GSE274103, then merged and verified the tar size and listing.
- Fixed `find_visium_dirs()` to discover GSE274103 `*_spatial` directories.
- Fixed flat GSE274103 spatial/H5 pairing by reading each `GSM..._spatial.tar.gz` top-level directory and matching it to its paired `GSM..._filtered_feature_bc_matrix.h5`.
- Adjusted the PC runner so R warnings/messages on stderr do not terminate the PowerShell process.

## Checksums

- PNG SHA256: `713f6f37789e7671ed00e53fef78d4dfe2cea43378ac021adc50970cb9d27640`
- PDF SHA256: `c983b0b3f2944047c7eb33dc9fa2b6248bbac7c01de4a1d4ef599d39b7c5672f`
- GSE274103 RAW tar SHA256: `3d77df521640e60cf43480c3e3d92ddb8628f771fd42f28e6b4f26d19c8522a8`

See `upload_manifest.tsv` for all local handoff package file hashes.
