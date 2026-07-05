# PC Supp Fig. S6 Spatial Signatures Results

Status: complete on Windows desktop.

## Scope

Request: `desktop_exchange/requests/20260704_supp_fig_s6_spatial_signatures`

Workflow run: `scripts/run_supp_fig_s6_spatial_signatures_pc.ps1`

The PC run reused the vetted spatial inputs from the Supp Fig. S4/S5 run and generated spatial feature/signature panels plus SPP1+ TAM versus MDSC signature correlations.

## Data Loaded

| dataset_id | status | notes |
|---|---|---|
| `10x_GBM` | loaded | 10878 spots |
| `GSE274103_PDAC` | loaded | 4134 spots; standardized from GEO GSE274103 |
| `10x_CRC` | loaded | 9080 spots |
| `GSE226997_CRC` | skipped_by_env | GEO supplementary tar is 41.2 GB and was not downloaded for this desktop handoff |

## Main Results

| dataset_id | Spearman rho | p_value | spots |
|---|---:|---:|---:|
| `10x_GBM` | 0.568594053301822 | 0 | 10878 |
| `GSE274103_PDAC` | 0.222422516105431 | 1.66105695894019e-47 | 4134 |
| `10x_CRC` | 0.612900497126217 | 0 | 9080 |

Run summary:

| metric | value |
|---|---:|
| datasets_loaded | 3 |
| datasets_planned | 4 |
| spp1_signature_genes | 10 |
| mdsc_signature_genes | 19 |
| correlation_panels | 3 |

## Output Files

Key local files are listed in `upload_manifest.tsv`.

Large or binary outputs retained locally:

| file | bytes | sha256 |
|---|---:|---|
| `figures/supp_fig_s6_spatial_signatures.png` | 6870800 | `7064f958d8205ca7ea59c41b7b4675566636f4785975b0dd6b25c314f810f584` |
| `figures/supp_fig_s6_spatial_signatures.pdf` | 10317028 | `c3d3d18dd04c5034b95b461ca6b980a86390c5922792eba850eb62e45b429a4e` |
| `tables/supp_fig_s6_signature_scores.tsv.gz` | 1225600 | `42e839b66c4dac23ceb860e9ccc6453b85c618b7f9c06cf08dd94e5479397e26` |

## QA

- R workflow log ends with `save_figure complete` and `done complete`.
- Visual inspection of the regenerated PNG passed: title, labels, colorbars, scatter panels, and skipped-dataset placeholders are visible without obvious clipping or overlap.
- R package build-version warnings were observed and recorded in `logs/s6_runner2_console.log`; no workflow failure was recorded by the R log.
- Raw GEO inputs, extracted Visium folders, and BPCells-style caches were not included in this handoff package.
