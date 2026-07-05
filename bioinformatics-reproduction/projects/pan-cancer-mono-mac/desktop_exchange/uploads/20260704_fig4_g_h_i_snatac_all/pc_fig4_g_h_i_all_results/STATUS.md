# Fig4G-H-I snATAC all-scope PC run status

Completed on 2026-07-05 Asia/Shanghai.

## Request

Request directory:

`desktop_exchange/requests/20260704_fig4_g_h_i_snatac_all`

Output package:

`desktop_exchange/uploads/20260704_fig4_g_h_i_snatac_all/pc_fig4_g_h_i_all_results`

## Result

Status: complete.

Postrun checker passed:

- required files nonempty: 13
- run_scope: all
- requested_motifs_found: 6
- requested motif rows: 6
- motif activity rows: 54
- track gene rows: 4
- processed RDS policy: not_saved

Run summary:

- cells: 122218
- peaks: 1307563
- libraries: 15
- samples: 45
- transferred labels: 9
- transfer features: 2675
- differential peak rows: 30984
- requested motifs found exactly: 6

## Required outputs

Generated required Fig4G-H-I files:

- `figures/fig4g_gene_coverage_tracks.png`
- `figures/fig4g_gene_coverage_tracks.pdf`
- `figures/fig4h_motif_accessibility_umaps.png`
- `figures/fig4h_motif_accessibility_umaps.pdf`
- `figures/fig4i_motif_logos.png`
- `figures/fig4i_motif_logos.pdf`
- `tables/fig4_motif_lookup.tsv`
- `tables/fig4_requested_motif_ids.tsv`
- `tables/fig4e_chromvar_motif_activity_by_label.tsv`
- `tables/fig4g_track_region_summary.tsv`

Regenerated A-C files are also included in the package.

## PC adjustments

The original request logic was adapted for this PC environment:

- Targeted gene activity was computed for reference-overlapping transfer genes to avoid the previous all-gene memory failure.
- JASPAR motif matching was restricted to the six exact requested motif names after writing the full motif lookup table.
- `chromVAR::computeDeviations()` was used directly because local Signac 1.17.1 does not export `RunChromVAR`.
- `BiocParallel::SerialParam()` was used to avoid socket serialization failures with the large peak matrix.
- Fig4G coverage tracks were completed by a coverage-only resume script using the same local fragments and predicted labels from the completed all-scope run.

## Visual QA

PNG figures were visually inspected.

- Fig4G coverage tracks: readable; gene tracks and peak tracks visible. Some x-axis labels are dense in the four-panel layout.
- Fig4H motif UMAPs: readable; titles now use motif names rather than motif IDs.
- Fig4I motif logos: readable; no clipping observed.
- Regenerated Fig4A-C: readable; no legend clipping observed.

## Scientific caveats

- The all-scope workflow used all 122218 cells from the matched H5/barcode set and did not apply the formal Fig4A QC filters used in the separate formal Fig4A rerun.
- The UMAP contains long track-like structures, consistent with prior GSE306459 snATAC runs; interpretation should consider possible library or technical structure.
- CoveragePlot emitted warnings about a small number of gene annotation segments outside the displayed regions; the requested gene coverage tracks were still generated.

## Not uploaded

Raw H5 files, fragment files, tabix index files, the raw GEO archive, R packages, and processed RDS objects are not included.
