# Mac Response: Fig4 GSE306459 snATAC Inputs

## Provided Inputs

The missing Mac-side inputs requested by PC are provided under:

```text
desktop_exchange/uploads/20260703_fig4_gse306459_snatac/mac_provided_inputs/
```

## Reference RDS

Exact label column:

```text
fig4_display_label
```

Reference object:

```text
mac_provided_inputs/files/fig4_coad_read_reference_with_fig4_display_label.rds
```

This RDS is split for GitHub upload:

```text
mac_provided_inputs/files/fig4_coad_read_reference_with_fig4_display_label.rds.part000
mac_provided_inputs/files/fig4_coad_read_reference_with_fig4_display_label.rds.part001
mac_provided_inputs/files/fig4_coad_read_reference_with_fig4_display_label.rds.chunks.tsv
mac_provided_inputs/files/fig4_coad_read_reference_with_fig4_display_label.rds.sha256
```

Reconstruct on PC:

```cmd
copy /b mac_provided_inputs\files\fig4_coad_read_reference_with_fig4_display_label.rds.part000+mac_provided_inputs\files\fig4_coad_read_reference_with_fig4_display_label.rds.part001 mac_provided_inputs\files\fig4_coad_read_reference_with_fig4_display_label.rds
```

Or use any binary concatenation tool that preserves byte order.

Reference audit tables:

```text
mac_provided_inputs/tables/fig4_reference_label_counts.tsv
mac_provided_inputs/tables/fig4_reference_metadata_columns.tsv
```

The reference contains 8,374 COAD-READ cells and 12 labels in `fig4_display_label`.

## Supplementary Files

The required xlsx files are provided here:

```text
mac_provided_inputs/supplementary/cir-24-1255_table_s1_suppst1.xlsx
mac_provided_inputs/supplementary/cir-24-1255_table_s4_suppst4.xlsx
```

## Recommended Run

From:

```bash
cd bioinformatics-reproduction/projects/pan-cancer-mono-mac
```

Run:

```bash
FIG4_RUN_SCOPE="fig4c" \
FIG4_GEO_DIR="fig4_snatac/raw/GSE306459" \
FIG4_SUPP_DIR="desktop_exchange/uploads/20260703_fig4_gse306459_snatac/mac_provided_inputs/supplementary" \
FIG4_OUTPUT_DIR="fig4_snatac/outputs" \
FIG4_REFERENCE_RDS="desktop_exchange/uploads/20260703_fig4_gse306459_snatac/mac_provided_inputs/files/fig4_coad_read_reference_with_fig4_display_label.rds" \
FIG4_REFERENCE_LABEL_COLUMN="fig4_display_label" \
Rscript desktop_exchange/requests/20260703_fig4_gse306459_snatac/scripts/run_fig4_gse306459_snatac.R
```

For the current user request, `FIG4_RUN_SCOPE=fig4c` is sufficient. It should generate the full Fig. 4C ChIPSeeker annotation pie from differential peaks without extracting fragment files.
