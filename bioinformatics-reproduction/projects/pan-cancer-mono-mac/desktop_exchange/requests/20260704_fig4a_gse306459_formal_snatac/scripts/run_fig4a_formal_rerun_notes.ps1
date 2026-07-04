$ErrorActionPreference = "Stop"

# Adjust these paths if the PC checkout root differs.
$env:FIG4_GEO_DIR = "fig4_snatac/raw/GSE306459"
$env:FIG4_SUPP_DIR = "desktop_exchange/uploads/20260703_fig4_gse306459_snatac/mac_provided_inputs/supplementary"
$env:FIG4_OUTPUT_DIR = "desktop_exchange/uploads/20260704_fig4a_gse306459_formal_snatac/pc_formal_fig4a_results"
$env:FIG4_REFERENCE_RDS = "desktop_exchange/uploads/20260703_fig4_gse306459_snatac/mac_provided_inputs/files/fig4_coad_read_reference_with_fig4_display_label.rds"
$env:FIG4_REFERENCE_LABEL_COLUMN = "fig4_display_label"
$env:FIG4_RUN_SCOPE = "fig4a_formal"

New-Item -ItemType Directory -Force -Path $env:FIG4_OUTPUT_DIR | Out-Null

Write-Host "Formal Fig4A rerun request."
Write-Host "Do not use the fast-merge UMAP as final Fig4A."
Write-Host "Run a fragments-backed workflow with snATAC QC before LSI/UMAP."
Write-Host "R packages are not uploaded by Mac; install missing packages locally on PC."
Write-Host ""
Write-Host "Required upload directory:"
Write-Host $env:FIG4_OUTPUT_DIR
Write-Host ""
Write-Host "Read the request README before running:"
Write-Host "desktop_exchange/requests/20260704_fig4a_gse306459_formal_snatac/README.md"
