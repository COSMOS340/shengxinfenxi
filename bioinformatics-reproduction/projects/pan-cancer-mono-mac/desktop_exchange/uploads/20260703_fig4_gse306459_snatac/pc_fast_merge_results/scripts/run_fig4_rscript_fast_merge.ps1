$ErrorActionPreference = "Stop"

$repo = "I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac"
$requestDir = Join-Path $repo "desktop_exchange\requests\20260703_fig4_gse306459_snatac"
$logPath = Join-Path $repo "fig4_snatac\fig4_gse306459_snatac.fast_merge.console.log"
$tmpDir = "I:\codex-run-tmp\fig4_snatac_fast_merge"
$outDir = Join-Path $repo "fig4_snatac\outputs_fast_merge"
New-Item -ItemType Directory -Force -Path $tmpDir, $outDir | Out-Null
Set-Location -LiteralPath $repo

$env:R_LIBS_USER = "H:\R\library"
$env:R_LIBS = "H:\R\library"
$env:TMP = $tmpDir
$env:TEMP = $tmpDir
$env:TMPDIR = $tmpDir
$env:R_USER = "I:\codex-run-tmp"
$env:PATH = "H:\rtools44\usr\bin;H:\rtools44\x86_64-w64-mingw32.static.posix\bin;H:\R\bin;$env:PATH"
$env:FIG4_REQUEST_DIR = $requestDir
$env:FIG4_RUN_SCOPE = "fig4c"
$env:FIG4_GEO_DIR = Join-Path $repo "fig4_snatac\raw\GSE306459"
$env:FIG4_SUPP_DIR = Join-Path $repo "desktop_exchange\uploads\20260703_fig4_gse306459_snatac\mac_provided_inputs\supplementary"
$env:FIG4_OUTPUT_DIR = $outDir
$env:FIG4_REFERENCE_RDS = Join-Path $repo "desktop_exchange\uploads\20260703_fig4_gse306459_snatac\mac_provided_inputs\files\fig4_coad_read_reference_with_fig4_display_label.rds"
$env:FIG4_REFERENCE_LABEL_COLUMN = "fig4_display_label"

& "H:\R\bin\Rscript.exe" "I:\codex_tools\run_fig4_gse306459_snatac_fast_merge.R" *>&1 |
  Tee-Object -FilePath $logPath
exit $LASTEXITCODE
