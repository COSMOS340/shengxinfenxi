$ErrorActionPreference = "Stop"

$repo = "I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac"
$tmpDir = "I:\codex-run-tmp\fig4a_formal_snatac"
$outputDir = Join-Path $repo "desktop_exchange\uploads\20260704_fig4a_gse306459_formal_snatac\pc_formal_fig4a_results"
$requestDir = Join-Path $repo "desktop_exchange\requests\20260704_fig4a_gse306459_formal_snatac"
$logPath = Join-Path $repo "desktop_exchange\uploads\20260704_fig4a_gse306459_formal_snatac\pc_formal_fig4a_console.log"

New-Item -ItemType Directory -Force -Path $tmpDir, $outputDir, (Split-Path $logPath) | Out-Null
Set-Location -LiteralPath $repo

$env:R_LIBS_USER = "H:\R\library"
$env:R_LIBS = "H:\R\library"
$env:TMP = $tmpDir
$env:TEMP = $tmpDir
$env:TMPDIR = $tmpDir
$env:R_USER = "I:\codex-run-tmp"
$env:PATH = "H:\rtools44\usr\bin;H:\rtools44\x86_64-w64-mingw32.static.posix\bin;H:\R\bin;$env:PATH"

$env:FIG4_REPO_DIR = $repo
$env:FIG4_REQUEST_DIR = $requestDir
$env:FIG4_GEO_DIR = Join-Path $repo "fig4_snatac\raw\GSE306459"
$env:FIG4_SUPP_DIR = Join-Path $repo "desktop_exchange\uploads\20260703_fig4_gse306459_snatac\mac_provided_inputs\supplementary"
$env:FIG4_OUTPUT_DIR = $outputDir
$env:FIG4_REFERENCE_RDS = Join-Path $repo "desktop_exchange\uploads\20260703_fig4_gse306459_snatac\mac_provided_inputs\files\fig4_coad_read_reference_with_fig4_display_label.rds"
$env:FIG4_REFERENCE_LABEL_COLUMN = "fig4_display_label"
$env:FIG4_RUN_SCOPE = "fig4a_formal"
$env:FIG4_NPCS = "30"
$env:FIG4_RESOLUTION = "0.4"

& "H:\R\bin\Rscript.exe" "I:\codex_tools\run_fig4a_formal_snatac.R" *>&1 |
  Tee-Object -FilePath $logPath
exit $LASTEXITCODE
