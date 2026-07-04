$ErrorActionPreference = "Stop"

$repo = "I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac"
$requestDir = Join-Path $repo "desktop_exchange\requests\20260704_supp_fig_s6_spatial_signatures"
$defaultDataDir = Join-Path $repo "desktop_exchange\requests\20260704_supp_fig_s4_s5_spatial_featureplots\data"
$outDir = Join-Path $requestDir "outputs"
$tmpDir = "I:\codex-run-tmp\supp_fig_s6_spatial_signatures"
$logPath = Join-Path $requestDir "supp_fig_s6_spatial_signatures.console.log"
New-Item -ItemType Directory -Force -Path $tmpDir, $outDir | Out-Null
Set-Location -LiteralPath $repo

$env:R_LIBS_USER = "H:\R\library"
$env:R_LIBS = "H:\R\library"
$env:TMP = $tmpDir
$env:TEMP = $tmpDir
$env:TMPDIR = $tmpDir
$env:R_USER = "I:\codex-run-tmp"
$env:PATH = "H:\rtools44\usr\bin;H:\rtools44\x86_64-w64-mingw32.static.posix\bin;H:\R\bin;$env:PATH"
$env:S6_REQUEST_DIR = $requestDir
if (-not $env:S6_DATA_DIR) {
  $env:S6_DATA_DIR = $defaultDataDir
}
$env:S6_OUTPUT_DIR = $outDir

& "H:\R\bin\Rscript.exe" "desktop_exchange\requests\20260704_supp_fig_s6_spatial_signatures\scripts\run_supp_fig_s6_spatial_signatures.R" *>&1 |
  Tee-Object -FilePath $logPath
exit $LASTEXITCODE
