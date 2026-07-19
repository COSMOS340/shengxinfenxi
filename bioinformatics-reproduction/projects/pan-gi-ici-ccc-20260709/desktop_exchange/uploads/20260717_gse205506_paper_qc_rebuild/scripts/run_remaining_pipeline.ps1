param(
    [int]$FirstRoundProcessId = 25416
)

$ErrorActionPreference = 'Stop'
$uploadDirectory = 'I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\uploads\20260717_gse205506_paper_qc_rebuild'
$scriptDirectory = Join-Path $uploadDirectory 'scripts'
$tempDirectory = 'I:\temp\gse205506_paper_qc_20260717'
$rscript = 'H:\R\bin\Rscript.exe'

$env:TEMP = $tempDirectory
$env:TMP = $tempDirectory
$env:TMPDIR = $tempDirectory
$env:R_USER = 'I:\codex-config\r-user'
$env:R_LIBS_USER = 'H:\R\library'
New-Item -ItemType Directory -Force -Path $tempDirectory, $env:R_USER | Out-Null

function Write-ProgressLog {
    param([string]$Message)
    Write-Output "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] $Message"
}

function Invoke-RStage {
    param(
        [Parameter(Mandatory = $true)][string]$ScriptName,
        [Parameter(Mandatory = $true)][string]$LogStem,
        [Parameter(Mandatory = $true)][string]$SuccessText
    )

    $scriptPath = Join-Path $scriptDirectory $ScriptName
    $stdoutPath = Join-Path $uploadDirectory "$LogStem.stdout.log"
    $stderrPath = Join-Path $uploadDirectory "$LogStem.stderr.log"
    Write-ProgressLog "Starting $ScriptName"
    $process = Start-Process -FilePath $rscript -ArgumentList @('--vanilla', $scriptPath) `
        -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath `
        -WindowStyle Hidden -Wait -PassThru
    if ($process.ExitCode -ne 0) {
        throw "$ScriptName failed with exit code $($process.ExitCode); inspect $stderrPath"
    }
    if (-not (Select-String -LiteralPath $stdoutPath -SimpleMatch $SuccessText -Quiet)) {
        throw "$ScriptName exited without its required success marker: $SuccessText"
    }
    if (Select-String -LiteralPath $stderrPath -SimpleMatch 'Execution halted' -Quiet) {
        throw "$ScriptName stderr contains Execution halted"
    }
    Write-ProgressLog "Completed $ScriptName"
}

$firstRoundProcess = Get-Process -Id $FirstRoundProcessId -ErrorAction SilentlyContinue
if ($firstRoundProcess) {
    Write-ProgressLog "Waiting for first-round R process $FirstRoundProcessId"
    Wait-Process -Id $FirstRoundProcessId
}
$firstStdout = Join-Path $uploadDirectory '02_first_round_rpca_broad_mt_qc.stdout.log'
$firstStderr = Join-Path $uploadDirectory '02_first_round_rpca_broad_mt_qc.stderr.log'
if (-not (Select-String -LiteralPath $firstStdout -SimpleMatch 'First-round RPCA broad classification and compartment mitochondrial QC complete' -Quiet)) {
    throw 'First-round stage exited without its required success marker'
}
if (Select-String -LiteralPath $firstStderr -SimpleMatch 'Execution halted' -Quiet) {
    throw 'First-round stage stderr contains Execution halted'
}
Write-ProgressLog 'First-round stage verified'

Invoke-RStage `
    -ScriptName '03_rebuild_final_rpca_atlas.R' `
    -LogStem '03_rebuild_final_rpca_atlas' `
    -SuccessText 'Final post-mitochondrial-QC RPCA atlas complete'

Invoke-RStage `
    -ScriptName '04_recluster_compartments.R' `
    -LogStem '04_recluster_compartments' `
    -SuccessText 'All six compartment objects complete; no whole cluster deleted'

Invoke-RStage `
    -ScriptName '05_render_paper_qc_figures.R' `
    -LogStem '05_render_paper_qc_figures' `
    -SuccessText 'All required paper-QC figures rendered'

Write-ProgressLog 'Automated analysis stages complete; manual figure inspection is next'
