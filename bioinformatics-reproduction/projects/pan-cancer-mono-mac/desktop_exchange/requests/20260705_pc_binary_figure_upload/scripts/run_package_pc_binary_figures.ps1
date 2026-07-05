$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$packager = Join-Path $scriptDir "package_pc_binary_figures.py"

if ($env:PC_BINARY_FIGURE_PYTHON) {
    $python = $env:PC_BINARY_FIGURE_PYTHON
} elseif (Get-Command python -ErrorAction SilentlyContinue) {
    $python = "python"
} elseif (Get-Command py -ErrorAction SilentlyContinue) {
    $python = "py"
} else {
    throw "No Python executable found. Set PC_BINARY_FIGURE_PYTHON to a Python executable path."
}

& $python $packager
if ($LASTEXITCODE -ne 0) {
    throw "package_pc_binary_figures.py failed with exit code $LASTEXITCODE"
}

