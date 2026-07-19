param(
    [Parameter(Mandatory = $true)]
    [string]$UploadDirectory
)

$ErrorActionPreference = 'Stop'
$uploadDirectory = (Resolve-Path -LiteralPath $UploadDirectory).Path
$manifestPath = Join-Path $uploadDirectory 'output_manifest.tsv'
$excludedNames = @('output_manifest.tsv', 'github_upload_receipt.tsv')

$uploadFiles = Get-ChildItem -LiteralPath $uploadDirectory -Recurse -File |
    Where-Object { $_.Name -notin $excludedNames } |
    Sort-Object FullName

$rows = foreach ($item in $uploadFiles) {
    $relativePath = $item.FullName.Substring($uploadDirectory.Length).TrimStart('\').Replace('\', '/')
    $notes = switch -Regex ($relativePath) {
        '^figures/.+\.png$' { 'visually inspected 300-dpi raster figure'; break }
        '^figures/.+\.pdf$' { 'vector-container PDF with rasterized UMAP points where applicable'; break }
        '^scripts/' { 'reproducibility script'; break }
        '\.(stdout|stderr)\.log$' { 'analysis execution log'; break }
        '^sessionInfo_.+\.txt$' { 'R session information'; break }
        '^STATUS\.md$' { 'completion status'; break }
        '^MAC_REVIEW_REQUEST\.md$' { 'Mac review checklist'; break }
        default { 'lightweight reproducibility output' }
    }
    [pscustomobject]@{
        relative_path = $relativePath
        absolute_path = $item.FullName
        size_bytes = $item.Length
        sha256 = (Get-FileHash -LiteralPath $item.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
        delivery = 'github_upload'
        notes = $notes
    }
}

$localInventoryPath = Join-Path $uploadDirectory 'gse205506_paper_qc_local_object_inventory.tsv'
if (-not (Test-Path -LiteralPath $localInventoryPath -PathType Leaf)) {
    throw "Local object inventory is missing: $localInventoryPath"
}
$localInventory = Import-Csv -LiteralPath $localInventoryPath -Delimiter "`t"
foreach ($entry in $localInventory) {
    $rows += [pscustomobject]@{
        relative_path = ''
        absolute_path = $entry.absolute_path
        size_bytes = [int64]$entry.size_bytes
        sha256 = $entry.sha256.ToLowerInvariant()
        delivery = 'local_only'
        notes = "$($entry.object_type): $($entry.status); $($entry.notes)"
    }
}

$duplicatePaths = @($rows | Where-Object delivery -eq 'github_upload' | Group-Object relative_path | Where-Object Count -gt 1)
if ($duplicatePaths.Count -gt 0) {
    throw "Duplicate GitHub paths in manifest: $($duplicatePaths.Name -join ', ')"
}
$oversized = @($rows | Where-Object { $_.delivery -eq 'github_upload' -and $_.size_bytes -ge 100MB })
if ($oversized.Count -gt 0) {
    throw "GitHub upload contains files at or above 100 MB: $($oversized.relative_path -join ', ')"
}

$rows | Export-Csv -LiteralPath $manifestPath -Delimiter "`t" -NoTypeInformation -Encoding utf8
$readback = Import-Csv -LiteralPath $manifestPath -Delimiter "`t"
if ($readback.Count -ne $rows.Count) {
    throw 'Manifest row-count verification failed'
}

[pscustomobject]@{
    manifest_path = $manifestPath
    github_files = @($rows | Where-Object delivery -eq 'github_upload').Count
    local_only_files = @($rows | Where-Object delivery -eq 'local_only').Count
    manifest_rows = $rows.Count
}

