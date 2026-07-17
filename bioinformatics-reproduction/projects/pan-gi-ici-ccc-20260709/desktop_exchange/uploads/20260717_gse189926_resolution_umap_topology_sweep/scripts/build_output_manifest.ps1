param([Parameter(Mandatory = $true)][string]$UploadDirectory)

$ErrorActionPreference = 'Stop'
$uploadDirectory = (Resolve-Path -LiteralPath $UploadDirectory).Path
$required = @(
    'STATUS.md', 'MAC_REVIEW_REQUEST.md',
    'gse189926_sweep_input_object_audit.tsv',
    'gse189926_resolution_sweep_summary.tsv',
    'gse189926_resolution_cluster_counts.tsv',
    'gse189926_resolution_sample_patient_dominance.tsv',
    'gse189926_resolution_transition_tables.tsv',
    'gse189926_resolution_top10_markers.tsv',
    'gse189926_resolution_top50_markers.tsv',
    'gse189926_umap_route_parameters.tsv',
    'gse189926_umap_route_mixing_metrics.tsv',
    'gse189926_umap_local_linearity_summary.tsv',
    'gse189926_umap_local_linearity_by_cluster.tsv',
    'gse189926_filament_region_cell_audit.tsv.gz',
    'gse189926_filament_region_review.tsv',
    'gse189926_resolution_selection_evidence.tsv',
    'gse189926_umap_selection_evidence.tsv'
)
foreach ($relativePath in $required) {
    if (-not (Test-Path -LiteralPath (Join-Path $uploadDirectory $relativePath) -PathType Leaf)) {
        throw "Required output is missing: $relativePath"
    }
}

$topFiles = Get-ChildItem -LiteralPath $uploadDirectory -File | Where-Object {
    $_.Name -notin @('output_manifest.tsv', 'github_upload_receipt.tsv')
}
$childFiles = @(
    $(Get-ChildItem -LiteralPath (Join-Path $uploadDirectory 'scripts') -File)
    $(Get-ChildItem -LiteralPath (Join-Path $uploadDirectory 'figures') -File)
)
$githubItems = @($topFiles) + @($childFiles)

$rows = foreach ($item in $githubItems) {
    $relativePath = $item.FullName.Substring($uploadDirectory.Length + 1).Replace('\', '/')
    [pscustomobject]@{
        relative_path = $relativePath
        absolute_path = $item.FullName
        size_bytes = $item.Length
        sha256 = (Get-FileHash -LiteralPath $item.FullName -Algorithm SHA256).Hash
        delivery_status = 'github_upload'
        notes = 'GSE189926 topology sweep lightweight output'
    }
}

$localOnly = @(
    @{ Path = 'I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\04_objects\20260717_gse189926_resolution_umap_topology_sweep\gse189926_resolution_umap_topology_sweep.rds'; Note = 'final sweep object; local only; original object not overwritten' },
    @{ Path = 'I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\04_objects\20260717_gse189926_resolution_umap_topology_sweep\gse189926_sweep_pre_umap_checkpoint.rds'; Note = 'pre-UMAP checkpoint; local only' },
    @{ Path = 'I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\04_objects\20260717_gse189926_resolution_umap_topology_sweep\gse189926_sweep_post_umap_checkpoint.rds'; Note = 'post-UMAP checkpoint; local only' }
)
foreach ($entry in $localOnly) {
    $item = Get-Item -LiteralPath $entry.Path
    $rows += [pscustomobject]@{
        relative_path = ''
        absolute_path = $item.FullName
        size_bytes = $item.Length
        sha256 = (Get-FileHash -LiteralPath $item.FullName -Algorithm SHA256).Hash
        delivery_status = 'local_only'
        notes = $entry.Note
    }
}

$manifestPath = Join-Path $uploadDirectory 'output_manifest.tsv'
$rows | Sort-Object delivery_status, relative_path, absolute_path |
    Export-Csv -LiteralPath $manifestPath -Delimiter "`t" -NoTypeInformation -Encoding utf8

$manifest = Import-Csv -LiteralPath $manifestPath -Delimiter "`t"
if (($manifest | Where-Object delivery_status -eq 'github_upload').Count -ne $githubItems.Count) {
    throw 'Manifest GitHub file count verification failed'
}
if (($manifest | Where-Object { -not $_.sha256 }).Count -ne 0) {
    throw 'Manifest contains a blank SHA-256 value'
}

[pscustomobject]@{
    manifest_path = $manifestPath
    github_files = $githubItems.Count
    local_only_files = $localOnly.Count
    manifest_rows = $manifest.Count
}
