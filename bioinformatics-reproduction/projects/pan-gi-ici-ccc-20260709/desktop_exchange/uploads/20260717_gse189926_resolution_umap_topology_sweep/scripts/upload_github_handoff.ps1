param(
    [Parameter(Mandatory = $true)][string]$UploadDirectory,
    [Parameter(Mandatory = $true)][string]$ExpectedHead,
    [Parameter(Mandatory = $true)][string]$ExpectedRequestBlob,
    [Parameter(Mandatory = $true)][string]$ExpectedRequestManifestBlob
)

$ErrorActionPreference = 'Stop'
$repository = 'COSMOS340/shengxinfenxi'
$branch = 'pan-gi-ici-metadata-request-20260709'
$requestPath = 'bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/20260717_gse189926_resolution_umap_topology_sweep_request.md'
$requestManifestPath = 'bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/20260717_gse189926_resolution_umap_topology_sweep_manifest.tsv'
$remoteBase = 'bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260717_gse189926_resolution_umap_topology_sweep'
$uploadDirectory = (Resolve-Path -LiteralPath $UploadDirectory).Path

$secure = Get-Content -LiteralPath 'I:\codex-config\github_pat.securestring' | ConvertTo-SecureString
$credential = [pscredential]::new('github', $secure)
$token = $credential.GetNetworkCredential().Password
$headers = @{
    Authorization = "Bearer $token"
    Accept = 'application/vnd.github+json'
    'X-GitHub-Api-Version' = '2022-11-28'
    'User-Agent' = 'Codex'
}
$apiBase = "https://api.github.com/repos/$repository"

try {
    $ref = Invoke-RestMethod -Headers $headers -Uri "$apiBase/git/ref/heads/$branch"
    if ($ref.object.sha -ne $ExpectedHead) {
        throw "Remote branch changed before upload. Expected $ExpectedHead but found $($ref.object.sha)"
    }
    foreach ($request in @(
        @{ Path = $requestPath; Sha = $ExpectedRequestBlob },
        @{ Path = $requestManifestPath; Sha = $ExpectedRequestManifestBlob }
    )) {
        $encodedPath = [uri]::EscapeDataString($request.Path).Replace('%2F', '/')
        $remoteFile = Invoke-RestMethod -Headers $headers -Uri "$apiBase/contents/$encodedPath`?ref=$branch"
        if ($remoteFile.sha -ne $request.Sha) {
            throw "Mac request changed before upload: $($request.Path)"
        }
    }

    $manifestPath = Join-Path $uploadDirectory 'output_manifest.tsv'
    $manifest = Import-Csv -LiteralPath $manifestPath -Delimiter "`t"
    $uploadRows = @($manifest | Where-Object delivery_status -eq 'github_upload')
    $uploadRows += [pscustomobject]@{ relative_path = 'output_manifest.tsv'; absolute_path = $manifestPath }
    $duplicatePaths = $uploadRows.relative_path | Group-Object | Where-Object Count -gt 1
    if ($duplicatePaths) { throw "Duplicate remote paths: $($duplicatePaths.Name -join ', ')" }

    $baseCommit = Invoke-RestMethod -Headers $headers -Uri "$apiBase/git/commits/$ExpectedHead"
    $treeEntries = foreach ($row in $uploadRows) {
        if (-not (Test-Path -LiteralPath $row.absolute_path -PathType Leaf)) { throw "Upload file missing: $($row.absolute_path)" }
        $blobBody = @{
            content = [Convert]::ToBase64String([IO.File]::ReadAllBytes($row.absolute_path))
            encoding = 'base64'
        } | ConvertTo-Json -Compress
        $blob = Invoke-RestMethod -Method Post -Headers $headers -ContentType 'application/json' -Uri "$apiBase/git/blobs" -Body $blobBody
        [pscustomobject]@{ path = "$remoteBase/$($row.relative_path)"; mode = '100644'; type = 'blob'; sha = $blob.sha }
    }

    $treeBody = @{ base_tree = $baseCommit.tree.sha; tree = @($treeEntries) } | ConvertTo-Json -Depth 6 -Compress
    $tree = Invoke-RestMethod -Method Post -Headers $headers -ContentType 'application/json' -Uri "$apiBase/git/trees" -Body $treeBody
    $commitBody = @{
        message = 'Complete GSE189926 resolution and UMAP topology sweep'
        tree = $tree.sha
        parents = @($ExpectedHead)
    } | ConvertTo-Json -Depth 4 -Compress
    $commit = Invoke-RestMethod -Method Post -Headers $headers -ContentType 'application/json' -Uri "$apiBase/git/commits" -Body $commitBody
    $updateBody = @{ sha = $commit.sha; force = $false } | ConvertTo-Json -Compress
    Invoke-RestMethod -Method Patch -Headers $headers -ContentType 'application/json' -Uri "$apiBase/git/refs/heads/$branch" -Body $updateBody | Out-Null

    $verifiedRef = Invoke-RestMethod -Headers $headers -Uri "$apiBase/git/ref/heads/$branch"
    if ($verifiedRef.object.sha -ne $commit.sha) { throw 'Remote head verification failed' }
    foreach ($relativePath in @('STATUS.md', 'MAC_REVIEW_REQUEST.md', 'output_manifest.tsv')) {
        $remotePath = "$remoteBase/$relativePath"
        $encodedPath = [uri]::EscapeDataString($remotePath).Replace('%2F', '/')
        $remoteFile = Invoke-RestMethod -Headers $headers -Uri "$apiBase/contents/$encodedPath`?ref=$branch"
        $expectedBlob = $treeEntries | Where-Object path -eq $remotePath
        if ($remoteFile.sha -ne $expectedBlob.sha) { throw "Read-back failed for $relativePath" }
    }

    $receiptPath = Join-Path $uploadDirectory 'github_upload_receipt.tsv'
    @(
        [pscustomobject]@{ metric = 'previous_head'; value = $ExpectedHead },
        [pscustomobject]@{ metric = 'new_head'; value = $commit.sha },
        [pscustomobject]@{ metric = 'request_blob'; value = $ExpectedRequestBlob },
        [pscustomobject]@{ metric = 'request_manifest_blob'; value = $ExpectedRequestManifestBlob },
        [pscustomobject]@{ metric = 'uploaded_files'; value = $uploadRows.Count },
        [pscustomobject]@{ metric = 'readback_verified_files'; value = 3 }
    ) | Export-Csv -LiteralPath $receiptPath -Delimiter "`t" -NoTypeInformation -Encoding utf8
    [pscustomobject]@{ previous_head = $ExpectedHead; new_head = $commit.sha; uploaded_files = $uploadRows.Count; receipt_path = $receiptPath }
}
finally {
    $token = $null
    Remove-Variable secure, credential, headers -ErrorAction SilentlyContinue
}
