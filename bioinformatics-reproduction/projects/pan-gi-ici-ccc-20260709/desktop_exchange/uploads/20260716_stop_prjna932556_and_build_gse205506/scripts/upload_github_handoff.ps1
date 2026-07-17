param(
    [Parameter(Mandatory = $true)]
    [string]$UploadDirectory,
    [Parameter(Mandatory = $true)]
    [string]$ExpectedHead,
    [Parameter(Mandatory = $true)]
    [string]$ExpectedRequestBlob
)

$ErrorActionPreference = 'Stop'
$repository = 'COSMOS340/shengxinfenxi'
$branch = 'pan-gi-ici-metadata-request-20260709'
$requestPath = 'bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/20260716_stop_prjna932556_and_build_gse205506_request.md'
$remoteBase = 'bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260716_stop_prjna932556_and_build_gse205506'
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

    $encodedRequestPath = [uri]::EscapeDataString($requestPath).Replace('%2F', '/')
    $requestFile = Invoke-RestMethod -Headers $headers -Uri "$apiBase/contents/$encodedRequestPath`?ref=$branch"
    if ($requestFile.sha -ne $ExpectedRequestBlob) {
        throw "Request file changed before upload. Expected $ExpectedRequestBlob but found $($requestFile.sha)"
    }

    $manifestPath = Join-Path $uploadDirectory 'output_manifest.tsv'
    $manifest = Import-Csv -LiteralPath $manifestPath -Delimiter "`t"
    $uploadRows = @($manifest | Where-Object delivery -eq 'github_upload')
    $uploadRows += [pscustomobject]@{
        relative_path = 'output_manifest.tsv'
        absolute_path = $manifestPath
    }

    $duplicatePaths = $uploadRows.relative_path | Group-Object | Where-Object Count -gt 1
    if ($duplicatePaths) {
        throw "Duplicate remote paths in upload list: $($duplicatePaths.Name -join ', ')"
    }

    $baseCommit = Invoke-RestMethod -Headers $headers -Uri "$apiBase/git/commits/$ExpectedHead"
    $treeEntries = foreach ($row in $uploadRows) {
        $absolutePath = $row.absolute_path
        if (-not (Test-Path -LiteralPath $absolutePath -PathType Leaf)) {
            throw "Upload file missing: $absolutePath"
        }
        $bytes = [IO.File]::ReadAllBytes($absolutePath)
        $blobBody = @{
            content = [Convert]::ToBase64String($bytes)
            encoding = 'base64'
        } | ConvertTo-Json -Compress
        $blob = Invoke-RestMethod -Method Post -Headers $headers -ContentType 'application/json' -Uri "$apiBase/git/blobs" -Body $blobBody
        [pscustomobject]@{
            path = "$remoteBase/$($row.relative_path)"
            mode = '100644'
            type = 'blob'
            sha = $blob.sha
        }
    }

    $treeBody = @{
        base_tree = $baseCommit.tree.sha
        tree = @($treeEntries)
    } | ConvertTo-Json -Depth 6 -Compress
    $tree = Invoke-RestMethod -Method Post -Headers $headers -ContentType 'application/json' -Uri "$apiBase/git/trees" -Body $treeBody

    $commitBody = @{
        message = 'Complete GSE205506 R Seurat build and marker review'
        tree = $tree.sha
        parents = @($ExpectedHead)
    } | ConvertTo-Json -Depth 4 -Compress
    $commit = Invoke-RestMethod -Method Post -Headers $headers -ContentType 'application/json' -Uri "$apiBase/git/commits" -Body $commitBody

    $updateBody = @{ sha = $commit.sha; force = $false } | ConvertTo-Json -Compress
    Invoke-RestMethod -Method Patch -Headers $headers -ContentType 'application/json' -Uri "$apiBase/git/refs/heads/$branch" -Body $updateBody | Out-Null

    $verifiedRef = Invoke-RestMethod -Headers $headers -Uri "$apiBase/git/ref/heads/$branch"
    if ($verifiedRef.object.sha -ne $commit.sha) {
        throw "Remote head verification failed after upload"
    }

    $readback = foreach ($relativePath in @('STATUS.md', 'MAC_REVIEW_REQUEST.md', 'output_manifest.tsv')) {
        $remotePath = "$remoteBase/$relativePath"
        $encodedPath = [uri]::EscapeDataString($remotePath).Replace('%2F', '/')
        $remoteFile = Invoke-RestMethod -Headers $headers -Uri "$apiBase/contents/$encodedPath`?ref=$branch"
        $expectedBlob = $treeEntries | Where-Object path -eq $remotePath
        if ($remoteFile.sha -ne $expectedBlob.sha) {
            throw "Read-back blob verification failed for $relativePath"
        }
        [pscustomobject]@{
            relative_path = $relativePath
            remote_blob_sha = $remoteFile.sha
            verified = $true
        }
    }

    $receiptPath = Join-Path $uploadDirectory 'github_upload_receipt.tsv'
    @(
        [pscustomobject]@{ metric = 'previous_head'; value = $ExpectedHead },
        [pscustomobject]@{ metric = 'new_head'; value = $commit.sha },
        [pscustomobject]@{ metric = 'request_blob'; value = $ExpectedRequestBlob },
        [pscustomobject]@{ metric = 'uploaded_files'; value = $uploadRows.Count },
        [pscustomobject]@{ metric = 'readback_verified_files'; value = $readback.Count }
    ) | Export-Csv -LiteralPath $receiptPath -Delimiter "`t" -NoTypeInformation -Encoding utf8

    [pscustomobject]@{
        previous_head = $ExpectedHead
        new_head = $commit.sha
        uploaded_files = $uploadRows.Count
        readback_verified_files = $readback.Count
        receipt_path = $receiptPath
    }
}
finally {
    $token = $null
    Remove-Variable secure, credential, headers -ErrorAction SilentlyContinue
}
