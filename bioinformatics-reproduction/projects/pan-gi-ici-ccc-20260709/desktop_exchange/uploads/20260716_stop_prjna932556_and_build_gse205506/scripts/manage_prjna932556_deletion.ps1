param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('Prepare', 'DeleteHost', 'DeleteGuest', 'Audit')]
    [string]$Mode
)

$ErrorActionPreference = 'Stop'

$OutputDir = 'I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\uploads\20260716_stop_prjna932556_and_build_gse205506'
$ManifestPath = Join-Path $OutputDir 'prjna932556_deletion_manifest_before.tsv'
$AuditPath = Join-Path $OutputDir 'prjna932556_deletion_audit_after.tsv'
$SummaryPath = Join-Path $OutputDir 'prjna932556_deletion_summary.tsv'
$LocalOnlyRoot = 'I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only'
$Key = 'I:\codex-config\tools\ubuntu-vm\ssh\id_ed25519'
$KnownHosts = 'I:\codex-config\tools\ubuntu-vm\ssh\known_hosts_cellranger72_resume'

$Targets = @(
    [pscustomobject]@{
        absolute_path = "$LocalOnlyRoot\20260713_prjna932556_bamtofastq_feasibility"
        path_type = 'windows_directory'
        evidence_linking_path_to_prjna932556 = 'Exact path recorded by 20260713_prjna932556_bamtofastq_feasibility STATUS and BAM-to-FASTQ logs for SRR23490337/PRJNA932556.'
    },
    [pscustomobject]@{
        absolute_path = "$LocalOnlyRoot\20260713_prjna932556_ena_bam_pair_download"
        path_type = 'windows_directory'
        evidence_linking_path_to_prjna932556 = 'Exact path recorded by 20260713_prjna932556_ena_bam_pair_download_manifest.tsv and verified inventory for SRR23490337/PRJNA932556 BAM and BAI.'
    },
    [pscustomobject]@{
        absolute_path = "$LocalOnlyRoot\20260714_prjna932556_srr23490337_full_bamtofastq"
        path_type = 'windows_directory'
        evidence_linking_path_to_prjna932556 = 'Exact path recorded by 20260714_prjna932556_srr23490337_full_bamtofastq manifest and SHA256 inventory for SRR23490337/PRJNA932556 FASTQs.'
    },
    [pscustomobject]@{
        absolute_path = '/work/input/W6B_0_1_HGGWHDSX2'
        path_type = 'qemu_guest_directory'
        evidence_linking_path_to_prjna932556 = 'Exact Cell Ranger input path from the active R3_SRR23490337 command and live guest inventory; source FASTQs are audited as SRR23490337/PRJNA932556.'
    },
    [pscustomobject]@{
        absolute_path = '/work/run/R3_SRR23490337'
        path_type = 'qemu_guest_directory'
        evidence_linking_path_to_prjna932556 = 'Exact Cell Ranger pipestance path named by --id=R3_SRR23490337 and confirmed by live guest inventory.'
    },
    [pscustomobject]@{
        absolute_path = '/work/run/__R3_SRR23490337.mro'
        path_type = 'qemu_guest_file'
        evidence_linking_path_to_prjna932556 = 'Exact Martian run definition generated for Cell Ranger --id=R3_SRR23490337 and confirmed by live guest inventory.'
    }
)

function Get-WindowsPathBytes([string]$Path) {
    $item = Get-Item -LiteralPath $Path -Force
    if (-not $item.PSIsContainer) {
        return [int64]$item.Length
    }
    $measure = Get-ChildItem -LiteralPath $Path -File -Recurse -Force | Measure-Object Length -Sum
    return [int64]$measure.Sum
}

function Invoke-Guest([string]$Command) {
    $result = & ssh.exe -i $Key -o BatchMode=yes -o StrictHostKeyChecking=no -o "UserKnownHostsFile=$KnownHosts" -p 2222 ubuntu@127.0.0.1 $Command
    if ($LASTEXITCODE -ne 0) {
        throw "Guest command failed with exit code $LASTEXITCODE"
    }
    return $result
}

function Get-GuestPathBytes([string]$Path) {
    $line = Invoke-Guest "sudo du -sb -- '$Path'"
    return [int64](($line | Select-Object -First 1) -split "`t")[0]
}

function ConvertTo-Tsv([object[]]$Rows, [string]$Path) {
    $content = ($Rows | ConvertTo-Csv -Delimiter "`t" -NoTypeInformation) -join "`r`n"
    [System.IO.File]::WriteAllText($Path, $content + "`r`n", [System.Text.UTF8Encoding]::new($false))
}

if ($Mode -eq 'Prepare') {
    $rows = foreach ($target in $Targets) {
        $bytes = if ($target.path_type -like 'windows_*') {
            Get-WindowsPathBytes $target.absolute_path
        } else {
            Get-GuestPathBytes $target.absolute_path
        }
        [pscustomobject]@{
            absolute_path = $target.absolute_path
            path_type = $target.path_type
            size_bytes = $bytes
            evidence_linking_path_to_prjna932556 = $target.evidence_linking_path_to_prjna932556
            deletion_status = 'planned_exact_path_delete'
            error_message = ''
        }
    }
    ConvertTo-Tsv $rows $ManifestPath
    $rows | ConvertTo-Json -Depth 3
    exit 0
}

if (-not (Test-Path -LiteralPath $ManifestPath)) {
    throw "Deletion manifest is missing: $ManifestPath"
}

if ($Mode -eq 'DeleteHost') {
    $allowed = @($Targets | Where-Object { $_.path_type -like 'windows_*' } | ForEach-Object { $_.absolute_path })
    foreach ($path in $allowed) {
        $resolved = (Resolve-Path -LiteralPath $path).Path
        if ($resolved -ne $path -or -not $resolved.StartsWith($LocalOnlyRoot + '\', [System.StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing non-whitelisted path: $resolved"
        }
        Remove-Item -LiteralPath $resolved -Recurse -Force
    }
    exit 0
}

if ($Mode -eq 'DeleteGuest') {
    $command = "set -e; test -d /work/input/W6B_0_1_HGGWHDSX2; test -d /work/run/R3_SRR23490337; test -f /work/run/__R3_SRR23490337.mro; sudo rm -rf -- /work/input/W6B_0_1_HGGWHDSX2 /work/run/R3_SRR23490337; sudo rm -f -- /work/run/__R3_SRR23490337.mro; sync"
    Invoke-Guest $command | Out-Null
    exit 0
}

$before = Import-Csv -LiteralPath $ManifestPath -Delimiter "`t"
$auditRows = foreach ($row in $before) {
    $existsAfter = $false
    $remaining = [int64]0
    $errorMessage = ''
    try {
        if ($row.path_type -like 'windows_*') {
            $existsAfter = Test-Path -LiteralPath $row.absolute_path
            if ($existsAfter) { $remaining = Get-WindowsPathBytes $row.absolute_path }
        } else {
            $probe = Invoke-Guest "if sudo test -e '$($row.absolute_path)'; then echo EXISTS; sudo du -sb -- '$($row.absolute_path)'; else echo ABSENT; fi"
            $existsAfter = (($probe | Select-Object -First 1) -eq 'EXISTS')
            if ($existsAfter) { $remaining = [int64](($probe | Select-Object -Skip 1 -First 1) -split "`t")[0] }
        }
    } catch {
        $errorMessage = $_.Exception.Message
    }
    $beforeBytes = [int64]$row.size_bytes
    [pscustomobject]@{
        absolute_path = $row.absolute_path
        path_type = $row.path_type
        size_bytes_before = $beforeBytes
        exists_after = $existsAfter
        size_bytes_remaining = $remaining
        bytes_deleted = [math]::Max([int64]0, $beforeBytes - $remaining)
        deletion_status = if ($errorMessage) { 'audit_error' } elseif ($existsAfter) { 'remaining' } else { 'deleted_verified_absent' }
        error_message = $errorMessage
    }
}
ConvertTo-Tsv $auditRows $AuditPath
$summary = [pscustomobject]@{
    total_bytes_before = [int64](($auditRows | Measure-Object size_bytes_before -Sum).Sum)
    total_bytes_deleted = [int64](($auditRows | Measure-Object bytes_deleted -Sum).Sum)
    total_bytes_remaining = [int64](($auditRows | Measure-Object size_bytes_remaining -Sum).Sum)
    paths_total = $auditRows.Count
    paths_deleted_verified = @($auditRows | Where-Object deletion_status -eq 'deleted_verified_absent').Count
    paths_remaining = @($auditRows | Where-Object deletion_status -eq 'remaining').Count
    paths_audit_error = @($auditRows | Where-Object deletion_status -eq 'audit_error').Count
}
ConvertTo-Tsv @($summary) $SummaryPath
$auditRows | ConvertTo-Json -Depth 3
