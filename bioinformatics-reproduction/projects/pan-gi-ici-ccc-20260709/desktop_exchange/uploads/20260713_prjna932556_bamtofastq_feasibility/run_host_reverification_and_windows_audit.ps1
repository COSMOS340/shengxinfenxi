$ErrorActionPreference = 'Stop'

$exchange = 'I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange'
$inputDir = Join-Path $exchange 'local_only\20260713_prjna932556_ena_bam_pair_download'
$uploadDir = Join-Path $exchange 'uploads\20260713_prjna932556_bamtofastq_feasibility'
$bam = Join-Path $inputDir 'R3_possorted_genome_bam.bam.1'
$bai = Join-Path $inputDir 'R3_possorted_genome_bam.bam.1.bai'
$auditPath = Join-Path $uploadDir 'runtime_tool_audit.txt'
$reverificationPath = Join-Path $uploadDir 'input_bam_bai_reverification.tsv'
$utf8 = New-Object System.Text.UTF8Encoding($false)

New-Item -ItemType Directory -Path $uploadDir -Force | Out-Null

function Invoke-WslAudit {
    param(
        [Parameter(Mandatory = $true)][string]$Label,
        [Parameter(Mandatory = $true)][string]$Arguments
    )

    $safeLabel = $Label -replace '[^A-Za-z0-9_.-]', '_'
    $tempPath = Join-Path $uploadDir ('.' + $safeLabel + '.utf16.txt')
    $commandLine = 'chcp 65001 > nul & wsl.exe ' + $Arguments + ' > "' + $tempPath + '" 2>&1'
    & cmd.exe /d /c $commandLine | Out-Null
    $exitCode = $LASTEXITCODE
    $output = if (Test-Path -LiteralPath $tempPath) {
        [System.IO.File]::ReadAllText($tempPath, [System.Text.Encoding]::Unicode).TrimEnd()
    } else {
        ''
    }
    if (Test-Path -LiteralPath $tempPath) {
        Remove-Item -LiteralPath $tempPath -Force
    }
    if ($exitCode -ne 0) {
        $output = 'Command returned a nonzero exit code; localized usage text omitted. Windows feature states are recorded below.'
    }

    return @(
        "COMMAND: wsl.exe $Arguments",
        "EXIT_CODE: $exitCode",
        'OUTPUT:',
        $(if ($output) { $output } else { '<no output>' }),
        ''
    )
}

$audit = New-Object System.Collections.Generic.List[string]
$audit.Add('SRR23490337 bamtofastq feasibility runtime and tool audit')
$audit.Add("AUDIT_START: $((Get-Date).ToString('o'))")
$audit.Add("HOST_OS: $([System.Environment]::OSVersion.VersionString)")
$audit.Add("POWERSHELL: $($PSVersionTable.PSVersion)")
$audit.Add('')

$audit.Add('COMMAND: chcp 65001')
$chcpOutput = (& cmd.exe /d /c 'chcp 65001' 2>&1 | Out-String).TrimEnd()
$audit.Add("EXIT_CODE: $LASTEXITCODE")
$audit.Add('OUTPUT:')
$audit.Add($(if ($chcpOutput) { $chcpOutput } else { '<no output>' }))
$audit.Add('')

foreach ($line in (Invoke-WslAudit -Label 'wsl_status' -Arguments '--status')) {
    $audit.Add([string]$line)
}
foreach ($line in (Invoke-WslAudit -Label 'wsl_list_verbose' -Arguments '-l -v')) {
    $audit.Add([string]$line)
}

$audit.Add('COMMAND: docker.exe --version')
$docker = Get-Command docker.exe -ErrorAction SilentlyContinue
if ($null -eq $docker) {
    $audit.Add('EXIT_CODE: 127')
    $audit.Add('OUTPUT:')
    $audit.Add('docker.exe not found on PATH')
} else {
    $dockerOutput = (& $docker.Source --version 2>&1 | Out-String).TrimEnd()
    $audit.Add("EXIT_CODE: $LASTEXITCODE")
    $audit.Add('OUTPUT:')
    $audit.Add($(if ($dockerOutput) { $dockerOutput } else { '<no output>' }))
}
$audit.Add('')

foreach ($featureName in @('Microsoft-Windows-Subsystem-Linux', 'VirtualMachinePlatform')) {
    $feature = Get-WindowsOptionalFeature -Online -FeatureName $featureName
    $audit.Add("WINDOWS_FEATURE: $featureName")
    $audit.Add("STATE: $($feature.State)")
}
$audit.Add('')

foreach ($toolName in @('bash', 'samtools', 'bamtofastq')) {
    $tool = Get-Command $toolName -ErrorAction SilentlyContinue
    if ($null -eq $tool) {
        $audit.Add("WINDOWS_PATH_TOOL: $toolName = NOT_FOUND")
    } else {
        $audit.Add("WINDOWS_PATH_TOOL: $toolName = $($tool.Source)")
    }
}
$audit.Add("AUDIT_END: $((Get-Date).ToString('o'))")
[System.IO.File]::WriteAllLines($auditPath, $audit, $utf8)

$expected = @(
    [pscustomobject]@{ Role = 'BAM'; Path = $bam; ExpectedBytes = [int64]21494429186; ExpectedMD5 = '056a5c996212a64af708d48c5b585e48' },
    [pscustomobject]@{ Role = 'BAI'; Path = $bai; ExpectedBytes = [int64]10390008; ExpectedMD5 = 'f391d287338dbf9d3856630338957bb8' }
)

$rows = foreach ($item in $expected) {
    $file = Get-Item -LiteralPath $item.Path
    $actualMD5 = (Get-FileHash -LiteralPath $item.Path -Algorithm MD5).Hash.ToLowerInvariant()
    [pscustomobject]@{
        role = $item.Role
        pc_local_path = $item.Path
        expected_bytes = $item.ExpectedBytes
        actual_bytes = $file.Length
        bytes_match = ($file.Length -eq $item.ExpectedBytes)
        expected_md5 = $item.ExpectedMD5
        actual_md5 = $actualMD5
        md5_match = ($actualMD5 -eq $item.ExpectedMD5)
        checked_at = (Get-Date).ToString('o')
    }
}

$header = "role`tpc_local_path`texpected_bytes`tactual_bytes`tbytes_match`texpected_md5`tactual_md5`tmd5_match`tchecked_at"
$lines = New-Object System.Collections.Generic.List[string]
$lines.Add($header)
foreach ($row in $rows) {
    $lines.Add((@(
        $row.role,
        $row.pc_local_path,
        $row.expected_bytes,
        $row.actual_bytes,
        $row.bytes_match.ToString().ToLowerInvariant(),
        $row.expected_md5,
        $row.actual_md5,
        $row.md5_match.ToString().ToLowerInvariant(),
        $row.checked_at
    ) -join "`t"))
}
[System.IO.File]::WriteAllLines($reverificationPath, $lines, $utf8)

if (($rows | Where-Object { -not $_.bytes_match -or -not $_.md5_match }).Count -gt 0) {
    throw 'BAM/BAI reverification failed. See input_bam_bai_reverification.tsv.'
}

Write-Output "PHASE1_OK`t$reverificationPath`t$auditPath"
