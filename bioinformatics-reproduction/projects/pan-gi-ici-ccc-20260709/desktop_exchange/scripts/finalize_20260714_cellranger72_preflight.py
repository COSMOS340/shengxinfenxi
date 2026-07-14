from __future__ import annotations

import csv
import ctypes
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from urllib.error import HTTPError, URLError
from urllib.request import Request, build_opener, HTTPRedirectHandler


LOCAL_ROOT = Path(r"I:\shengxinfenxi")
EXCHANGE_DIR = (
    LOCAL_ROOT
    / "bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange"
)
REQUEST_NAME = "20260714_prjna932556_cellranger72_preflight"
UPLOAD_DIR = EXCHANGE_DIR / "uploads" / REQUEST_NAME
FASTQ_ROOT = (
    EXCHANGE_DIR
    / "local_only/20260714_prjna932556_srr23490337_full_bamtofastq"
)
GUEST_AUDIT = Path(
    r"I:\codex-config\tools\ubuntu-vm\runs\20260714_cellranger72_preflight_audit\guest_resource_audit.txt"
)
QEMU_EXE = Path(r"I:\codex-config\tools\qemu\20260501\qemu-system-x86_64.exe")
QEMU_SERIAL = Path(
    r"I:\codex-config\tools\ubuntu-vm\runs\20260714_cellranger72_preflight_audit\serial.log"
)
PREVIOUS_UPLOAD = (
    EXCHANGE_DIR
    / "uploads/20260714_prjna932556_srr23490337_full_bamtofastq"
)
CELLRANGER_DOWNLOAD_PAGE = "https://www.10xgenomics.com/support/software/cell-ranger/downloads"
REFDATA_URL = "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


CHECKED_AT = now_iso()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [
            {key.lstrip("\ufeff").strip('"'): value for key, value in row.items()}
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def run_command(args: list[str], timeout: int = 60) -> tuple[int, str]:
    try:
        result = subprocess.run(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    except (FileNotFoundError, OSError) as exc:
        return 127, str(exc)
    raw = result.stdout
    if raw.count(b"\x00") > len(raw) // 8:
        text = raw.decode("utf-16le", errors="replace")
    else:
        text = raw.decode("utf-8", errors="replace")
    return result.returncode, text.strip()


def powershell_json(script: str) -> object:
    code, text = run_command(
        [
            "powershell",
            "-NoProfile",
            "-ExecutionPolicy",
            "Bypass",
            "-Command",
            script,
        ],
        timeout=60,
    )
    if code != 0:
        raise RuntimeError(f"PowerShell command failed with exit {code}: {text}")
    return json.loads(text)


def cpuid_features() -> dict[str, object]:
    if os.name != "nt":
        return {"host_avx_supported": "not_windows", "host_avx2_supported": "not_windows"}
    kernel32 = ctypes.windll.kernel32
    mem_commit = 0x1000
    mem_reserve = 0x2000
    page_execute_readwrite = 0x40
    cpuid_code = bytes.fromhex(
        "53 89 C8 89 D1 0F A2 41 89 00 41 89 58 04 41 89 48 08 41 89 50 0C 5B C3"
    )
    kernel32.VirtualAlloc.restype = ctypes.c_void_p
    addr = kernel32.VirtualAlloc(
        None,
        len(cpuid_code),
        mem_commit | mem_reserve,
        page_execute_readwrite,
    )
    ctypes.memmove(addr, cpuid_code, len(cpuid_code))
    func = ctypes.WINFUNCTYPE(None, ctypes.c_uint32, ctypes.c_uint32, ctypes.c_void_p)(addr)

    def cpuid(leaf: int, subleaf: int = 0) -> tuple[int, int, int, int]:
        out = (ctypes.c_uint32 * 4)()
        func(leaf, subleaf, ctypes.addressof(out))
        return int(out[0]), int(out[1]), int(out[2]), int(out[3])

    _, _, ecx1, _ = cpuid(1, 0)
    _, ebx7, _, _ = cpuid(7, 0)
    xsave = bool(ecx1 & (1 << 26))
    osxsave = bool(ecx1 & (1 << 27))
    avx_cpu = bool(ecx1 & (1 << 28))
    avx2_cpu = bool(ebx7 & (1 << 5))
    xcr0 = None
    os_avx_state = False
    if osxsave:
        xgetbv_code = bytes.fromhex("0F 01 D0 48 C1 E2 20 48 09 D0 C3")
        xaddr = kernel32.VirtualAlloc(
            None,
            len(xgetbv_code),
            mem_commit | mem_reserve,
            page_execute_readwrite,
        )
        ctypes.memmove(xaddr, xgetbv_code, len(xgetbv_code))
        xfunc = ctypes.WINFUNCTYPE(ctypes.c_uint64, ctypes.c_uint32)(xaddr)
        xcr0 = int(xfunc(0))
        os_avx_state = (xcr0 & 0x6) == 0x6
    return {
        "host_xsave_supported": str(xsave),
        "host_osxsave_supported": str(osxsave),
        "host_xcr0_hex": "" if xcr0 is None else hex(xcr0),
        "host_avx_supported": str(avx_cpu and osxsave and os_avx_state),
        "host_avx2_supported": str(avx2_cpu and avx_cpu and osxsave and os_avx_state),
    }


def fastq_role(path: Path) -> tuple[str, str, str] | None:
    match = re.match(r"^(.+_L\d{3})_([RI][12])_(\d{3})\.fastq\.gz$", path.name)
    if not match:
        return None
    return match.group(1), match.group(2), match.group(3)


def first_100_lengths(path: Path) -> dict[int, int]:
    lengths: dict[int, int] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        for index in range(100):
            header = handle.readline()
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not (header and sequence and plus and quality):
                break
            length = len(sequence.rstrip("\r\n"))
            lengths[length] = lengths.get(length, 0) + 1
    return lengths


def format_lengths(lengths: dict[int, int]) -> str:
    return ";".join(f"{length}:{lengths[length]}" for length in sorted(lengths))


def parse_guest_audit(path: Path) -> dict[str, list[str]]:
    values: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values.setdefault(key, []).append(value)
    return values


class RedirectRecorder(HTTPRedirectHandler):
    def __init__(self) -> None:
        super().__init__()
        self.redirects: list[str] = []

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        self.redirects.append(f"{code}:{newurl}")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def http_metadata(url: str, method: str = "HEAD") -> dict[str, object]:
    recorder = RedirectRecorder()
    opener = build_opener(recorder)
    request = Request(url, method=method, headers={"User-Agent": "Codex-PC-preflight"})
    try:
        with opener.open(request, timeout=60) as response:
            return {
                "http_status": response.status,
                "final_url": response.geturl(),
                "content_length": response.headers.get("Content-Length", ""),
                "redirect_chain": ";".join(recorder.redirects),
                "reachable": "True",
                "error": "",
            }
    except HTTPError as exc:
        return {
            "http_status": exc.code,
            "final_url": exc.geturl(),
            "content_length": exc.headers.get("Content-Length", ""),
            "redirect_chain": ";".join(recorder.redirects),
            "reachable": "False",
            "error": str(exc),
        }
    except URLError as exc:
        return {
            "http_status": "",
            "final_url": "",
            "content_length": "",
            "redirect_chain": ";".join(recorder.redirects),
            "reachable": "False",
            "error": str(exc.reason),
        }


def find_cellranger_observations() -> tuple[list[Path], list[Path]]:
    executable_hits: list[Path] = []
    r_package_hits: list[Path] = []
    roots = [
        Path(r"I:\codex-config\tools"),
        EXCHANGE_DIR.parent,
        Path(r"H:\\"),
    ]
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("cellranger*"):
            if len(executable_hits) + len(r_package_hits) >= 100:
                break
            path_text = str(path).lower()
            if "library\\cellranger" in path_text or "library/cellranger" in path_text:
                r_package_hits.append(path)
            elif path.is_file() and path.name.lower() in {"cellranger", "cellranger.exe"}:
                executable_hits.append(path)
    return executable_hits, r_package_hits


def find_reference_dirs() -> list[Path]:
    hits: list[Path] = []
    roots = [
        Path(r"I:\codex-config\tools"),
        EXCHANGE_DIR.parent,
        Path(r"H:\\"),
    ]
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("refdata-gex-GRCh38-2020-A"):
            if path.is_dir():
                hits.append(path)
                if len(hits) >= 20:
                    return hits
    return hits


UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

# 1. FASTQ audit.
fastqs = sorted(path for path in FASTQ_ROOT.rglob("*.fastq.gz") if path.is_file())
role_counts = {"R1": 0, "R2": 0, "I1": 0, "I2": 0}
pair_keys = {"R1": set(), "R2": set()}
all_gzip = True
total_fastq_bytes = 0
for fastq in fastqs:
    total_fastq_bytes += fastq.stat().st_size
    with fastq.open("rb") as handle:
        all_gzip = all_gzip and handle.read(2) == b"\x1f\x8b"
    parsed = fastq_role(fastq)
    if parsed:
        prefix, role, chunk = parsed
        role_counts[role] += 1
        if role in pair_keys:
            pair_keys[role].add(f"{fastq.relative_to(FASTQ_ROOT).parent.as_posix()}|{prefix}|{chunk}")
missing_r1 = sorted(pair_keys["R2"] - pair_keys["R1"])
missing_r2 = sorted(pair_keys["R1"] - pair_keys["R2"])
r1_first = next(path for path in fastqs if "_R1_" in path.name)
r2_first = next(path for path in fastqs if "_R2_" in path.name)
previous_copy = read_tsv(PREVIOUS_UPLOAD / "host_fastq_copy_verification.tsv")[0]
fastq_passed = (
    len(fastqs) == 804
    and role_counts["R1"] == 402
    and role_counts["R2"] == 402
    and total_fastq_bytes == 29268262419
    and all_gzip
    and not missing_r1
    and not missing_r2
    and previous_copy["sha256_mismatches"] == "0"
)
write_tsv(
    UPLOAD_DIR / "fastq_reverification_summary.tsv",
    [
        "checked_at",
        "fastq_root",
        "fastq_files",
        "r1_files",
        "r2_files",
        "i1_files",
        "i2_files",
        "total_fastq_bytes",
        "all_files_gzip",
        "all_r1_chunks_have_matching_r2",
        "missing_r1_chunks",
        "missing_r2_chunks",
        "r1_first100_length_distribution",
        "r2_first100_length_distribution",
        "existing_host_sha256_mismatches",
        "fastq_reverification_passed",
    ],
    [
        {
            "checked_at": CHECKED_AT,
            "fastq_root": str(FASTQ_ROOT),
            "fastq_files": len(fastqs),
            "r1_files": role_counts["R1"],
            "r2_files": role_counts["R2"],
            "i1_files": role_counts["I1"],
            "i2_files": role_counts["I2"],
            "total_fastq_bytes": total_fastq_bytes,
            "all_files_gzip": str(all_gzip),
            "all_r1_chunks_have_matching_r2": str(not missing_r1 and not missing_r2),
            "missing_r1_chunks": len(missing_r1),
            "missing_r2_chunks": len(missing_r2),
            "r1_first100_length_distribution": format_lengths(first_100_lengths(r1_first)),
            "r2_first100_length_distribution": format_lengths(first_100_lengths(r2_first)),
            "existing_host_sha256_mismatches": previous_copy["sha256_mismatches"],
            "fastq_reverification_passed": str(fastq_passed),
        }
    ],
)

# 2. Windows host resources.
host_info = powershell_json(
    "$os=Get-CimInstance Win32_OperatingSystem; "
    "$cpu=Get-CimInstance Win32_Processor; "
    "$cs=Get-CimInstance Win32_ComputerSystem; "
    "[pscustomobject]@{"
    "windows_caption=$os.Caption;"
    "windows_version=$os.Version;"
    "windows_build=$os.BuildNumber;"
    "cpu_model=($cpu | Select-Object -First 1 -ExpandProperty Name);"
    "physical_core_count=($cpu | Measure-Object NumberOfCores -Sum).Sum;"
    "logical_processor_count=($cpu | Measure-Object NumberOfLogicalProcessors -Sum).Sum;"
    "total_physical_ram_bytes=[int64]$cs.TotalPhysicalMemory;"
    "available_physical_ram_bytes=[int64]($os.FreePhysicalMemory*1024)"
    "} | ConvertTo-Json -Compress"
)
assert isinstance(host_info, dict)
host_info.update(cpuid_features())
write_tsv(
    UPLOAD_DIR / "host_system_resources.tsv",
    [
        "checked_at",
        "windows_caption",
        "windows_version",
        "windows_build",
        "cpu_model",
        "physical_core_count",
        "logical_processor_count",
        "total_physical_ram_bytes",
        "available_physical_ram_bytes",
        "host_xsave_supported",
        "host_osxsave_supported",
        "host_xcr0_hex",
        "host_avx_supported",
        "host_avx2_supported",
    ],
    [{"checked_at": CHECKED_AT, **host_info}],
)
volume_rows_raw = powershell_json(
    "Get-CimInstance Win32_LogicalDisk -Filter \"DriveType=3\" | "
    "Select-Object DeviceID,VolumeName,FileSystem,@{n='Size';e={[int64]$_.Size}},@{n='FreeSpace';e={[int64]$_.FreeSpace}} | "
    "ConvertTo-Json -Compress"
)
if isinstance(volume_rows_raw, dict):
    volume_rows_raw = [volume_rows_raw]
volume_rows = []
max_free = 0
for row in volume_rows_raw:
    free_space = int(row["FreeSpace"])
    max_free = max(max_free, free_space)
    volume_rows.append(
        {
            "checked_at": CHECKED_AT,
            "device_id": row["DeviceID"],
            "volume_name": row.get("VolumeName", ""),
            "file_system": row.get("FileSystem", ""),
            "total_bytes": int(row["Size"]),
            "free_bytes": free_space,
            "meets_500000000000_free_bytes": str(free_space >= 500000000000),
        }
    )
write_tsv(
    UPLOAD_DIR / "host_fixed_volume_audit.tsv",
    [
        "checked_at",
        "device_id",
        "volume_name",
        "file_system",
        "total_bytes",
        "free_bytes",
        "meets_500000000000_free_bytes",
    ],
    volume_rows,
)

# 3. Linux/QEMU runtime audit.
guest = parse_guest_audit(GUEST_AUDIT)
qemu_code, qemu_version = run_command([str(QEMU_EXE), "--version"], timeout=30)
wsl_status_code, wsl_status = run_command(["wsl.exe", "--status"], timeout=30)
wsl_list_code, wsl_list = run_command(["wsl.exe", "--list", "--verbose"], timeout=30)
docker_version_code, docker_version = run_command(["docker", "--version"], timeout=30)
docker_info_code, docker_info = run_command(["docker", "info"], timeout=30)
docker_service_code, docker_service = run_command(
    [
        "powershell",
        "-NoProfile",
        "-Command",
        "Get-Service -Name Docker,com.docker.service -ErrorAction SilentlyContinue | "
        "Select-Object Name,Status,StartType | Format-Table -AutoSize | Out-String",
    ],
    timeout=30,
)
guest_cpu_count = int(guest.get("NPROC", ["0"])[0])
guest_ram_bytes = int(guest.get("MEMTOTAL_BYTES", ["0"])[0])
guest_avx_visible = guest.get("AVX_VISIBLE", ["False"])[0]
guest_work_free_bytes = ""
df_text = GUEST_AUDIT.read_text(encoding="utf-8", errors="replace")
df_matches = re.findall(r"/dev/vda1\s+\d+\s+\d+\s+(\d+)\s+\d+%\s+/", df_text)
if df_matches:
    guest_work_free_bytes = df_matches[-1]
linux_audit = f"""# Linux runtime resource audit

Generated: {CHECKED_AT}

## WSL

- `wsl.exe --status` exit code: {wsl_status_code}

```text
{wsl_status}
```

- `wsl.exe --list --verbose` exit code: {wsl_list_code}

```text
{wsl_list}
```

## Docker

- `docker --version` exit code: {docker_version_code}

```text
{docker_version}
```

- `docker info` exit code: {docker_info_code}

```text
{docker_info}
```

- Docker service query exit code: {docker_service_code}

```text
{docker_service}
```

## QEMU Ubuntu audit VM

- QEMU executable: `{QEMU_EXE}`
- QEMU version exit code: {qemu_code}

```text
{qemu_version}
```

- Audit run directory: `I:\\codex-config\\tools\\ubuntu-vm\\runs\\20260714_cellranger72_preflight_audit`
- Snapshot mode: true
- FASTQs copied into guest: false
- Guest architecture: `{guest.get('UNAME_M', [''])[0]}`
- Allocated virtual CPU count: `4`
- Guest-visible CPU count: `{guest_cpu_count}`
- Allocated guest RAM bytes: `4294967296`
- Guest MemTotal bytes: `{guest_ram_bytes}`
- Guest MemAvailable bytes: `{guest.get('MEMAVAILABLE_BYTES', [''])[0]}`
- Guest-visible AVX: `{guest_avx_visible}`
- Guest-visible AVX2: `{guest.get('AVX2_VISIBLE', ['False'])[0]}`
- Guest free bytes on audited root/work volume: `{guest_work_free_bytes}`
- `ulimit -n`: `{guest.get('ULIMIT_N', [''])[0]}`
- `ulimit -u`: `{guest.get('ULIMIT_U', [''])[0]}`
- Guest resource audit raw file: `{GUEST_AUDIT}`
- Guest serial log: `{QEMU_SERIAL}`

```text
{GUEST_AUDIT.read_text(encoding='utf-8', errors='replace')}
```
"""
(UPLOAD_DIR / "linux_runtime_resource_audit.txt").write_text(linux_audit, encoding="utf-8", newline="\n")

# 4. Software and source availability.
exe_hits, r_package_hits = find_cellranger_observations()
ref_hits = find_reference_dirs()
cellranger_observed = bool(exe_hits)
cellranger_path = exe_hits[0] if exe_hits else None
cellranger_version = ""
cellranger_help_exit = ""
cellranger_help_observed = "False"
if cellranger_path:
    version_code, version_text = run_command([str(cellranger_path), "--version"], timeout=60)
    help_code, help_text = run_command([str(cellranger_path), "count", "--help"], timeout=60)
    cellranger_version = version_text
    cellranger_help_exit = str(help_code)
    cellranger_help_observed = str(help_code == 0)
    (UPLOAD_DIR / "cellranger_count_help.txt").write_text(help_text, encoding="utf-8", newline="\n")
cellranger_rows = [
    {
        "tool_or_reference": "10x_cellranger_executable",
        "required": "Cell Ranger 7.2.0",
        "observed": str(cellranger_observed),
        "observed_path": "" if cellranger_path is None else str(cellranger_path),
        "observed_version": cellranger_version,
        "bytes": "" if cellranger_path is None else cellranger_path.stat().st_size,
        "sha256": "" if cellranger_path is None else sha256_file(cellranger_path),
        "cellranger_count_help_exit_code": cellranger_help_exit,
        "cellranger_count_help_observed": cellranger_help_observed,
        "notes": "No 10x Cell Ranger executable found under known roots; R package cellranger observations are not the 10x pipeline executable: "
        + ";".join(str(path) for path in r_package_hits[:10]),
    },
    {
        "tool_or_reference": "refdata-gex-GRCh38-2020-A",
        "required": "refdata-gex-GRCh38-2020-A",
        "observed": str(bool(ref_hits)),
        "observed_path": "" if not ref_hits else str(ref_hits[0]),
        "observed_version": "",
        "bytes": "",
        "sha256": "",
        "cellranger_count_help_exit_code": "",
        "cellranger_count_help_observed": "",
        "notes": "No local reference directory observed under known roots" if not ref_hits else "",
    },
]
write_tsv(
    UPLOAD_DIR / "cellranger72_presence_audit.tsv",
    [
        "tool_or_reference",
        "required",
        "observed",
        "observed_path",
        "observed_version",
        "bytes",
        "sha256",
        "cellranger_count_help_exit_code",
        "cellranger_count_help_observed",
        "notes",
    ],
    cellranger_rows,
)

cellranger_page = http_metadata(CELLRANGER_DOWNLOAD_PAGE, "GET")
refdata_head = http_metadata(REFDATA_URL, "HEAD")
source_rows = [
    {
        "source_role": "cellranger_download_page",
        "url": CELLRANGER_DOWNLOAD_PAGE,
        "method": "GET",
        "http_status": cellranger_page["http_status"],
        "final_url": cellranger_page["final_url"],
        "content_length": cellranger_page["content_length"],
        "redirect_chain": cellranger_page["redirect_chain"],
        "reachable": cellranger_page["reachable"],
        "exact_archive_url_obtained_noninteractively": "False",
        "manual_eula_required": "True",
        "notes": "Official download page is accessible, but no exact Cell Ranger 7.2.0 archive URL was obtained without manual EULA flow.",
        "error": cellranger_page["error"],
    },
    {
        "source_role": "refdata-gex-GRCh38-2020-A_archive",
        "url": REFDATA_URL,
        "method": "HEAD",
        "http_status": refdata_head["http_status"],
        "final_url": refdata_head["final_url"],
        "content_length": refdata_head["content_length"],
        "redirect_chain": refdata_head["redirect_chain"],
        "reachable": refdata_head["reachable"],
        "exact_archive_url_obtained_noninteractively": "True",
        "manual_eula_required": "False",
        "notes": "Metadata request only; archive was not downloaded.",
        "error": refdata_head["error"],
    },
]
write_tsv(
    UPLOAD_DIR / "official_source_access_audit.tsv",
    [
        "source_role",
        "url",
        "method",
        "http_status",
        "final_url",
        "content_length",
        "redirect_chain",
        "reachable",
        "exact_archive_url_obtained_noninteractively",
        "manual_eula_required",
        "notes",
        "error",
    ],
    source_rows,
)

# 5. Resource decision and planned command.
storage_gate_passed = max_free >= 500000000000
criteria = [
    {
        "criterion": "fastq_reverification_passed",
        "observed_value": str(fastq_passed),
        "required_value": "True",
        "passed": str(fastq_passed),
        "notes": "",
    },
    {
        "criterion": "guest_visible_cpu_cores",
        "observed_value": guest_cpu_count,
        "required_value": ">=8",
        "passed": str(guest_cpu_count >= 8),
        "notes": "",
    },
    {
        "criterion": "guest_ram_bytes",
        "observed_value": guest_ram_bytes,
        "required_value": ">=64000000000",
        "passed": str(guest_ram_bytes >= 64000000000),
        "notes": "",
    },
    {
        "criterion": "guest_visible_avx",
        "observed_value": guest_avx_visible,
        "required_value": "True",
        "passed": str(guest_avx_visible == "True"),
        "notes": "",
    },
    {
        "criterion": "fixed_volume_free_bytes",
        "observed_value": max_free,
        "required_value": ">=500000000000",
        "passed": str(storage_gate_passed),
        "notes": "maximum fixed-volume free bytes observed",
    },
    {
        "criterion": "cellranger_7_2_0_acquisition_route",
        "observed_value": f"page_reachable={cellranger_page['reachable']};manual_eula_required=True",
        "required_value": "exact Cell Ranger 7.2.0 acquisition route available",
        "passed": "False",
        "notes": "manual EULA flow required; no non-interactive archive URL recorded",
    },
    {
        "criterion": "refdata_gex_grch38_2020_a_source_reachable",
        "observed_value": f"reachable={refdata_head['reachable']};status={refdata_head['http_status']};content_length={refdata_head['content_length']}",
        "required_value": "reachable",
        "passed": str(refdata_head["reachable"] == "True"),
        "notes": "",
    },
]
ready = all(row["passed"] == "True" for row in criteria)
decision = "READY_FOR_CELLRANGER72_ACQUISITION_AND_COUNT_REQUEST" if ready else "LOCAL_CELLRANGER_COUNT_NOT_READY"
criteria.append(
    {
        "criterion": "preflight_decision",
        "observed_value": decision,
        "required_value": "READY_FOR_CELLRANGER72_ACQUISITION_AND_COUNT_REQUEST",
        "passed": str(ready),
        "notes": "No download, installation, count, or quantification was run.",
    }
)
write_tsv(
    UPLOAD_DIR / "preflight_decision.tsv",
    ["criterion", "observed_value", "required_value", "passed", "notes"],
    criteria,
)

planned_command = f"""# Planning record only. Do not execute from this file.
# Cell Ranger 7.2.0 executable: not observed locally.
# refdata-gex-GRCh38-2020-A local path: not observed locally, so the transcriptome placeholder is intentionally retained.
# FASTQ root below is the observed PC-local path containing W6B_0_1_HGGWHDSX2.

cellranger count \\
  --id=R3_SRR23490337 \\
  --transcriptome=<absolute_path_to_refdata-gex-GRCh38-2020-A> \\
  --fastqs={FASTQ_ROOT} \\
  --sample=bamtofastq \\
  --chemistry=auto \\
  --localcores={guest_cpu_count} \\
  --localmem={max(1, guest_ram_bytes // 1000000000)} \\
  --disable-ui
"""
(UPLOAD_DIR / "planned_cellranger_count_command.txt").write_text(
    planned_command, encoding="utf-8", newline="\n"
)

failed = [
    row["criterion"]
    for row in criteria
    if row["passed"] != "True" and row["criterion"] != "preflight_decision"
]
status = f"""# Cell Ranger 7.2.0 preflight status

Status: {decision}
Completed: {CHECKED_AT}

## Scope

- Request: `20260714_prjna932556_cellranger72_preflight_request.md`
- Sample: `R3`
- Run: `SRR23490337`
- Audit only: no downloads, no installation, no Cell Ranger count, no alternative quantification
- Full FASTQ upload to GitHub: no

## FASTQ Reverification

- FASTQ root: `{FASTQ_ROOT}`
- Gzip FASTQ files: `{len(fastqs)}`
- R1 files: `{role_counts['R1']}`
- R2 files: `{role_counts['R2']}`
- Total FASTQ bytes: `{total_fastq_bytes}`
- Pairing passed: `{not missing_r1 and not missing_r2}`
- Existing host SHA256 mismatches: `{previous_copy['sha256_mismatches']}`

## Resource Gate

- Guest-visible CPU cores: `{guest_cpu_count}` / required `8`
- Guest RAM bytes: `{guest_ram_bytes}` / required `64000000000`
- Guest AVX visible: `{guest_avx_visible}`
- Maximum fixed-volume free bytes: `{max_free}` / required `500000000000`
- Cell Ranger 7.2.0 local executable: `{cellranger_observed}`
- Cell Ranger archive route: manual EULA flow required
- GRCh38-2020-A source reachable: `{refdata_head['reachable']}`

## Decision

`{decision}`

Failed criteria: `{", ".join(failed)}`
"""
(UPLOAD_DIR / "STATUS.md").write_text(status, encoding="utf-8", newline="\n")

root_status = f"""# Pan-GI ICI CCC Desktop Exchange

Updated: {datetime.now().date().isoformat()}

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

- Request: `20260714_prjna932556_cellranger72_preflight_request.md`
- Manifest: `20260714_prjna932556_cellranger72_preflight_manifest.tsv`
- PC result: `{decision}`
- Audit upload: `uploads/{REQUEST_NAME}`
- Stop boundary honored: no large downloads, no installation, no Cell Ranger count, no other PRJNA932556 runs

## Most Recent Completed Request

The current request, `20260714_prjna932556_cellranger72_preflight_request.md`, is complete as an audit. The PC-local `SRR23490337` FASTQs were reverified. The current QEMU audit environment does not satisfy the requested Cell Ranger count gate because it exposes {guest_cpu_count} guest CPU cores and {guest_ram_bytes} guest RAM bytes, and the largest fixed-volume free space observed is {max_free} bytes. The decision is `{decision}`.

## Completed Upload Directories

- `uploads/{REQUEST_NAME}` (Cell Ranger 7.2.0 preflight audit)
- `uploads/20260714_prjna932556_srr23490337_full_bamtofastq` (full bamtofastq success audit; full FASTQs PC-local only)
- `uploads/20260713_prjna932556_bamtofastq_feasibility`
- `uploads/20260713_prjna932556_ena_bam_pair_download`
- `uploads/20260711_prjna932556_ena_bam_acquisition`
- `uploads/20260711_prjna932556_windows_sra_pilot`
- `uploads/20260711_gse235863_timepoint_prjna932556`
- `uploads/20260710_gse189926_r_only_rerun`
- `uploads/20260710_gse189926_annotation_qc_refinement`
- `uploads/20260709_metadata_download`
- `uploads/20260709_small_file_download`
- `uploads/20260709_matrix_download`
- `uploads/20260709_preprocess_priority_response`
- `uploads/20260709_response_object_construction`
"""
(EXCHANGE_DIR / "STATUS.md").write_text(root_status, encoding="utf-8-sig", newline="\n")

manifest_rows = []
repo_upload_root = (
    "bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/"
    f"desktop_exchange/uploads/{REQUEST_NAME}"
)
for path in sorted(p for p in UPLOAD_DIR.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt"):
    manifest_rows.append(
        {
            "item_role": "lightweight_upload_file",
            "local_path": str(path),
            "repo_path": f"{repo_upload_root}/{path.name}",
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "uploaded_to_github": "True",
            "notes": "",
        }
    )
for path in [
    EXCHANGE_DIR / f"{REQUEST_NAME}_request.md",
    EXCHANGE_DIR / f"{REQUEST_NAME}_manifest.tsv",
    EXCHANGE_DIR / "STATUS.md",
    EXCHANGE_DIR / "scripts/finalize_20260714_cellranger72_preflight.py",
]:
    manifest_rows.append(
        {
            "item_role": "repo_context_file",
            "local_path": str(path),
            "repo_path": path.relative_to(LOCAL_ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "uploaded_to_github": "True",
            "notes": "",
        }
    )
manifest_rows.append(
    {
        "item_role": "local_only_full_fastq_directory",
        "local_path": str(FASTQ_ROOT),
        "repo_path": "",
        "bytes": total_fastq_bytes,
        "sha256": "",
        "uploaded_to_github": "False",
        "notes": f"{len(fastqs)} gzip FASTQ files; not uploaded",
    }
)
write_tsv(
    UPLOAD_DIR / "upload_manifest.tsv",
    ["item_role", "local_path", "repo_path", "bytes", "sha256", "uploaded_to_github", "notes"],
    manifest_rows,
)
sha_lines = []
for path in sorted(p for p in UPLOAD_DIR.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt"):
    sha_lines.append(f"{sha256_file(path)}  {path.name}")
(UPLOAD_DIR / "SHA256SUMS.txt").write_text("\n".join(sha_lines) + "\n", encoding="utf-8", newline="\n")

print(json.dumps({"decision": decision, "upload_dir": str(UPLOAD_DIR), "failed": failed}, indent=2))
