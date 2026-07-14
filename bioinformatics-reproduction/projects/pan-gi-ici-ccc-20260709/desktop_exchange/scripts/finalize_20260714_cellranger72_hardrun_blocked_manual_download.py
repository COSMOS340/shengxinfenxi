from __future__ import annotations

import base64
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from datetime import datetime, timezone
from urllib.error import HTTPError, URLError
from urllib.request import Request, build_opener, HTTPRedirectHandler


LOCAL_ROOT = Path(r"I:\shengxinfenxi")
EXCHANGE_DIR = (
    LOCAL_ROOT
    / "bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange"
)
REQUEST_NAME = "20260714_prjna932556_cellranger72_single_sample_hardrun"
UPLOAD_DIR = EXCHANGE_DIR / "uploads" / REQUEST_NAME
FASTQ_DIR = (
    EXCHANGE_DIR
    / "local_only/20260714_prjna932556_srr23490337_full_bamtofastq/W6B_0_1_HGGWHDSX2"
)
EXPECTED_SHA_TSV = (
    EXCHANGE_DIR
    / "uploads/20260714_prjna932556_srr23490337_full_bamtofastq/full_bamtofastq_local_fastq_sha256.tsv"
)
QEMU_EXE = Path(r"I:\codex-config\tools\qemu\20260501\qemu-system-x86_64.exe")
REQUEST_PATH = EXCHANGE_DIR / f"{REQUEST_NAME}_request.md"
MANIFEST_PATH = EXCHANGE_DIR / f"{REQUEST_NAME}_manifest.tsv"
CELLRANGER_DOWNLOAD_PAGE = "https://www.10xgenomics.com/support/software/cell-ranger/downloads"
REFDATA_URL = "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"
STATUS_VALUE = "BLOCKED_CELLRANGER72_MANUAL_DOWNLOAD"


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
        timeout=120,
    )
    if code != 0:
        raise RuntimeError(f"PowerShell command failed with exit {code}: {text}")
    return json.loads(text)


def gzip_lengths_first100(path: Path) -> dict[int, int]:
    lengths: dict[int, int] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        for _ in range(100):
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


def fastq_parse(path: Path) -> tuple[str, str, str] | None:
    match = re.match(r"^(.+_L\d{3})_([RI][12])_(\d{3})\.fastq\.gz$", path.name)
    if not match:
        return None
    return match.group(1), match.group(2), match.group(3)


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
    request = Request(url, method=method, headers={"User-Agent": "Codex-PC-hardrun"})
    try:
        with opener.open(request, timeout=90) as response:
            body_preview_sha256 = ""
            if method == "GET":
                body_preview_sha256 = hashlib.sha256(response.read(1024 * 1024)).hexdigest()
            return {
                "http_status": response.status,
                "final_url": response.geturl(),
                "content_length": response.headers.get("Content-Length", ""),
                "redirect_chain": ";".join(recorder.redirects),
                "reachable": "True",
                "body_preview_sha256": body_preview_sha256,
                "error": "",
            }
    except HTTPError as exc:
        return {
            "http_status": exc.code,
            "final_url": exc.geturl(),
            "content_length": exc.headers.get("Content-Length", ""),
            "redirect_chain": ";".join(recorder.redirects),
            "reachable": "False",
            "body_preview_sha256": "",
            "error": str(exc),
        }
    except URLError as exc:
        return {
            "http_status": "",
            "final_url": "",
            "content_length": "",
            "redirect_chain": ";".join(recorder.redirects),
            "reachable": "False",
            "body_preview_sha256": "",
            "error": str(exc.reason),
        }


def find_cellranger_local() -> tuple[list[Path], list[Path], list[Path]]:
    archives: list[Path] = []
    executables: list[Path] = []
    r_package_hits: list[Path] = []
    roots = [
        Path(r"I:\codex-config\tools"),
        Path(r"I:\shengxinfenxi"),
        Path(r"F:\\"),
        Path(r"H:\\"),
        Path(r"C:\Users\Administrator\Downloads"),
    ]
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            path_text = str(path).lower()
            name = path.name.lower()
            if "library\\cellranger" in path_text or "library/cellranger" in path_text:
                r_package_hits.append(path)
                continue
            if path.is_file() and "cellranger-7.2.0" in name and (
                name.endswith(".tar.gz") or name.endswith(".tgz") or name.endswith(".zip")
            ):
                archives.append(path)
            if path.is_file() and name in {"cellranger", "cellranger.exe"}:
                executables.append(path)
    return archives, executables, r_package_hits


UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

# FASTQ checksum and pairing reverification.
expected_rows = read_tsv(EXPECTED_SHA_TSV)
expected_by_path = {row["relative_fastq_path"].split("/", 1)[-1]: row for row in expected_rows}
fastqs = sorted(path for path in FASTQ_DIR.rglob("*.fastq.gz") if path.is_file())
role_counts = {"R1": 0, "R2": 0, "I1": 0, "I2": 0}
pair_keys = {"R1": set(), "R2": set()}
all_gzip = True
total_bytes = 0
missing_expected = 0
extra_observed = 0
size_mismatches = 0
sha256_mismatches = 0
verified_files = 0
detail_rows: list[dict[str, object]] = []
expected_name_set = set(expected_by_path)
observed_name_set = {path.name for path in fastqs}

for missing_name in sorted(expected_name_set - observed_name_set):
    row = expected_by_path[missing_name]
    missing_expected += 1
    detail_rows.append(
        {
            "fastq_name": missing_name,
            "observed_path": "",
            "expected_size_bytes": row["file_size_bytes"],
            "observed_size_bytes": "",
            "expected_sha256": row["sha256"],
            "observed_sha256": "",
            "verification_status": "missing",
        }
    )

for index, fastq in enumerate(fastqs, 1):
    with fastq.open("rb") as handle:
        gzip_ok = handle.read(2) == b"\x1f\x8b"
    all_gzip = all_gzip and gzip_ok
    parsed = fastq_parse(fastq)
    if parsed:
        prefix, role, chunk = parsed
        role_counts[role] += 1
        if role in pair_keys:
            pair_keys[role].add(f"{prefix}|{chunk}")
    size = fastq.stat().st_size
    total_bytes += size
    observed_sha = sha256_file(fastq)
    expected = expected_by_path.get(fastq.name)
    status_parts: list[str] = []
    expected_size = ""
    expected_sha = ""
    if expected is None:
        extra_observed += 1
        status_parts.append("extra")
    else:
        expected_size = expected["file_size_bytes"]
        expected_sha = expected["sha256"]
        if int(expected_size) != size:
            size_mismatches += 1
            status_parts.append("size_mismatch")
        if expected_sha != observed_sha:
            sha256_mismatches += 1
            status_parts.append("sha256_mismatch")
        if not status_parts:
            verified_files += 1
            status_parts.append("ok")
    detail_rows.append(
        {
            "fastq_name": fastq.name,
            "observed_path": str(fastq),
            "expected_size_bytes": expected_size,
            "observed_size_bytes": size,
            "expected_sha256": expected_sha,
            "observed_sha256": observed_sha,
            "verification_status": ";".join(status_parts),
        }
    )
    if index % 25 == 0:
        print(f"FASTQ_SHA_PROGRESS {index}/{len(fastqs)} bytes={total_bytes}", file=sys.stderr, flush=True)

write_tsv(
    UPLOAD_DIR / "fastq_sha256_reverification_detail.tsv",
    [
        "fastq_name",
        "observed_path",
        "expected_size_bytes",
        "observed_size_bytes",
        "expected_sha256",
        "observed_sha256",
        "verification_status",
    ],
    detail_rows,
)
missing_r1 = sorted(pair_keys["R2"] - pair_keys["R1"])
missing_r2 = sorted(pair_keys["R1"] - pair_keys["R2"])
r1_first = next(path for path in fastqs if "_R1_" in path.name)
r2_first = next(path for path in fastqs if "_R2_" in path.name)
fastq_passed = (
    len(fastqs) == 804
    and role_counts["R1"] == 402
    and role_counts["R2"] == 402
    and total_bytes == 29268262419
    and all_gzip
    and not missing_r1
    and not missing_r2
    and missing_expected == 0
    and extra_observed == 0
    and size_mismatches == 0
    and sha256_mismatches == 0
)
write_tsv(
    UPLOAD_DIR / "fastq_reverification_summary.tsv",
    [
        "checked_at",
        "fastq_directory",
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
        "verified_sha256_files",
        "missing_expected_files",
        "extra_observed_files",
        "size_mismatches",
        "sha256_mismatches",
        "r1_first100_length_distribution",
        "r2_first100_length_distribution",
        "fastq_reverification_passed",
    ],
    [
        {
            "checked_at": CHECKED_AT,
            "fastq_directory": str(FASTQ_DIR),
            "fastq_files": len(fastqs),
            "r1_files": role_counts["R1"],
            "r2_files": role_counts["R2"],
            "i1_files": role_counts["I1"],
            "i2_files": role_counts["I2"],
            "total_fastq_bytes": total_bytes,
            "all_files_gzip": str(all_gzip),
            "all_r1_chunks_have_matching_r2": str(not missing_r1 and not missing_r2),
            "missing_r1_chunks": len(missing_r1),
            "missing_r2_chunks": len(missing_r2),
            "verified_sha256_files": verified_files,
            "missing_expected_files": missing_expected,
            "extra_observed_files": extra_observed,
            "size_mismatches": size_mismatches,
            "sha256_mismatches": sha256_mismatches,
            "r1_first100_length_distribution": format_lengths(gzip_lengths_first100(r1_first)),
            "r2_first100_length_distribution": format_lengths(gzip_lengths_first100(r2_first)),
            "fastq_reverification_passed": str(fastq_passed),
        }
    ],
)

# Host resources and disk snapshot.
host = powershell_json(
    "$os=Get-CimInstance Win32_OperatingSystem; "
    "$cpu=Get-CimInstance Win32_Processor; "
    "$cs=Get-CimInstance Win32_ComputerSystem; "
    "$vols=Get-CimInstance Win32_LogicalDisk -Filter \"DriveType=3\" | "
    "Select-Object DeviceID,VolumeName,FileSystem,@{n='Size';e={[int64]$_.Size}},@{n='FreeSpace';e={[int64]$_.FreeSpace}}; "
    "[pscustomobject]@{"
    "windows_caption=$os.Caption;"
    "windows_version=$os.Version;"
    "windows_build=$os.BuildNumber;"
    "cpu_model=($cpu | Select-Object -First 1 -ExpandProperty Name);"
    "physical_core_count=($cpu | Measure-Object NumberOfCores -Sum).Sum;"
    "logical_processor_count=($cpu | Measure-Object NumberOfLogicalProcessors -Sum).Sum;"
    "total_physical_ram_bytes=[int64]$cs.TotalPhysicalMemory;"
    "available_physical_ram_bytes=[int64]($os.FreePhysicalMemory*1024);"
    "volumes=$vols"
    "} | ConvertTo-Json -Depth 4 -Compress"
)
assert isinstance(host, dict)
volumes = host["volumes"]
if isinstance(volumes, dict):
    volumes = [volumes]
volume_summary = ";".join(
    f"{row['DeviceID']} total={row['Size']} free={row['FreeSpace']}" for row in volumes
)
f_free = next(int(row["FreeSpace"]) for row in volumes if row["DeviceID"] == "F:")
f_virtual_capacity = max(0, f_free - 20000000000)
write_tsv(
    UPLOAD_DIR / "host_resource_before_run.tsv",
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
        "fixed_volumes",
        "authorized_new_file_volume",
        "authorized_volume_free_bytes",
        "authorized_volume_reserve_bytes",
        "planned_sparse_qemu_data_disk_virtual_capacity_bytes",
    ],
    [
        {
            "checked_at": CHECKED_AT,
            "windows_caption": host["windows_caption"],
            "windows_version": host["windows_version"],
            "windows_build": host["windows_build"],
            "cpu_model": host["cpu_model"],
            "physical_core_count": host["physical_core_count"],
            "logical_processor_count": host["logical_processor_count"],
            "total_physical_ram_bytes": host["total_physical_ram_bytes"],
            "available_physical_ram_bytes": host["available_physical_ram_bytes"],
            "fixed_volumes": volume_summary,
            "authorized_new_file_volume": "F:",
            "authorized_volume_free_bytes": f_free,
            "authorized_volume_reserve_bytes": 20000000000,
            "planned_sparse_qemu_data_disk_virtual_capacity_bytes": f_virtual_capacity,
        }
    ],
)

qemu_accel_code, qemu_accel_text = run_command([str(QEMU_EXE), "-accel", "help"], timeout=30)
qemu_version_code, qemu_version_text = run_command([str(QEMU_EXE), "--version"], timeout=30)
qemu_audit = f"""# QEMU acceleration audit

Generated: {CHECKED_AT}

- QEMU executable: `{QEMU_EXE}`
- `qemu-system-x86_64.exe --version` exit code: {qemu_version_code}

```text
{qemu_version_text}
```

- `qemu-system-x86_64.exe -accel help` exit code: {qemu_accel_code}

```text
{qemu_accel_text}
```

- WHPX listed by QEMU: `{str('whpx' in qemu_accel_text.lower())}`
- WHPX short boot test: `not_attempted`
- Reason: `{STATUS_VALUE}` before VM setup, because exact Cell Ranger 7.2.0 archive acquisition requires manual official 10x download flow.
- QEMU data disk created: `False`
- New files written to F: `False`
"""
(UPLOAD_DIR / "qemu_acceleration_audit.txt").write_text(qemu_audit, encoding="utf-8", newline="\n")

write_tsv(
    UPLOAD_DIR / "guest_resource_allocation.tsv",
    [
        "checked_at",
        "requested_first_vcpus",
        "requested_first_ram_gib",
        "retry_vcpus",
        "retry_ram_gib",
        "vm_started",
        "guest_visible_vcpus",
        "guest_memtotal_bytes",
        "localcores_for_cellranger",
        "localmem_for_cellranger_gib",
        "ulimit_n",
        "ulimit_u",
        "reason_not_started",
    ],
    [
        {
            "checked_at": CHECKED_AT,
            "requested_first_vcpus": 6,
            "requested_first_ram_gib": 24,
            "retry_vcpus": 6,
            "retry_ram_gib": 20,
            "vm_started": "False",
            "guest_visible_vcpus": "",
            "guest_memtotal_bytes": "",
            "localcores_for_cellranger": "",
            "localmem_for_cellranger_gib": "",
            "ulimit_n": "",
            "ulimit_u": "",
            "reason_not_started": STATUS_VALUE,
        }
    ],
)

# Cell Ranger local search and official page metadata.
archives, executables, r_package_hits = find_cellranger_local()
download_page = http_metadata(CELLRANGER_DOWNLOAD_PAGE, "GET")
cellranger_archive_observed = bool(archives)
cellranger_executable_observed = bool(executables)
software_rows = [
    {
        "checked_at": CHECKED_AT,
        "component": "Cell Ranger 7.2.0 archive",
        "required": "cellranger-7.2.0 official archive",
        "observed": str(cellranger_archive_observed),
        "path": "" if not archives else str(archives[0]),
        "bytes": "" if not archives else archives[0].stat().st_size,
        "sha256": "" if not archives else sha256_file(archives[0]),
        "version_stdout": "",
        "count_help_exit_code": "",
        "manual_download_required": str(not cellranger_archive_observed),
        "notes": "No legitimate Cell Ranger 7.2.0 archive found in search roots.",
    },
    {
        "checked_at": CHECKED_AT,
        "component": "Cell Ranger executable",
        "required": "cellranger cellranger-7.2.0",
        "observed": str(cellranger_executable_observed),
        "path": "" if not executables else str(executables[0]),
        "bytes": "" if not executables else executables[0].stat().st_size,
        "sha256": "" if not executables else sha256_file(executables[0]),
        "version_stdout": "",
        "count_help_exit_code": "",
        "manual_download_required": str(not cellranger_archive_observed),
        "notes": "R package cellranger observations are not the 10x pipeline executable: "
        + ";".join(str(path) for path in r_package_hits[:10]),
    },
    {
        "checked_at": CHECKED_AT,
        "component": "Official Cell Ranger download page",
        "required": CELLRANGER_DOWNLOAD_PAGE,
        "observed": download_page["reachable"],
        "path": download_page["final_url"],
        "bytes": download_page["content_length"],
        "sha256": download_page["body_preview_sha256"],
        "version_stdout": "",
        "count_help_exit_code": "",
        "manual_download_required": "True",
        "notes": "Download page reachable, but no exact Cell Ranger 7.2.0 archive URL was obtained without manual official 10x EULA/download flow. No identity, email, institution, or EULA response was invented.",
    },
]
write_tsv(
    UPLOAD_DIR / "software_install_audit.tsv",
    [
        "checked_at",
        "component",
        "required",
        "observed",
        "path",
        "bytes",
        "sha256",
        "version_stdout",
        "count_help_exit_code",
        "manual_download_required",
        "notes",
    ],
    software_rows,
)

ref_head = http_metadata(REFDATA_URL, "HEAD")
write_tsv(
    UPLOAD_DIR / "reference_install_audit.tsv",
    [
        "checked_at",
        "component",
        "source_url",
        "http_status",
        "final_url",
        "content_length",
        "redirect_chain",
        "reachable",
        "archive_downloaded",
        "archive_path",
        "archive_bytes",
        "archive_sha256",
        "extracted_path",
        "extracted_bytes",
        "metadata_files",
        "notes",
    ],
    [
        {
            "checked_at": CHECKED_AT,
            "component": "refdata-gex-GRCh38-2020-A",
            "source_url": REFDATA_URL,
            "http_status": ref_head["http_status"],
            "final_url": ref_head["final_url"],
            "content_length": ref_head["content_length"],
            "redirect_chain": ref_head["redirect_chain"],
            "reachable": ref_head["reachable"],
            "archive_downloaded": "False",
            "archive_path": "",
            "archive_bytes": "",
            "archive_sha256": "",
            "extracted_path": "",
            "extracted_bytes": "",
            "metadata_files": "",
            "notes": f"Reference download not started because status is {STATUS_VALUE}. Metadata HEAD request only.",
        }
    ],
)

component_rows = [
    {
        "checked_at": CHECKED_AT,
        "component": "verified_fastq_directory",
        "path": str(FASTQ_DIR),
        "bytes": total_bytes,
        "sha256": "",
        "stored_on": "I:",
        "uploaded_to_github": "False",
        "notes": "Existing verified input; not copied.",
    },
    {
        "checked_at": CHECKED_AT,
        "component": "cellranger_archive",
        "path": "",
        "bytes": 0,
        "sha256": "",
        "stored_on": "F:",
        "uploaded_to_github": "False",
        "notes": STATUS_VALUE,
    },
    {
        "checked_at": CHECKED_AT,
        "component": "reference_archive",
        "path": "",
        "bytes": 0,
        "sha256": "",
        "stored_on": "F:",
        "uploaded_to_github": "False",
        "notes": "not_downloaded",
    },
    {
        "checked_at": CHECKED_AT,
        "component": "qemu_data_disk",
        "path": "",
        "bytes": 0,
        "sha256": "",
        "stored_on": "F:",
        "uploaded_to_github": "False",
        "notes": "not_created",
    },
    {
        "checked_at": CHECKED_AT,
        "component": "cellranger_pipestance",
        "path": "",
        "bytes": 0,
        "sha256": "",
        "stored_on": "F:",
        "uploaded_to_github": "False",
        "notes": "not_created",
    },
]
write_tsv(
    UPLOAD_DIR / "disk_component_sizes.tsv",
    ["checked_at", "component", "path", "bytes", "sha256", "stored_on", "uploaded_to_github", "notes"],
    component_rows,
)
write_tsv(
    UPLOAD_DIR / "disk_usage_timeseries.tsv",
    [
        "checked_at",
        "phase",
        "host_f_free_bytes",
        "host_f_required_reserve_bytes",
        "guest_work_free_bytes",
        "qemu_host_memory_bytes",
        "notes",
    ],
    [
        {
            "checked_at": CHECKED_AT,
            "phase": "before_setup",
            "host_f_free_bytes": f_free,
            "host_f_required_reserve_bytes": 20000000000,
            "guest_work_free_bytes": "",
            "qemu_host_memory_bytes": "",
            "notes": STATUS_VALUE,
        }
    ],
)

(UPLOAD_DIR / "cellranger_count_command.txt").write_text(
    f"""# Cell Ranger count was not executed.
# Status: {STATUS_VALUE}
# Reason: no legitimate Cell Ranger 7.2.0 archive/executable was available locally, and the official download page requires manual 10x EULA/download flow.

cellranger count \\
  --id=R3_SRR23490337 \\
  --transcriptome=<not_observed_refdata-gex-GRCh38-2020-A> \\
  --fastqs=<guest_path_not_created> \\
  --sample=bamtofastq \\
  --chemistry=auto \\
  --localcores=<not_started> \\
  --localmem=<not_started> \\
  --disable-ui
""",
    encoding="utf-8",
    newline="\n",
)
write_tsv(
    UPLOAD_DIR / "cellranger_count_runtime.tsv",
    ["start_time", "end_time", "exit_code", "status", "stdout_path", "stderr_path", "notes"],
    [
        {
            "start_time": "",
            "end_time": CHECKED_AT,
            "exit_code": "NOT_RUN",
            "status": STATUS_VALUE,
            "stdout_path": "",
            "stderr_path": "",
            "notes": "Stopped before installation/count per request boundary.",
        }
    ],
)
(UPLOAD_DIR / "cellranger_count_exit_code.txt").write_text(
    f"NOT_RUN\n{STATUS_VALUE}\n", encoding="utf-8", newline="\n"
)
(UPLOAD_DIR / "chemistry_detection.txt").write_text(
    f"NOT_RUN\n{STATUS_VALUE}\n", encoding="utf-8", newline="\n"
)
(UPLOAD_DIR / "cellranger_count_log_tail.txt").write_text(
    f"No Cell Ranger logs exist because count was not started. Status: {STATUS_VALUE}\n",
    encoding="utf-8",
    newline="\n",
)

status_text = f"""# SRR23490337 Cell Ranger 7.2.0 hard run status

Status: {STATUS_VALUE}
Completed: {CHECKED_AT}

## Scope

- Request: `20260714_prjna932556_cellranger72_single_sample_hardrun_request.md`
- Run: `SRR23490337`
- Sample: `R3`
- FASTQ directory: `{FASTQ_DIR}`
- New F: writes performed: `False`
- Cell Ranger count executed: `False`
- Other PRJNA932556 runs processed: `False`
- Large outputs uploaded to GitHub: `False`

## Input Reverification

- gzip FASTQ files: `{len(fastqs)}`
- R1 files: `{role_counts['R1']}`
- R2 files: `{role_counts['R2']}`
- total FASTQ bytes: `{total_bytes}`
- verified SHA256 files: `{verified_files}`
- SHA256 mismatches: `{sha256_mismatches}`
- pairing passed: `{not missing_r1 and not missing_r2}`
- FASTQ reverification passed: `{fastq_passed}`

## Block Reason

No legitimate Cell Ranger `7.2.0` archive or 10x Cell Ranger executable was observed in the searched local roots. The official 10x download page is reachable, but the exact archive acquisition requires the official manual EULA/download flow. Per request, I did not invent an archive URL, identity, email address, institution, or EULA response.

## Resource Snapshot

- Authorized new-file volume: `F:`
- F: free bytes before setup: `{f_free}`
- F: reserve bytes: `20000000000`
- Sparse QEMU data disk virtual capacity that would have been used: `{f_virtual_capacity}`
- QEMU acceleration options: see `qemu_acceleration_audit.txt`
"""
(UPLOAD_DIR / "STATUS.md").write_text(status_text, encoding="utf-8", newline="\n")

root_status = f"""# Pan-GI ICI CCC Desktop Exchange

Updated: {datetime.now().date().isoformat()}

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

- Request: `20260714_prjna932556_cellranger72_single_sample_hardrun_request.md`
- Manifest: `20260714_prjna932556_cellranger72_single_sample_hardrun_manifest.tsv`
- PC result: `{STATUS_VALUE}`
- Upload: `uploads/{REQUEST_NAME}`
- Stop boundary honored: no Cell Ranger archive URL was invented, no manual EULA response was fabricated, no installation/count was started, no F: task files were created, and no large outputs were uploaded.

## Most Recent Completed Request

The current request stopped at the authorized manual-download boundary. The `SRR23490337` FASTQs were fully reverified first: {len(fastqs)} gzip files, {role_counts['R1']} R1, {role_counts['R2']} R2, {total_bytes} bytes, {sha256_mismatches} SHA256 mismatches. The official 10x download page is reachable, but Cell Ranger 7.2.0 archive acquisition requires manual official 10x EULA/download flow.

## Completed Upload Directories

- `uploads/{REQUEST_NAME}` ({STATUS_VALUE})
- `uploads/20260714_prjna932556_cellranger72_preflight` (Cell Ranger 7.2.0 preflight audit)
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

manifest_rows: list[dict[str, object]] = []
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
    REQUEST_PATH,
    MANIFEST_PATH,
    EXCHANGE_DIR / "STATUS.md",
    EXCHANGE_DIR / "scripts/finalize_20260714_cellranger72_hardrun_blocked_manual_download.py",
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
        "item_role": "local_only_fastq_directory",
        "local_path": str(FASTQ_DIR),
        "repo_path": "",
        "bytes": total_bytes,
        "sha256": "",
        "uploaded_to_github": "False",
        "notes": f"{len(fastqs)} gzip FASTQs; not uploaded",
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

print(
    json.dumps(
        {
            "status": STATUS_VALUE,
            "fastq_reverification_passed": fastq_passed,
            "verified_sha256_files": verified_files,
            "sha256_mismatches": sha256_mismatches,
            "upload_dir": str(UPLOAD_DIR),
        },
        indent=2,
    )
)
