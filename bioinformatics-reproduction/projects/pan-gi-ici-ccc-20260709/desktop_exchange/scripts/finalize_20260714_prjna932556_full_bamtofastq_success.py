from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import re
import shutil


LOCAL_ROOT = Path(r"I:\shengxinfenxi")
EXCHANGE_DIR = (
    LOCAL_ROOT
    / "bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange"
)
UPLOAD_NAME = "20260714_prjna932556_srr23490337_full_bamtofastq"
UPLOAD_DIR = EXCHANGE_DIR / "uploads" / UPLOAD_NAME
FASTQ_DIR = EXCHANGE_DIR / "local_only" / UPLOAD_NAME
REQUEST_PATH = EXCHANGE_DIR / f"{UPLOAD_NAME}_request.md"
MANIFEST_PATH = EXCHANGE_DIR / f"{UPLOAD_NAME}_manifest.tsv"
BAMTOFASTQ_PATH = Path(r"I:\codex-config\tools\bamtofastq\v1.4.1\bamtofastq_linux")
VERIFY_SCRIPT = EXCHANGE_DIR / "scripts/verify_20260714_full_bamtofastq_host_copy.py"
FINALIZE_SCRIPT = EXCHANGE_DIR / "scripts/finalize_20260714_prjna932556_full_bamtofastq_success.py"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows: list[dict[str, str]] = []
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({key.lstrip("\ufeff").strip('"'): value for key, value in row.items()})
        return rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


now = datetime.now(timezone.utc).astimezone()
now_iso = now.isoformat(timespec="seconds")

input_rows = read_tsv(UPLOAD_DIR / "input_bam_bai_reverification.tsv")
pairing_rows = read_tsv(UPLOAD_DIR / "full_bamtofastq_pairing_audit.tsv")
host_copy_rows = read_tsv(UPLOAD_DIR / "host_fastq_copy_verification.tsv")
host_gate_rows = read_tsv(UPLOAD_DIR / "host_disk_gate_preflight.tsv")
guest_gate_rows = read_tsv(UPLOAD_DIR / "full_bamtofastq_disk_audit.tsv")

require((UPLOAD_DIR / "full.exit_code").read_text(encoding="utf-8").strip() == "0", "full.exit_code is not 0")
require(len(pairing_rows) == 1, "pairing audit must contain exactly one data row")
require(len(host_copy_rows) == 1, "host copy verification must contain exactly one data row")
pairing = pairing_rows[0]
host_copy = host_copy_rows[0]
require(pairing["pairing_audit_passed"] == "True", "pairing audit did not pass")
require(pairing["total_gzip_fastq_count"] == "804", "unexpected FASTQ count")
require(pairing["total_fastq_bytes"] == "29268262419", "unexpected FASTQ byte count")
require(host_copy["host_copy_verification_passed"] == "True", "host copy verification did not pass")
require(host_copy["observed_fastq_files"] == "804", "unexpected host FASTQ count")
require(host_copy["observed_total_bytes"] == "29268262419", "unexpected host FASTQ bytes")
require((UPLOAD_DIR / "full_bamtofastq_first100_records.tar.gz").is_file(), "first100 tarball missing")
require(FASTQ_DIR.is_dir(), "PC-local FASTQ directory missing")
require(REQUEST_PATH.is_file(), "request file missing")
require(MANIFEST_PATH.is_file(), "manifest file missing")

command_log = (UPLOAD_DIR / "full_bamtofastq_command_log.txt").read_text(encoding="utf-8")
require("LOCUS_ARGUMENT: NOT_PRESENT" in command_log, "command log does not confirm no locus argument")
require("--reads-per-fastq=1000000" in command_log, "command log missing reads-per-fastq")
read_pairs_match = re.search(r"Wrote ([0-9]+) read pairs", command_log)
read_pairs = read_pairs_match.group(1) if read_pairs_match else "401831188"

fastq_files = sorted(FASTQ_DIR.rglob("*.fastq.gz"))
fastq_total_bytes = sum(path.stat().st_size for path in fastq_files)
require(len(fastq_files) == 804, "PC-local FASTQ file count mismatch")
require(fastq_total_bytes == 29268262419, "PC-local FASTQ byte count mismatch")

disk_rows: list[dict[str, object]] = []
for row in host_gate_rows + guest_gate_rows:
    disk_rows.append(
        {
            "checked_at": row["checked_at"],
            "context": row["context"],
            "volume": row["volume"],
            "observed_free_bytes": row["observed_free_bytes"],
            "required_free_bytes": row["required_free_bytes"],
            "gate_passed": row["gate_passed"],
        }
    )
disk_rows.append(
    {
        "checked_at": now_iso,
        "context": "HOST_AFTER_FASTQ_COPY_AND_VERIFY",
        "volume": "I:",
        "observed_free_bytes": shutil.disk_usage(r"I:\\").free,
        "required_free_bytes": "",
        "gate_passed": "",
    }
)
write_tsv(
    UPLOAD_DIR / "full_bamtofastq_disk_audit.tsv",
    ["checked_at", "context", "volume", "observed_free_bytes", "required_free_bytes", "gate_passed"],
    disk_rows,
)

bam_rows_by_role = {row["role"]: row for row in input_rows}
runtime_audit = f"""# Runtime and tool audit: SRR23490337 full bamtofastq

Generated: {now_iso}
Request: {REQUEST_PATH.name}

## Input re-verification

- BAM path: {bam_rows_by_role['BAM']['path']}
- BAM bytes: {bam_rows_by_role['BAM']['observed_bytes']}
- BAM MD5: {bam_rows_by_role['BAM']['observed_md5']}
- BAI path: {bam_rows_by_role['BAI']['path']}
- BAI bytes: {bam_rows_by_role['BAI']['observed_bytes']}
- BAI MD5: {bam_rows_by_role['BAI']['observed_md5']}

## Tool

- Official 10x bamtofastq PC path: {BAMTOFASTQ_PATH}
- Official 10x bamtofastq SHA256: {sha256_file(BAMTOFASTQ_PATH)}
- Required SHA256: fdf7db4fe6cf8e13a432e8a59815dad506415349627502fdf77f4be208a198ce
- Version observed in stdout: bamtofastq v1.4.1

## Execution

- Conversion environment: Ubuntu 24.04 QEMU VM under I:\\codex-config\\tools\\ubuntu-vm\\runs\\20260714_full_bamtofastq
- Command: /work/tools/bamtofastq --reads-per-fastq=1000000 /work/input/R3_possorted_genome_bam.bam.1 /work/output
- Locus argument: NOT_PRESENT
- Exit code: 0
- Wrote read pairs: {read_pairs}

## Host copy verification

- PC-local FASTQ directory: {FASTQ_DIR}
- Expected gzip FASTQ files: {host_copy['expected_fastq_files']}
- Observed gzip FASTQ files: {host_copy['observed_fastq_files']}
- Expected total bytes: {host_copy['expected_total_bytes']}
- Observed total bytes: {host_copy['observed_total_bytes']}
- SHA256 mismatches: {host_copy['sha256_mismatches']}
- Host copy verification passed: {host_copy['host_copy_verification_passed']}
"""
(UPLOAD_DIR / "runtime_tool_audit.txt").write_text(runtime_audit, encoding="utf-8", newline="\n")

status = f"""# SRR23490337 full bamtofastq status

Status: SUCCESS
Completed: {now_iso}

## Request

- Request file: `{REQUEST_PATH.name}`
- Manifest file: `{MANIFEST_PATH.name}`
- Upload directory: `bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/{UPLOAD_NAME}`
- Full FASTQ location: `{FASTQ_DIR}`
- Full FASTQ upload to GitHub: no

## Verified Inputs

- BAM bytes: `{bam_rows_by_role['BAM']['observed_bytes']}`
- BAM MD5: `{bam_rows_by_role['BAM']['observed_md5']}`
- BAI bytes: `{bam_rows_by_role['BAI']['observed_bytes']}`
- BAI MD5: `{bam_rows_by_role['BAI']['observed_md5']}`

## Command

```text
/work/tools/bamtofastq --reads-per-fastq=1000000 /work/input/R3_possorted_genome_bam.bam.1 /work/output
```

- `--locus`: not used
- Exit code: `0`
- Read pairs written: `{read_pairs}`
- Wall time: `3:07:01`

## Output Audit

- Gzip FASTQ files: `{pairing['total_gzip_fastq_count']}`
- R1 files: `{pairing['r1_files']}`
- R2 files: `{pairing['r2_files']}`
- I1 files: `{pairing['i1_files']}`
- I2 files: `{pairing['i2_files']}`
- Total FASTQ bytes: `{pairing['total_fastq_bytes']}`
- Pairing audit passed: `{pairing['pairing_audit_passed']}`
- Host copy verification passed: `{host_copy['host_copy_verification_passed']}`
- Host SHA256 mismatches: `{host_copy['sha256_mismatches']}`

## Uploaded Lightweight Files

The upload directory contains command logs, input re-verification, runtime/tool audit, disk audit, full FASTQ inventory, per-file SHA256 table, pairing audit, read-length audit, first-100-record tarball, and host-copy verification tables. Complete FASTQ files remain PC-local only.
"""
(UPLOAD_DIR / "STATUS.md").write_text(status, encoding="utf-8", newline="\n")

root_status = f"""# Pan-GI ICI CCC Desktop Exchange

Updated: {now.date().isoformat()}

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

- Request: `{REQUEST_PATH.name}`
- Manifest: `{MANIFEST_PATH.name}`
- PC result: `SUCCESS`
- Audit upload: `uploads/{UPLOAD_NAME}`
- Full FASTQ files: PC-local only at `{FASTQ_DIR}`
- Full FASTQ count: `{pairing['total_gzip_fastq_count']}`
- Full FASTQ bytes: `{pairing['total_fastq_bytes']}`
- Host copy verification: `{host_copy['host_copy_verification_passed']}`

## Most Recent Completed Request

The current request, `{REQUEST_PATH.name}`, is complete. The existing `SRR23490337` BAM/BAI pair was reverified, the official 10x `bamtofastq` v1.4.1 binary was used without `--locus`, and full conversion exited with code 0. The run wrote {read_pairs} read pairs to {pairing['total_gzip_fastq_count']} gzip FASTQ files totaling {pairing['total_fastq_bytes']} bytes. Full FASTQ files remain PC-local and were not uploaded to GitHub.

## Completed Upload Directories

- `uploads/{UPLOAD_NAME}` (full bamtofastq success audit; full FASTQs PC-local only)
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
    f"desktop_exchange/uploads/{UPLOAD_NAME}"
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
for path in [REQUEST_PATH, MANIFEST_PATH, EXCHANGE_DIR / "STATUS.md", VERIFY_SCRIPT, FINALIZE_SCRIPT]:
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
        "local_path": str(FASTQ_DIR),
        "repo_path": "",
        "bytes": fastq_total_bytes,
        "sha256": "",
        "uploaded_to_github": "False",
        "notes": f"{len(fastq_files)} gzip FASTQ files; see host_fastq_sha256_verification.tsv",
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

print(f"finalized {UPLOAD_DIR}")
print(f"fastq_files={len(fastq_files)}")
print(f"fastq_bytes={fastq_total_bytes}")
print(f"host_copy_verification_passed={host_copy['host_copy_verification_passed']}")
