from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import sys
import time


if len(sys.argv) != 4:
    raise SystemExit(
        "usage: verify_20260714_full_bamtofastq_host_copy.py "
        "<fastq_output_dir> <expected_sha256_tsv> <audit_dir>"
    )

fastq_output_dir = Path(sys.argv[1])
expected_sha256_tsv = Path(sys.argv[2])
audit_dir = Path(sys.argv[3])
detail_path = audit_dir / "host_fastq_sha256_verification.tsv"
summary_path = audit_dir / "host_fastq_copy_verification.tsv"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


expected_rows: list[dict[str, str]] = []
with expected_sha256_tsv.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {"relative_fastq_path", "file_size_bytes", "sha256"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"expected SHA256 table missing columns: {sorted(missing)}")
    for row in reader:
        expected_rows.append(row)

observed_files = sorted(
    path.relative_to(fastq_output_dir).as_posix()
    for path in fastq_output_dir.rglob("*.fastq.gz")
    if path.is_file()
)
expected_paths = [row["relative_fastq_path"] for row in expected_rows]
expected_path_set = set(expected_paths)
observed_path_set = set(observed_files)

missing_paths = sorted(expected_path_set - observed_path_set)
extra_paths = sorted(observed_path_set - expected_path_set)
size_mismatches = 0
sha256_mismatches = 0
verified_files = 0
expected_total_bytes = sum(int(row["file_size_bytes"]) for row in expected_rows)
observed_total_bytes = 0

with detail_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "relative_fastq_path",
            "expected_size_bytes",
            "observed_size_bytes",
            "expected_sha256",
            "observed_sha256",
            "verification_status",
        ]
    )
    last_progress = time.monotonic()
    for row in expected_rows:
        relative_path = row["relative_fastq_path"]
        expected_size = int(row["file_size_bytes"])
        expected_sha256 = row["sha256"]
        host_path = fastq_output_dir / Path(relative_path)
        if not host_path.is_file():
            writer.writerow(
                [relative_path, expected_size, "", expected_sha256, "", "missing"]
            )
            continue
        observed_size = host_path.stat().st_size
        observed_total_bytes += observed_size
        observed_sha256 = sha256_file(host_path)
        status_parts: list[str] = []
        if observed_size != expected_size:
            size_mismatches += 1
            status_parts.append("size_mismatch")
        if observed_sha256 != expected_sha256:
            sha256_mismatches += 1
            status_parts.append("sha256_mismatch")
        if not status_parts:
            verified_files += 1
            status_parts.append("ok")
        writer.writerow(
            [
                relative_path,
                expected_size,
                observed_size,
                expected_sha256,
                observed_sha256,
                ";".join(status_parts),
            ]
        )
        if verified_files % 25 == 0 or time.monotonic() - last_progress >= 30:
            handle.flush()
            print(
                f"verified={verified_files}/{len(expected_rows)} "
                f"bytes={observed_total_bytes}/{expected_total_bytes}",
                file=sys.stderr,
                flush=True,
            )
            last_progress = time.monotonic()
    for relative_path in extra_paths:
        host_path = fastq_output_dir / Path(relative_path)
        observed_size = host_path.stat().st_size
        observed_total_bytes += observed_size
        writer.writerow([relative_path, "", observed_size, "", sha256_file(host_path), "extra"])

all_passed = (
    not missing_paths
    and not extra_paths
    and size_mismatches == 0
    and sha256_mismatches == 0
    and verified_files == len(expected_rows)
)

with summary_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "expected_fastq_files",
            "observed_fastq_files",
            "verified_files",
            "expected_total_bytes",
            "observed_total_bytes",
            "missing_files",
            "extra_files",
            "size_mismatches",
            "sha256_mismatches",
            "host_copy_verification_passed",
        ]
    )
    writer.writerow(
        [
            len(expected_rows),
            len(observed_files),
            verified_files,
            expected_total_bytes,
            observed_total_bytes,
            len(missing_paths),
            len(extra_paths),
            size_mismatches,
            sha256_mismatches,
            str(all_passed),
        ]
    )

if not all_passed:
    raise SystemExit(
        "host copy verification failed: "
        f"missing={len(missing_paths)} extra={len(extra_paths)} "
        f"size_mismatches={size_mismatches} sha256_mismatches={sha256_mismatches}"
    )
