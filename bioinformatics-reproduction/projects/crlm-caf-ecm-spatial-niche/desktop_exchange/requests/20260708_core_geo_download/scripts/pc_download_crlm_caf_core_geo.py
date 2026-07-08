#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import shutil
import sys
import tarfile
import time
import urllib.request
from pathlib import Path


REQUEST_DIR = Path(__file__).resolve().parents[1]
INPUT_TSV = REQUEST_DIR / "inputs" / "core_geo_urls.tsv"
RAW_DIR = REQUEST_DIR / "raw_core_geo"
OUTPUT_DIR = REQUEST_DIR / "outputs"


def find_repo_root(start: Path) -> Path:
    for path in [start, *start.parents]:
        if (path / ".git").exists():
            return path
    raise RuntimeError(f"Cannot find repo root from {start}")


REPO_ROOT = find_repo_root(REQUEST_DIR)


def sha256sum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".part")
    start = partial.stat().st_size if partial.exists() else 0
    headers = {"Range": f"bytes={start}-"} if start else {}
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response, partial.open("ab") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    partial.replace(output)


def copy_or_download(row: dict[str, str]) -> tuple[Path, str]:
    target = RAW_DIR / row["dataset_id"] / row["file_name"]
    existing_rel = row["pc_existing_relative_path"]
    if existing_rel and existing_rel != "none":
        existing = REPO_ROOT / existing_rel
        if existing.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(existing, target)
            return target, "copied_from_previous_pc_raw"
    if target.exists():
        return target, "already_present"
    for attempt in range(1, 6):
        try:
            download(row["url"], target)
            return target, f"downloaded_attempt_{attempt}"
        except Exception as exc:
            if attempt == 5:
                raise
            print(f"retry {attempt} for {row['file_name']}: {exc}", file=sys.stderr)
            time.sleep(10)
    return target, "unreachable"


def integrity_check(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
        return "gzip_pass"
    if path.suffix == ".tar":
        with tarfile.open(path, "r") as tar:
            tar.getmembers()[:1]
        return "tar_pass"
    return "present"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(INPUT_TSV.open(), delimiter="\t"))
    status_path = OUTPUT_DIR / "core_geo_download_status.tsv"
    integrity_path = OUTPUT_DIR / "core_geo_file_integrity.tsv"
    with status_path.open("w", newline="") as status_handle, integrity_path.open("w", newline="") as integrity_handle:
        status_fields = [
            "dataset_id",
            "file_name",
            "url",
            "local_path",
            "status",
            "message",
            "bytes",
            "sha256",
            "expected_bytes",
            "known_sha256",
            "checked_at",
        ]
        integrity_fields = ["dataset_id", "file_name", "integrity_check", "checked_at"]
        status_writer = csv.DictWriter(status_handle, fieldnames=status_fields, delimiter="\t")
        integrity_writer = csv.DictWriter(integrity_handle, fieldnames=integrity_fields, delimiter="\t")
        status_writer.writeheader()
        integrity_writer.writeheader()
        for row in rows:
            checked_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")
            try:
                path, message = copy_or_download(row)
                size = path.stat().st_size
                digest = sha256sum(path)
                expected_bytes = row["expected_bytes"]
                known_sha = row["known_sha256"]
                status = "ok"
                if expected_bytes and size != int(expected_bytes):
                    status = "size_mismatch"
                if known_sha and digest != known_sha:
                    status = "sha256_mismatch"
                check = integrity_check(path)
                status_writer.writerow(
                    {
                        "dataset_id": row["dataset_id"],
                        "file_name": row["file_name"],
                        "url": row["url"],
                        "local_path": str(path),
                        "status": status,
                        "message": message,
                        "bytes": size,
                        "sha256": digest,
                        "expected_bytes": expected_bytes,
                        "known_sha256": known_sha,
                        "checked_at": checked_at,
                    }
                )
                integrity_writer.writerow(
                    {
                        "dataset_id": row["dataset_id"],
                        "file_name": row["file_name"],
                        "integrity_check": check,
                        "checked_at": checked_at,
                    }
                )
            except Exception as exc:
                status_writer.writerow(
                    {
                        "dataset_id": row["dataset_id"],
                        "file_name": row["file_name"],
                        "url": row["url"],
                        "local_path": str(RAW_DIR / row["dataset_id"] / row["file_name"]),
                        "status": "error",
                        "message": repr(exc),
                        "bytes": "",
                        "sha256": "",
                        "expected_bytes": row["expected_bytes"],
                        "known_sha256": row["known_sha256"],
                        "checked_at": checked_at,
                    }
                )
    print(status_path)
    print(integrity_path)


if __name__ == "__main__":
    main()
