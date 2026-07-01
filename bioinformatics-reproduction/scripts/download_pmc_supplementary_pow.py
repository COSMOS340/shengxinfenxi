#!/usr/bin/env python3
"""Download public PMC supplementary files that use the CloudPMC proof-of-work gate.

This script is for public PMC article assets where the browser receives a
short proof-of-work page before the real PDF/XLSX payload. It does not bypass
publisher authentication, paywalls, institutional access controls, or robots
rules. Use it only for public PMC file URLs that the user is allowed to access.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urlparse
from urllib.request import Request, build_opener


FIELDNAMES = [
    "url",
    "filename",
    "output_path",
    "status",
    "bytes",
    "sha256",
    "pow_difficulty",
    "pow_nonce",
    "started_at",
    "finished_at",
    "error",
]


@dataclass(frozen=True)
class DownloadTarget:
    url: str
    filename: str


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def filename_from_url(url: str) -> str:
    parsed = urlparse(url)
    name = Path(unquote(parsed.path)).name
    if not name:
        raise ValueError(f"Cannot derive filename from URL: {url}")
    return name


def read_targets(path: Path | None, urls: list[str]) -> list[DownloadTarget]:
    targets: list[DownloadTarget] = []
    if path:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise SystemExit(f"URL TSV is empty: {path}")
            if "url" not in reader.fieldnames:
                raise SystemExit(f"URL TSV must contain a 'url' column: {path}")
            for row in reader:
                url = (row.get("url") or "").strip()
                if not url:
                    continue
                filename = (row.get("filename") or "").strip() or filename_from_url(url)
                targets.append(DownloadTarget(url=url, filename=filename))
    for url in urls:
        targets.append(DownloadTarget(url=url, filename=filename_from_url(url)))
    if not targets:
        raise SystemExit("No URLs provided. Use --url-file or --url.")
    return targets


def fetch(url: str, cookie: str | None = None, timeout: int = 60) -> bytes:
    headers = {"User-Agent": "Mozilla/5.0"}
    if cookie:
        headers["Cookie"] = cookie
    request = Request(url, headers=headers)
    return build_opener().open(request, timeout=timeout).read()


def parse_pow(html: str) -> tuple[str, int, str] | None:
    match = re.search(
        r'POW_CHALLENGE = "([^"]+)".*?POW_DIFFICULTY = "([^"]+)".*?POW_COOKIE_NAME = "([^"]+)"',
        html,
        re.S,
    )
    if not match:
        return None
    challenge, difficulty, cookie_name = match.groups()
    return challenge, int(difficulty), cookie_name


def solve_pow(challenge: str, difficulty: int) -> int:
    prefix = "0" * difficulty
    nonce = 0
    while True:
        digest = hashlib.sha256(f"{challenge}{nonce}".encode()).hexdigest()
        if digest.startswith(prefix):
            return nonce
        nonce += 1


def looks_like_payload(data: bytes, filename: str) -> bool:
    lower = filename.lower()
    if lower.endswith(".pdf"):
        return data.startswith(b"%PDF")
    if lower.endswith(".xlsx"):
        return data.startswith(b"PK")
    if lower.endswith(".csv") or lower.endswith(".tsv") or lower.endswith(".txt"):
        return len(data) > 0 and b"<html" not in data[:512].lower()
    return len(data) > 4096 and b"POW_CHALLENGE" not in data[:4096]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def append_manifest(path: Path, row: dict[str, str | int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def download_target(target: DownloadTarget, output_dir: Path, timeout: int, overwrite: bool) -> dict[str, str | int]:
    output_path = output_dir / target.filename
    started_at = now_iso()
    if output_path.exists() and output_path.stat().st_size > 0 and not overwrite:
        data = output_path.read_bytes()
        return {
            "url": target.url,
            "filename": target.filename,
            "output_path": str(output_path),
            "status": "skipped_existing",
            "bytes": len(data),
            "sha256": sha256_bytes(data),
            "pow_difficulty": "",
            "pow_nonce": "",
            "started_at": started_at,
            "finished_at": now_iso(),
            "error": "",
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    pow_difficulty: str | int = ""
    pow_nonce: str | int = ""
    try:
        first = fetch(target.url, timeout=timeout)
        if looks_like_payload(first, target.filename):
            output_path.write_bytes(first)
            status = "downloaded_without_pow"
            payload = first
        else:
            pow_parts = parse_pow(first.decode("utf-8", errors="ignore"))
            if not pow_parts:
                blocked_path = output_path.with_suffix(output_path.suffix + ".blocked.html")
                blocked_path.write_bytes(first)
                raise RuntimeError(f"proof-of-work page not detected; blocked HTML written to {blocked_path}")

            challenge, difficulty, cookie_name = pow_parts
            pow_difficulty = difficulty
            nonce = solve_pow(challenge, difficulty)
            pow_nonce = nonce
            cookie = f"{cookie_name}={challenge},{nonce}"
            payload = fetch(target.url, cookie=cookie, timeout=timeout)
            if not looks_like_payload(payload, target.filename):
                blocked_path = output_path.with_suffix(output_path.suffix + ".blocked_after_pow.html")
                blocked_path.write_bytes(payload)
                raise RuntimeError(f"payload validation failed after proof-of-work; response written to {blocked_path}")
            output_path.write_bytes(payload)
            status = "downloaded_with_pow"

        return {
            "url": target.url,
            "filename": target.filename,
            "output_path": str(output_path),
            "status": status,
            "bytes": len(payload),
            "sha256": sha256_bytes(payload),
            "pow_difficulty": pow_difficulty,
            "pow_nonce": pow_nonce,
            "started_at": started_at,
            "finished_at": now_iso(),
            "error": "",
        }
    except (HTTPError, URLError, TimeoutError, OSError, RuntimeError) as exc:
        return {
            "url": target.url,
            "filename": target.filename,
            "output_path": str(output_path),
            "status": "failed",
            "bytes": 0,
            "sha256": "",
            "pow_difficulty": pow_difficulty,
            "pow_nonce": pow_nonce,
            "started_at": started_at,
            "finished_at": now_iso(),
            "error": str(exc),
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url-file", type=Path, help="TSV with columns 'url' and optional 'filename'.")
    parser.add_argument("--url", action="append", default=[], help="Single PMC file URL. May be repeated.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for downloaded files.")
    parser.add_argument("--manifest", type=Path, help="Download manifest TSV. Defaults to <output-dir>/download_manifest.tsv.")
    parser.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds.")
    parser.add_argument("--overwrite", action="store_true", help="Re-download existing non-empty output files.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.expanduser().resolve()
    manifest = args.manifest.expanduser().resolve() if args.manifest else output_dir / "download_manifest.tsv"
    targets = read_targets(args.url_file.expanduser().resolve() if args.url_file else None, args.url)

    failures = 0
    for index, target in enumerate(targets, start=1):
        print(f"[{index}/{len(targets)}] {target.filename}", flush=True)
        row = download_target(target, output_dir=output_dir, timeout=args.timeout, overwrite=args.overwrite)
        append_manifest(manifest, row)
        print(f"  {row['status']} {row['bytes']} bytes", flush=True)
        if row["status"] == "failed":
            failures += 1
            print(f"  ERROR: {row['error']}", file=sys.stderr, flush=True)

    print(f"Done. Failed files: {failures}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
