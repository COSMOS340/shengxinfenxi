#!/usr/bin/env python3
"""Download and verify GSE169246 TNBC RNA inputs from GEO."""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DownloadFile:
    name: str
    url: str
    bytes: int
    sha256: str


FILES = [
    DownloadFile(
        name="GSE169246_TNBC_RNA.counts.mtx.gz",
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169246/suppl/GSE169246_TNBC_RNA.counts.mtx.gz",
        bytes=2006028402,
        sha256="ad6c784ec3f4d965e2ad15e22467c66aa185eceead8808b956e75b33cfb7c76d",
    ),
    DownloadFile(
        name="GSE169246_TNBC_RNA.barcode.tsv.gz",
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169246/suppl/GSE169246_TNBC_RNA.barcode.tsv.gz",
        bytes=2149386,
        sha256="c61d02185e54905a1f7ed658d69fd588db2526a9130951272240338c083a9141",
    ),
    DownloadFile(
        name="GSE169246_TNBC_RNA.feature.tsv.gz",
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169246/suppl/GSE169246_TNBC_RNA.feature.tsv.gz",
        bytes=100694,
        sha256="83f6d996d1a539ec02161884a3ee29920ab3c56207a4242dd93ff62be2702b60",
    ),
    DownloadFile(
        name="GSE169246_series_matrix.txt.gz",
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169246/matrix/GSE169246_series_matrix.txt.gz",
        bytes=6354,
        sha256="80d720ca21d28916d26fc31ec7cfb93fcbe4f9445dc16ea2f0382cf36c8f3c88",
    ),
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_file(path: Path, spec: DownloadFile) -> bool:
    if not path.exists():
        return False
    size = path.stat().st_size
    if size != spec.bytes:
        print(f"[check] {spec.name}: size mismatch {size} != {spec.bytes}")
        return False
    digest = sha256_file(path)
    if digest != spec.sha256:
        print(f"[check] {spec.name}: sha256 mismatch {digest} != {spec.sha256}")
        return False
    print(f"[check] {spec.name}: OK")
    return True


def download_with_resume(spec: DownloadFile, out_path: Path, retries: int) -> None:
    part_path = out_path.with_suffix(out_path.suffix + ".part")
    if verify_file(out_path, spec):
        return

    for attempt in range(1, retries + 1):
        existing = part_path.stat().st_size if part_path.exists() else 0
        mode = "ab" if existing else "wb"
        request = urllib.request.Request(spec.url)
        if existing:
            request.add_header("Range", f"bytes={existing}-")
            print(f"[download] {spec.name}: resuming at {existing:,} bytes")
        else:
            print(f"[download] {spec.name}: starting")

        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                status = response.getcode()
                if existing and status == 200:
                    print(f"[download] {spec.name}: server ignored resume; restarting")
                    existing = 0
                    mode = "wb"
                elif existing and status not in (206, 200):
                    raise RuntimeError(f"unexpected HTTP status {status}")

                downloaded = existing
                last_report = time.time()
                with part_path.open(mode) as handle:
                    while True:
                        chunk = response.read(8 * 1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
                        downloaded += len(chunk)
                        now = time.time()
                        if now - last_report >= 10 or downloaded == spec.bytes:
                            pct = 100 * downloaded / spec.bytes
                            print(f"[download] {spec.name}: {downloaded:,}/{spec.bytes:,} bytes ({pct:.1f}%)")
                            last_report = now
        except (urllib.error.URLError, TimeoutError, RuntimeError) as exc:
            print(f"[download] {spec.name}: attempt {attempt}/{retries} failed: {exc}")
            if attempt == retries:
                raise
            time.sleep(min(60, 5 * attempt))
            continue

        if part_path.stat().st_size != spec.bytes:
            print(f"[download] {spec.name}: incomplete size {part_path.stat().st_size} != {spec.bytes}")
            if attempt == retries:
                raise RuntimeError(f"incomplete download for {spec.name}")
            continue

        part_path.replace(out_path)
        if not verify_file(out_path, spec):
            raise RuntimeError(f"verification failed for {spec.name}")
        return


def write_checksum_file(out_dir: Path) -> None:
    checksum_path = out_dir / "checksums.sha256"
    with checksum_path.open("w", encoding="utf-8") as handle:
        for spec in FILES:
            handle.write(f"{spec.sha256}  {spec.name}\n")
    print(f"[write] {checksum_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Download GSE169246 TNBC RNA inputs from GEO.")
    parser.add_argument("--out", default="inputs", help="Output directory. Default: inputs")
    parser.add_argument("--retries", type=int, default=5, help="Download attempts per file. Default: 5")
    args = parser.parse_args()

    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    for spec in FILES:
        download_with_resume(spec, out_dir / spec.name, args.retries)
    write_checksum_file(out_dir)
    print("[done] all TNBC inputs downloaded and verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
