#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Split a large output file into GitHub-sized chunks.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--chunk-size-mb", type=int, default=90)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = args.input
    if not source.exists():
        raise SystemExit(f"Input file not found: {source}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    chunk_bytes = args.chunk_size_mb * 1024 * 1024
    rows = []
    with source.open("rb") as handle:
        part_index = 0
        while True:
            data = handle.read(chunk_bytes)
            if not data:
                break
            part_name = f"{source.name}.part{part_index:03d}"
            part_path = args.out_dir / part_name
            part_path.write_bytes(data)
            rows.append({
                "part_index": part_index,
                "file_name": part_name,
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
            })
            part_index += 1

    manifest_path = args.out_dir / f"{source.name}.chunks.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["part_index", "file_name", "bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    checksum_path = args.out_dir / f"{source.name}.sha256"
    checksum_path.write_text(f"{sha256_file(source)}  {source.name}\n", encoding="utf-8")
    print(f"Wrote {len(rows)} chunks to {args.out_dir}")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Wrote source checksum: {checksum_path}")


if __name__ == "__main__":
    main()

