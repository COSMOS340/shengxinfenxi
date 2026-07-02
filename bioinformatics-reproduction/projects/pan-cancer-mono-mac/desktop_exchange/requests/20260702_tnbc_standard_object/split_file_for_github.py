#!/usr/bin/env python3

import argparse
import hashlib
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_file(input_path: Path, out_dir: Path, chunk_size_mb: int) -> None:
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    chunk_size = chunk_size_mb * 1024 * 1024
    if chunk_size <= 0:
        raise ValueError("chunk_size_mb must be positive")

    rows = []
    with input_path.open("rb") as source:
        index = 0
        while True:
            chunk = source.read(chunk_size)
            if not chunk:
                break
            part_name = f"{input_path.name}.part{index:03d}"
            part_path = out_dir / part_name
            with part_path.open("wb") as target:
                target.write(chunk)
            rows.append((index, part_name, part_path.stat().st_size, sha256_file(part_path)))
            index += 1

    manifest_path = out_dir / f"{input_path.name}.chunks.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("part_index\tpart_file\tbytes\tsha256\n")
        for row in rows:
            handle.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\n")

    checksum_path = out_dir / f"{input_path.name}.sha256"
    with checksum_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{sha256_file(input_path)}  {input_path.name}\n")

    print(f"Wrote {len(rows)} chunks to {out_dir}")
    print(f"Wrote chunk manifest: {manifest_path}")
    print(f"Wrote original checksum: {checksum_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Split a large binary file into GitHub-sized chunks.")
    parser.add_argument("--input", required=True, help="Input file path.")
    parser.add_argument("--out-dir", required=True, help="Output directory for chunks.")
    parser.add_argument("--chunk-size-mb", type=int, default=90, help="Chunk size in MiB.")
    args = parser.parse_args()
    split_file(Path(args.input), Path(args.out_dir), args.chunk_size_mb)


if __name__ == "__main__":
    main()
