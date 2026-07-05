#!/usr/bin/env python3
from __future__ import annotations

import base64
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "upload_manifest.tsv"
CHUNK_DIR = ROOT / "files" / "base64_chunks"
OUT_DIR = ROOT / "reassembled"


def sha256_bytes(data: bytes) -> str:
    h = hashlib.sha256()
    h.update(data)
    return h.hexdigest()


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["status"] != "packaged":
                continue
            file_id = row["file_id"]
            chunk_count = int(row["chunk_count"])
            pieces = []
            for idx in range(1, chunk_count + 1):
                chunk_path = CHUNK_DIR / f"{file_id}.part{idx:04d}.b64.txt"
                pieces.append(chunk_path.read_text(encoding="ascii").strip())
            data = base64.b64decode("".join(pieces).encode("ascii"))
            observed = sha256_bytes(data)
            expected = row["sha256"]
            if observed != expected:
                raise SystemExit(f"SHA256 mismatch for {file_id}: {observed} != {expected}")
            out_path = OUT_DIR / row["target_path"]
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_bytes(data)
            print(f"reassembled {out_path} ({len(data)} bytes)")


if __name__ == "__main__":
    main()
