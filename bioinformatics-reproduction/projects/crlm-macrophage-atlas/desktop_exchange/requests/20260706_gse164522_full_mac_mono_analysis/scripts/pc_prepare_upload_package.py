#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path


REQUEST_DIR = Path(__file__).resolve().parents[1]
DESKTOP_EXCHANGE_DIR = REQUEST_DIR.parents[1]
UPLOAD_DIR = DESKTOP_EXCHANGE_DIR / "uploads" / REQUEST_DIR.name
SOURCE_OUTPUTS = REQUEST_DIR / "outputs"
DEST_OUTPUTS = UPLOAD_DIR / "outputs"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    if not SOURCE_OUTPUTS.exists():
        raise SystemExit(f"Missing outputs directory: {SOURCE_OUTPUTS}")
    DEST_OUTPUTS.mkdir(parents=True, exist_ok=True)

    rows = []
    for source in sorted(SOURCE_OUTPUTS.rglob("*")):
        if not source.is_file():
            continue
        relative = source.relative_to(SOURCE_OUTPUTS)
        dest = DEST_OUTPUTS / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
        rows.append(
            {
                "scope": "upload",
                "path": str(Path("outputs") / relative).replace("\\", "/"),
                "bytes": str(dest.stat().st_size),
                "sha256": sha256_file(dest),
                "note": "committed_to_github",
            }
        )

    manifest = UPLOAD_DIR / "upload_manifest.tsv"
    with manifest.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["scope", "path", "bytes", "sha256", "note"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(UPLOAD_DIR)


if __name__ == "__main__":
    main()
