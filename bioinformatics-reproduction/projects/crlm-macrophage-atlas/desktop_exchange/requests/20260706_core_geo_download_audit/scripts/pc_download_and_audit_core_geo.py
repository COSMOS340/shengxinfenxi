#!/usr/bin/env python3
from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import os
import shutil
import subprocess
import tarfile
from pathlib import Path
from urllib.request import Request, urlopen


REQUEST_DIR = Path(__file__).resolve().parents[1]
INPUTS = REQUEST_DIR / "inputs" / "core_geo_urls.tsv"
RAW_DIR = REQUEST_DIR / "raw_core_geo"
OUT_DIR = REQUEST_DIR / "outputs"
RAW_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_with_curl(url: str, dest: Path) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    temp = dest.with_suffix(dest.suffix + ".part")
    curl = shutil.which("curl")
    if curl:
        command = [
            curl,
            "--location",
            "--fail",
            "--retry",
            "10",
            "--retry-delay",
            "5",
            "--continue-at",
            "-",
            "--output",
            str(temp),
            url,
        ]
        subprocess.run(command, check=True)
    else:
        req = Request(url, headers={"User-Agent": "CRLM-PC-download/1.0"})
        mode = "ab" if temp.exists() else "wb"
        with urlopen(req, timeout=120) as response, temp.open(mode) as out:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
    temp.replace(dest)
    return "downloaded"


def gzip_test(path: Path) -> str:
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1024 * 1024):
                pass
        return "pass"
    except Exception as exc:
        return f"fail:{repr(exc)}"


def expression_dimensions(path: Path) -> tuple[int, int]:
    with gzip.open(path, "rt", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        cells = len(header) - 1
        genes = sum(1 for _ in reader)
    return genes, cells


def mtx_dimensions(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            parts = line.strip().split()
            if len(parts) == 3:
                return int(parts[0]), int(parts[1]), int(parts[2])
    raise RuntimeError(f"No MatrixMarket dimension line found in {path}")


def main() -> None:
    now = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    with INPUTS.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    download_rows = []
    integrity_rows = []
    for row in rows:
        accession = row["accession"]
        file_name = row["file_name"]
        url = row["url"]
        dest = RAW_DIR / accession / file_name
        status = "ok"
        message = "already_present" if dest.exists() and dest.stat().st_size > 0 else ""
        try:
            if not message:
                message = download_with_curl(url, dest)
            bytes_size = dest.stat().st_size
            checksum = sha256_file(dest)
        except Exception as exc:
            status = "error"
            message = repr(exc)
            bytes_size = 0
            checksum = "NA"
        download_rows.append(
            {
                "accession": accession,
                "file_name": file_name,
                "url": url,
                "local_path": str(dest),
                "status": status,
                "message": message,
                "bytes": str(bytes_size),
                "sha256": checksum,
                "checked_at": now,
            }
        )
        if dest.suffix == ".gz" and dest.exists() and dest.stat().st_size > 0:
            integrity = gzip_test(dest)
        elif dest.exists() and dest.stat().st_size > 0:
            integrity = "non_gzip_present"
        else:
            integrity = "missing"
        integrity_rows.append(
            {
                "accession": accession,
                "file_name": file_name,
                "integrity_check": integrity,
                "checked_at": now,
            }
        )

    with (OUT_DIR / "core_geo_download_status.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(download_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(download_rows)

    with (OUT_DIR / "core_geo_file_integrity.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(integrity_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(integrity_rows)

    dim_rows = []
    for row in rows:
        if row["accession"] != "GSE164522" or "expression.csv.gz" not in row["file_name"]:
            continue
        path = RAW_DIR / row["accession"] / row["file_name"]
        if not path.exists() or path.stat().st_size == 0:
            continue
        genes, cells = expression_dimensions(path)
        dim_rows.append(
            {
                "accession": row["accession"],
                "file_name": row["file_name"],
                "genes": str(genes),
                "cells": str(cells),
                "checked_at": now,
            }
        )
    if dim_rows:
        with (OUT_DIR / "gse164522_expression_dimensions.tsv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(dim_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(dim_rows)

    mtx_path = RAW_DIR / "GSE178318" / "GSE178318_matrix.mtx.gz"
    if mtx_path.exists() and mtx_path.stat().st_size > 0:
        genes, cells, nonzero = mtx_dimensions(mtx_path)
        with (OUT_DIR / "gse178318_matrix_dimensions.tsv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["accession", "genes", "cells", "nonzero", "checked_at"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"accession": "GSE178318", "genes": genes, "cells": cells, "nonzero": nonzero, "checked_at": now})

    tar_path = RAW_DIR / "GSE225857" / "GSE225857_RAW.tar"
    if tar_path.exists() and tar_path.stat().st_size > 0:
        with tarfile.open(tar_path) as tf, (OUT_DIR / "gse225857_tar_filelist.tsv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["name", "size", "is_file"], delimiter="\t")
            writer.writeheader()
            for member in tf.getmembers():
                writer.writerow({"name": member.name, "size": member.size, "is_file": member.isfile()})

    status = OUT_DIR / "STATUS.md"
    with status.open("w") as f:
        f.write("# CRLM core GEO download audit\\n\\n")
        f.write(f"Completed: {now}\\n\\n")
        f.write("- Raw matrices were kept on PC and should not be committed.\\n")
        f.write("- Upload only the generated `outputs/` directory.\\n")
        f.write("- GSE178318 GEO data usage terms should be recorded before manuscript submission.\\n")

    print(OUT_DIR)


if __name__ == "__main__":
    main()
