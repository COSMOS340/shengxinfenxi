#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import tarfile
import time
from pathlib import Path


REQUEST_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = next(path for path in [REQUEST_DIR, *REQUEST_DIR.parents] if (path / ".git").exists())
PREVIOUS_RAW_DIR = (
    REPO_ROOT
    / "bioinformatics-reproduction"
    / "projects"
    / "crlm-caf-ecm-spatial-niche"
    / "desktop_exchange"
    / "requests"
    / "20260708_core_geo_download"
    / "raw_core_geo"
)
UPLOAD_DIR = (
    REPO_ROOT
    / "bioinformatics-reproduction"
    / "projects"
    / "crlm-caf-ecm-spatial-niche"
    / "desktop_exchange"
    / "uploads"
    / "20260708_raw_structure_audit"
)


def sha256sum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(paths: list[Path]) -> None:
    out = UPLOAD_DIR / "raw_file_manifest.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["dataset_id", "file_name", "path", "size_bytes", "sha256", "exists"],
            delimiter="\t",
        )
        writer.writeheader()
        for path in paths:
            writer.writerow(
                {
                    "dataset_id": path.parent.name,
                    "file_name": path.name,
                    "path": str(path),
                    "size_bytes": path.stat().st_size if path.exists() else "",
                    "sha256": sha256sum(path) if path.exists() else "",
                    "exists": str(path.exists()),
                }
            )


def write_tar_members(tar_path: Path, out_path: Path) -> list[tarfile.TarInfo]:
    members: list[tarfile.TarInfo] = []
    with tarfile.open(tar_path, "r") as tar, out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "size", "is_file"], delimiter="\t")
        writer.writeheader()
        for member in tar:
            members.append(member)
            writer.writerow({"name": member.name, "size": member.size, "is_file": member.isfile()})
    return members


def open_gzip_member(tar: tarfile.TarFile, member: tarfile.TarInfo):
    raw = tar.extractfile(member)
    if raw is None:
        raise RuntimeError(f"Cannot extract member {member.name}")
    return gzip.open(raw, "rt", encoding="utf-8", errors="replace")


def parse_mtx_dimensions_from_text(handle) -> tuple[int, int, int]:
    for line in handle:
        if line.startswith("%"):
            continue
        parts = line.strip().split()
        if len(parts) == 3:
            return int(parts[0]), int(parts[1]), int(parts[2])
        raise ValueError(f"Unexpected matrix dimension line: {line[:80]}")
    raise ValueError("Matrix dimension line not found")


def parse_mtx_dimensions_from_gzip(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return parse_mtx_dimensions_from_text(handle)


def count_text_table(handle) -> tuple[int, int, str]:
    first = handle.readline()
    if not first:
        return 0, 0, ""
    columns = first.rstrip("\n").split("\t")
    rows = 0
    for _ in handle:
        rows += 1
    return rows, len(columns), first.rstrip("\n")[:300]


def count_lines(handle) -> int:
    return sum(1 for _ in handle)


def summarize_tar_inner_files(tar_path: Path, out_path: Path) -> None:
    rows: list[dict[str, str | int]] = []
    with tarfile.open(tar_path, "r") as tar:
        for member in tar:
            if not member.isfile():
                continue
            name = member.name
            row: dict[str, str | int] = {
                "member_name": name,
                "member_size": member.size,
                "summary_type": "not_parsed",
                "primary_count": "",
                "secondary_count": "",
                "third_count": "",
                "notes": "",
            }
            try:
                if name.endswith(".mtx.gz"):
                    with open_gzip_member(tar, member) as handle:
                        n1, n2, n3 = parse_mtx_dimensions_from_text(handle)
                    row.update(
                        {
                            "summary_type": "matrix_market_dimensions",
                            "primary_count": n1,
                            "secondary_count": n2,
                            "third_count": n3,
                        }
                    )
                elif name.endswith(".barcodes.tsv.gz") or name.endswith(".features.tsv.gz"):
                    with open_gzip_member(tar, member) as handle:
                        row.update({"summary_type": "gzip_line_count", "primary_count": count_lines(handle)})
                elif name.endswith("_meta.txt.gz") or name.endswith("_counts.txt.gz"):
                    with open_gzip_member(tar, member) as handle:
                        rows_count, columns_count, header_preview = count_text_table(handle)
                    row.update(
                        {
                            "summary_type": "gzip_tabular_rows_columns",
                            "primary_count": rows_count,
                            "secondary_count": columns_count,
                            "notes": header_preview,
                        }
                    )
            except Exception as exc:
                row.update({"summary_type": "error", "notes": repr(exc)})
            rows.append(row)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "member_name",
                "member_size",
                "summary_type",
                "primary_count",
                "secondary_count",
                "third_count",
                "notes",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_gse178318_dimensions(matrix_path: Path, barcodes_path: Path, genes_path: Path) -> None:
    genes, cells, nonzero = parse_mtx_dimensions_from_gzip(matrix_path)
    with gzip.open(barcodes_path, "rt", encoding="utf-8", errors="replace") as handle:
        barcode_lines = count_lines(handle)
    with gzip.open(genes_path, "rt", encoding="utf-8", errors="replace") as handle:
        gene_lines = count_lines(handle)
    out = UPLOAD_DIR / "gse178318_matrix_dimensions.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["accession", "matrix_genes", "matrix_cells", "matrix_nonzero", "gene_lines", "barcode_lines"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "accession": "GSE178318",
                "matrix_genes": genes,
                "matrix_cells": cells,
                "matrix_nonzero": nonzero,
                "gene_lines": gene_lines,
                "barcode_lines": barcode_lines,
            }
        )


def write_status() -> None:
    out = UPLOAD_DIR / "STATUS.md"
    out.write_text(
        "# CRLM CAF/ECM raw structure audit\n\n"
        f"Completed: {time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n\n"
        f"Raw input directory: `{PREVIOUS_RAW_DIR}`\n\n"
        "Outputs are file manifests, tar member lists, and matrix/table dimensions only. "
        "No raw GEO files were committed.\n",
        encoding="utf-8",
    )


def main() -> None:
    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    gse225857_tar = PREVIOUS_RAW_DIR / "GSE225857" / "GSE225857_RAW.tar"
    gse245552_tar = PREVIOUS_RAW_DIR / "GSE245552" / "GSE245552_RAW.tar"
    gse178318_matrix = PREVIOUS_RAW_DIR / "GSE178318" / "GSE178318_matrix.mtx.gz"
    gse178318_barcodes = PREVIOUS_RAW_DIR / "GSE178318" / "GSE178318_barcodes.tsv.gz"
    gse178318_genes = PREVIOUS_RAW_DIR / "GSE178318" / "GSE178318_genes.tsv.gz"

    write_manifest([gse225857_tar, gse178318_matrix, gse178318_barcodes, gse178318_genes, gse245552_tar])
    write_tar_members(gse225857_tar, UPLOAD_DIR / "gse225857_tar_members.tsv")
    write_tar_members(gse245552_tar, UPLOAD_DIR / "gse245552_tar_members.tsv")
    summarize_tar_inner_files(gse225857_tar, UPLOAD_DIR / "gse225857_inner_file_summaries.tsv")
    summarize_tar_inner_files(gse245552_tar, UPLOAD_DIR / "gse245552_inner_file_summaries.tsv")
    write_gse178318_dimensions(gse178318_matrix, gse178318_barcodes, gse178318_genes)
    write_status()
    print(UPLOAD_DIR)


if __name__ == "__main__":
    main()

