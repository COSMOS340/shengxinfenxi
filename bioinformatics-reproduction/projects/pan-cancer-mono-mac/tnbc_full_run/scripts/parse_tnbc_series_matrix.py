#!/usr/bin/env python3
"""Parse GSE169246 series matrix sample metadata for TNBC."""

from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path


def find_project_root() -> Path:
    env_root = os.environ.get("TNBC_PROJECT_ROOT")
    if env_root:
        return Path(env_root)
    current = Path(__file__).resolve().parent
    for _ in range(8):
        package_matrix = current / "inputs" / "GSE169246_series_matrix.txt.gz"
        project_matrix = current / "03_data_raw" / "single_cell" / "GSE169246" / "GSE169246_series_matrix.txt.gz"
        if package_matrix.exists() or project_matrix.exists():
            return current
        if current.parent == current:
            break
        current = current.parent
    raise RuntimeError(f"Could not locate TNBC project/package root from {Path(__file__).resolve().parent}")


PROJECT_ROOT = find_project_root()
if (PROJECT_ROOT / "inputs" / "GSE169246_series_matrix.txt.gz").exists():
    SERIES_MATRIX = PROJECT_ROOT / "inputs" / "GSE169246_series_matrix.txt.gz"
    OUT_PATH = PROJECT_ROOT / "inputs" / "tnbc_gse169246_sample_metadata.tsv"
else:
    SERIES_MATRIX = PROJECT_ROOT / "03_data_raw" / "single_cell" / "GSE169246" / "GSE169246_series_matrix.txt.gz"
    OUT_PATH = PROJECT_ROOT / "00_metadata" / "tnbc_gse169246_sample_metadata.tsv"


def parse_quoted_fields(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t"))


def load_sample_rows() -> list[dict[str, str]]:
    sample_vectors: dict[str, list[str]] = {}
    characteristic_vectors: list[list[str]] = []
    with gzip.open(SERIES_MATRIX, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line.startswith("!Sample_"):
                continue
            fields = parse_quoted_fields(line)
            key = fields[0].lstrip("!")
            values = fields[1:]
            if key == "Sample_characteristics_ch1":
                characteristic_vectors.append(values)
            else:
                sample_vectors[key] = values

    titles = sample_vectors["Sample_title"]
    rows = [{"Sample_title": title} for title in titles]
    for key, values in sample_vectors.items():
        if key == "Sample_title":
            continue
        if len(values) != len(rows):
            raise ValueError(f"{key} length {len(values)} does not match Sample_title length {len(rows)}")
        for row, value in zip(rows, values, strict=True):
            row[key] = value

    for values in characteristic_vectors:
        if len(values) != len(rows):
            raise ValueError(f"Sample_characteristics_ch1 length {len(values)} does not match Sample_title length {len(rows)}")
        parsed_keys = []
        parsed_values = []
        for value in values:
            if ": " not in value:
                parsed_keys.append("characteristics_ch1")
                parsed_values.append(value)
                continue
            key, parsed_value = value.split(": ", 1)
            parsed_keys.append(key)
            parsed_values.append(parsed_value)
        unique_keys = set(parsed_keys)
        if len(unique_keys) != 1:
            raise ValueError(f"Mixed characteristic keys in one vector: {sorted(unique_keys)}")
        out_key = unique_keys.pop()
        for row, parsed_value in zip(rows, parsed_values, strict=True):
            row[out_key] = parsed_value

    for row in rows:
        title = row["Sample_title"]
        sample_title = title.replace(" (ATAC-Seq)", "")
        parts = sample_title.split("_")
        row["sample_id_from_title"] = sample_title
        row["timepoint_from_title"] = parts[0] if len(parts) >= 1 else ""
        row["patient_from_title"] = parts[1] if len(parts) >= 2 else ""
        row["site_suffix_from_title"] = parts[2] if len(parts) >= 3 else ""
        row["is_atac_title"] = "TRUE" if title.endswith(" (ATAC-Seq)") else "FALSE"
    return rows


def main() -> int:
    rows = load_sample_rows()
    fieldnames = [
        "Sample_title",
        "sample_id_from_title",
        "patient_from_title",
        "timepoint_from_title",
        "site_suffix_from_title",
        "is_atac_title",
        "Sample_geo_accession",
        "Sample_source_name_ch1",
        "tissue",
        "disease state",
        "group",
        "Sample_library_strategy",
        "Sample_library_source",
        "Sample_library_selection",
        "Sample_relation",
    ]
    extra = sorted({key for row in rows for key in row} - set(fieldnames))
    fieldnames.extend(extra)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT_PATH}")
    print(f"rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
