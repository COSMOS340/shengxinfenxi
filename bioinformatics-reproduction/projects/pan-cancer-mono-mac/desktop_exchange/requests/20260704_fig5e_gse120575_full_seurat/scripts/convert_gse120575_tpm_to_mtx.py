#!/usr/bin/env python3

import csv
import gzip
import urllib.request
from pathlib import Path


REQUEST_DIR = Path(__file__).resolve().parents[1]
RAW_DIR = REQUEST_DIR / "outputs" / "raw"
SPARSE_DIR = REQUEST_DIR / "outputs" / "sparse"

TPM_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
META_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz"

TPM_PATH = RAW_DIR / "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
META_PATH = RAW_DIR / "GSE120575_patient_ID_single_cells.txt.gz"

GENES_PATH = SPARSE_DIR / "genes.tsv"
CELLS_PATH = SPARSE_DIR / "cells.tsv"
CELL_METADATA_PATH = SPARSE_DIR / "cell_metadata.tsv"
MTX_PATH = SPARSE_DIR / "gse120575_tpm.mtx.gz"


def download_file(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return
    print(f"Downloading {url}")
    urllib.request.urlretrieve(url, path)


def parse_cell_metadata(matrix_cells, matrix_sample_ids):
    rows = []
    with gzip.open(META_PATH, "rt", encoding="latin-1") as handle:
        for line in handle:
            rows.append(line.rstrip("\n").split("\t"))

    header_index = None
    for index, row in enumerate(rows):
        if row and row[0] == "Sample name":
            header_index = index
            break
    if header_index is None:
        raise RuntimeError(f"Could not find metadata header in {META_PATH}")

    header = rows[header_index]
    patient_column = "characteristics: patinet ID (Pre=baseline; Post= on treatment)"
    required_columns = [
        "title",
        patient_column,
        "characteristics: response",
        "characteristics: therapy",
        "source name",
    ]
    missing_columns = [column for column in required_columns if column not in header]
    if missing_columns:
        raise RuntimeError("Missing metadata columns: " + ", ".join(missing_columns))

    records_by_cell = {}
    for row in rows[header_index + 1 :]:
        if not row or not row[0].startswith("Sample "):
            continue
        padded = row + [""] * (len(header) - len(row))
        record = dict(zip(header, padded))
        records_by_cell[record["title"]] = record

    if set(matrix_cells) != set(records_by_cell):
        raise RuntimeError(
            f"Cell mismatch: matrix_only={len(set(matrix_cells) - set(records_by_cell))}, "
            f"metadata_only={len(set(records_by_cell) - set(matrix_cells))}"
        )

    with CELL_METADATA_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "cell_id",
                "sample_id",
                "matrix_sample_id",
                "timepoint",
                "patient_id",
                "response",
                "response_group",
                "therapy",
                "source_name",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for cell_id, matrix_sample_id in zip(matrix_cells, matrix_sample_ids):
            record = records_by_cell[cell_id]
            sample_id = record[patient_column]
            normalized_matrix_sample_id = matrix_sample_id
            for suffix in ("_myeloid_enriched", "_T_enriched"):
                if normalized_matrix_sample_id.endswith(suffix):
                    normalized_matrix_sample_id = normalized_matrix_sample_id[: -len(suffix)]
            if normalized_matrix_sample_id != sample_id:
                raise RuntimeError(
                    f"Sample ID mismatch for {cell_id}: matrix={matrix_sample_id}, metadata={sample_id}"
                )

            if sample_id.startswith("Pre_"):
                timepoint = "Pretreatment"
                patient_id = sample_id[len("Pre_") :]
            elif sample_id.startswith("Post_"):
                timepoint = "Posttreatment"
                patient_id = sample_id[len("Post_") :]
            else:
                raise RuntimeError(f"Unexpected sample ID: {sample_id}")

            response = record["characteristics: response"]
            if response == "Responder":
                response_group = "R"
            elif response == "Non-responder":
                response_group = "NR"
            else:
                raise RuntimeError(f"Unexpected response value: {response}")

            writer.writerow(
                {
                    "cell_id": cell_id,
                    "sample_id": sample_id,
                    "matrix_sample_id": matrix_sample_id,
                    "timepoint": timepoint,
                    "patient_id": patient_id,
                    "response": response,
                    "response_group": response_group,
                    "therapy": record["characteristics: therapy"],
                    "source_name": record["source name"],
                }
            )


def read_header():
    with gzip.open(TPM_PATH, "rt", encoding="latin-1") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        sample_row = handle.readline().rstrip("\n").split("\t")
    cells = header[1:]
    sample_ids = sample_row[1:]
    if len(cells) != len(sample_ids):
        raise RuntimeError(f"Header/sample length mismatch: {len(cells)} versus {len(sample_ids)}")
    return cells, sample_ids


def count_genes_and_nonzero(n_cells):
    n_genes = 0
    nonzero = 0
    genes = []
    with gzip.open(TPM_PATH, "rt", encoding="latin-1") as handle:
        handle.readline()
        handle.readline()
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            while len(parts) - 1 > n_cells and parts[-1] == "":
                parts = parts[:-1]
            if len(parts) - 1 != n_cells:
                raise RuntimeError(f"Row length mismatch for {parts[0]}: {len(parts) - 1} versus {n_cells}")
            n_genes += 1
            genes.append(parts[0])
            for value in parts[1:]:
                if value and float(value) != 0.0:
                    nonzero += 1
    return n_genes, nonzero, genes


def write_sparse_matrix(n_cells, n_genes, nonzero):
    with gzip.open(MTX_PATH, "wt", encoding="utf-8") as out_handle:
        out_handle.write("%%MatrixMarket matrix coordinate real general\n")
        out_handle.write(f"{n_genes} {n_cells} {nonzero}\n")
        with gzip.open(TPM_PATH, "rt", encoding="latin-1") as in_handle:
            in_handle.readline()
            in_handle.readline()
            gene_index = 0
            for line in in_handle:
                parts = line.rstrip("\n").split("\t")
                while len(parts) - 1 > n_cells and parts[-1] == "":
                    parts = parts[:-1]
                gene_index += 1
                for cell_index, value in enumerate(parts[1:], start=1):
                    if value and float(value) != 0.0:
                        out_handle.write(f"{gene_index} {cell_index} {value}\n")


def main():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    SPARSE_DIR.mkdir(parents=True, exist_ok=True)
    download_file(TPM_URL, TPM_PATH)
    download_file(META_URL, META_PATH)

    cells, sample_ids = read_header()
    parse_cell_metadata(cells, sample_ids)
    n_genes, nonzero, genes = count_genes_and_nonzero(len(cells))

    with GENES_PATH.open("w", encoding="utf-8") as handle:
        for gene in genes:
            handle.write(gene + "\n")
    with CELLS_PATH.open("w", encoding="utf-8") as handle:
        for cell in cells:
            handle.write(cell + "\n")

    write_sparse_matrix(len(cells), n_genes, nonzero)
    print(f"genes={n_genes}")
    print(f"cells={len(cells)}")
    print(f"nonzero={nonzero}")
    print(f"matrix={MTX_PATH}")
    print(f"metadata={CELL_METADATA_PATH}")


if __name__ == "__main__":
    main()
