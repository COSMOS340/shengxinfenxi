import csv
import hashlib
from datetime import datetime
from pathlib import Path


PROJECT = Path("I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709")
EXCHANGE = PROJECT / "desktop_exchange"
OUT = EXCHANGE / "uploads/20260711_prjna932556_ena_bam_acquisition"
LOGDIR = EXCHANGE / "logs/20260711_prjna932556_ena_bam_acquisition"
LOCAL_DIR = EXCHANGE / "local_only/20260711_prjna932556_ena_bam_acquisition"
REPORT = OUT / "ena_srr23490337_file_report.tsv"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_text(path):
    if not path.exists():
        return ""
    data = path.read_bytes()
    for enc in ("utf-8", "utf-16-le", "utf-16", "gb18030"):
        try:
            return data.decode(enc)
        except UnicodeDecodeError:
            continue
    return data.decode("utf-8", errors="replace")


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_report():
    with open(REPORT, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"Expected one ENA row, observed {len(rows)}")
    return rows[0]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)
    row = parse_report()

    submitted_ftp_entries = [x for x in row["submitted_ftp"].split(";") if x]
    submitted_md5_entries = [x for x in row["submitted_md5"].split(";") if x]
    submitted_bytes_entries = [x for x in row["submitted_bytes"].split(";") if x]
    submitted_format_entries = [x for x in row["submitted_format"].split(";") if x]

    gate_conditions = {
        "run_accession_exact": row["run_accession"] == "SRR23490337",
        "study_accession_exact": row["study_accession"] == "PRJNA932556",
        "single_submitted_ftp": len(submitted_ftp_entries) == 1,
        "single_submitted_md5": len(submitted_md5_entries) == 1,
        "single_submitted_bytes": len(submitted_bytes_entries) == 1,
        "submitted_format_bam": row["submitted_format"] == "BAM",
    }
    gate_pass = all(gate_conditions.values())

    download_rows = [{
        "run_accession": row["run_accession"],
        "sample_accession": row["sample_accession"],
        "experiment_accession": row["experiment_accession"],
        "study_accession": row["study_accession"],
        "submitted_format": row["submitted_format"],
        "submitted_ftp": row["submitted_ftp"],
        "submitted_md5": row["submitted_md5"],
        "submitted_bytes": row["submitted_bytes"],
        "target_volume": "I:",
        "available_bytes_before_download": "",
        "required_bytes_if_single_bam_plus_20pct": "",
        "download_command": "not_run_phase_2_gate_failed",
        "download_start_time": "",
        "download_end_time": "",
        "exit_code": "",
        "final_byte_size": "",
        "download_status": "not_attempted_phase_2_gate_failed",
        "gate_run_accession_exact": str(gate_conditions["run_accession_exact"]),
        "gate_study_accession_exact": str(gate_conditions["study_accession_exact"]),
        "gate_single_submitted_ftp": str(gate_conditions["single_submitted_ftp"]),
        "gate_single_submitted_md5": str(gate_conditions["single_submitted_md5"]),
        "gate_single_submitted_bytes": str(gate_conditions["single_submitted_bytes"]),
        "gate_submitted_format_bam": str(gate_conditions["submitted_format_bam"]),
        "gate_pass": str(gate_pass),
        "failure_reason": "ENA returned submitted_format=BAM;BAI with two submitted_ftp, submitted_md5, and submitted_bytes entries; request forbids selecting another file when phase 2 gate fails",
    }]
    write_tsv(OUT / "ena_bam_download_audit.tsv", download_rows, list(download_rows[0].keys()))

    inventory_rows = [{
        "run_accession": row["run_accession"],
        "ena_sample_accession": row["sample_accession"],
        "ena_experiment_accession": row["experiment_accession"],
        "study_accession": row["study_accession"],
        "reported_submitted_file_name": row["submitted_ftp"],
        "reported_submitted_bytes": row["submitted_bytes"],
        "observed_local_bytes": "",
        "reported_md5": row["submitted_md5"],
        "observed_md5": "",
        "md5_match_result": "not_evaluable_download_not_attempted",
        "observed_sha256": "",
        "absolute_pc_local_path": str(LOCAL_DIR),
        "download_status": "not_attempted_phase_2_gate_failed",
    }]
    write_tsv(OUT / "ena_bam_file_inventory.tsv", inventory_rows, list(inventory_rows[0].keys()))

    linux_lines = []
    for command, filename in [
        ("wsl.exe --status", "wsl_status.txt"),
        ("wsl.exe -l -v", "wsl_list.txt"),
        ("docker.exe --version", "docker_version.txt"),
    ]:
        text = read_text(LOGDIR / filename)
        linux_lines.extend([
            f"command: {command}",
            f"log_file: {filename}",
            "output_begin",
            text,
            "output_end",
            "",
        ])
    (OUT / "linux_runtime_audit.txt").write_text("\n".join(linux_lines), encoding="utf-8")

    status = f"""# 20260711 PRJNA932556 ENA BAM Acquisition

## Status

Stopped before download because the ENA Phase 2 gate failed.

## ENA row

- run_accession: `{row['run_accession']}`
- study_accession: `{row['study_accession']}`
- sample_accession: `{row['sample_accession']}`
- experiment_accession: `{row['experiment_accession']}`
- submitted_format: `{row['submitted_format']}`
- submitted_ftp: `{row['submitted_ftp']}`
- submitted_md5: `{row['submitted_md5']}`
- submitted_bytes: `{row['submitted_bytes']}`

## Gate result

The request required exactly one submitted file, one checksum, one byte-size entry, and submitted format `BAM`. ENA returned `BAM;BAI` with two submitted file entries, two checksums, and two byte-size entries. No file was selected manually and no BAM was downloaded.

## Decision

`not_attempted_phase_2_gate_failed`
"""
    (OUT / "STATUS.md").write_text(status, encoding="utf-8")

    root_status = f"""# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-13

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

No open PC-side request is pending in this exchange directory.

The most recent request, `20260711_prjna932556_ena_bam_acquisition_request.md`, was completed as a gate-failed audit. The ENA read-run report for `SRR23490337` was retrieved exactly, but Phase 2 failed because ENA returned `submitted_format=BAM;BAI` with two submitted file, checksum, and byte-size entries. Per the request, no submitted file was manually selected and no BAM was downloaded.

## Completed Request Files

- `20260711_prjna932556_ena_bam_acquisition_request.md`
- `20260711_prjna932556_ena_bam_acquisition_manifest.tsv`
- `20260711_prjna932556_windows_sra_pilot_request.md`
- `20260711_prjna932556_windows_sra_pilot_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_request.md`
- `20260710_gse189926_r_only_rerun_manifest.tsv`
- `20260710_gse189926_r_only_rerun_request.md`
- `20260710_gse189926_annotation_qc_refinement_manifest.tsv`
- `20260710_gse189926_annotation_qc_refinement_request.md`
- `20260710_gse189926_umap_repair_policy.md`

## Completed Upload Directories

- `uploads/20260711_prjna932556_ena_bam_acquisition`
- `uploads/20260711_prjna932556_windows_sra_pilot`
- `uploads/20260711_gse235863_timepoint_prjna932556`
- `uploads/20260710_gse189926_r_only_rerun`
- `uploads/20260710_gse189926_annotation_qc_refinement`
- `uploads/20260709_metadata_download`
- `uploads/20260709_small_file_download`
- `uploads/20260709_matrix_download`
- `uploads/20260709_preprocess_priority_response`
- `uploads/20260709_response_object_construction`

Large RDS/h5ad objects and any full FASTQ/BAM outputs remain PC-local and are recorded in the relevant upload inventories when generated.
"""
    (EXCHANGE / "STATUS.md").write_text(root_status, encoding="utf-8")

    sha_rows = []
    for path in sorted(p for p in OUT.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt"):
        sha_rows.append({"sha256": sha256(path), "file_name": path.name})
    write_tsv(OUT / "SHA256SUMS.txt", sha_rows, ["sha256", "file_name"])

    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
