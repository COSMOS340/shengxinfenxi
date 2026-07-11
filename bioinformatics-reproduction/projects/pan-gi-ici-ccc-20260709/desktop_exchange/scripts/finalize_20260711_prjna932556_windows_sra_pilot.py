import csv
import gzip
import hashlib
import os
import tarfile
from pathlib import Path


PROJECT = Path("I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709")
EXCHANGE = PROJECT / "desktop_exchange"
OUT = EXCHANGE / "uploads/20260711_prjna932556_windows_sra_pilot"
LOGDIR = EXCHANGE / "logs/20260711_prjna932556_windows_sra_pilot"
FASTQ_DIR = EXCHANGE / "local_only/20260711_prjna932556_windows_sra_pilot/fastq_dump_100000"
TOOL_DIR = Path("I:/codex-config/tools/sratoolkit/sratoolkit.3.4.1-win64")
ARCHIVE = Path("I:/codex-config/tools/sratoolkit.3.4.1-win64.zip")
DOWNLOAD_URL = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.4.1/sratoolkit.3.4.1-win64.zip"
EXPECTED_MD5 = "adb911959f366f3e3de1f140128d4b23"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_text(path):
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FASTQ_DIR.mkdir(parents=True, exist_ok=True)

    install_rows = [{
        "tool_name": "NCBI SRA Toolkit",
        "download_url": DOWNLOAD_URL,
        "archive_name": ARCHIVE.name,
        "archive_path": str(ARCHIVE),
        "archive_size_bytes": ARCHIVE.stat().st_size if ARCHIVE.exists() else "",
        "archive_md5_observed": md5(ARCHIVE) if ARCHIVE.exists() else "",
        "archive_md5_expected": EXPECTED_MD5,
        "archive_md5_match": str(ARCHIVE.exists() and md5(ARCHIVE) == EXPECTED_MD5),
        "install_path": str(TOOL_DIR),
        "bin_path": str(TOOL_DIR / "bin"),
        "fastq_dump_version": "3.4.1",
        "prefetch_version": "3.4.1",
        "vdb_validate_version": "3.4.1",
        "install_status": "installed_user_writable_i_drive"
    }]
    write_tsv(
        OUT / "sra_toolkit_install_audit.tsv",
        install_rows,
        list(install_rows[0].keys()),
    )

    first_cmd = read_text(LOGDIR / "fastq_dump_command.txt").strip()
    first_stdout = read_text(LOGDIR / "fastq_dump_stdout.txt")
    first_stderr = read_text(LOGDIR / "fastq_dump_stderr.txt")
    second_cmd = read_text(LOGDIR / "fastq_dump_second_command.txt").strip()
    second_stdout = read_text(LOGDIR / "fastq_dump_second_stdout.txt")
    second_stderr = read_text(LOGDIR / "fastq_dump_second_stderr.txt")
    command_log = [
        "PRJNA932556 Windows SRA Toolkit read-structure pilot command log",
        "",
        "Tool versions:",
        "fastq-dump.exe : 3.4.1",
        "prefetch.exe : 3.4.1",
        "vdb-validate.exe : 3.4.1",
        "",
        "Attempt 1: requested command as written by Mac request",
        f"command: {first_cmd}",
        "exit_code: 1",
        "reason: official SRA Toolkit 3.4.1 fastq-dump does not accept --include-technical; usage text was printed and no FASTQ files were generated",
        "stdout_begin",
        first_stdout,
        "stdout_end",
        "stderr_begin",
        first_stderr,
        "stderr_end",
        "",
        "Attempt 2: explicit non-silent compatibility command",
        f"command: {second_cmd}",
        "exit_code: terminated_by_pc_timeout_after_no_output",
        "reason: --include-technical was omitted because fastq-dump 3.4.1 has no such option and does not skip technical reads unless --skip-technical is supplied",
        "monitoring: after more than 30 minutes plus an additional 3 minute check, no FASTQ file, stdout, or stderr output had been produced",
        "stdout_begin",
        second_stdout,
        "stdout_end",
        "stderr_begin",
        second_stderr,
        "stderr_end",
        "",
        "No complete SRR23490337 run was downloaded. No Cell Ranger or quantification step was run.",
    ]
    (OUT / "prjna932556_pilot_command_log.txt").write_text("\n".join(command_log), encoding="utf-8")

    fastq_files = sorted([p for p in FASTQ_DIR.iterdir() if p.is_file() and p.suffix.lower() in {".fastq", ".fq"}])
    if fastq_files:
        inventory_rows = []
        for path in fastq_files:
            inventory_rows.append({
                "file_name": path.name,
                "absolute_path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
                "upload_policy": "pc_local_full_fastq_not_uploaded",
                "status": "generated"
            })
    else:
        inventory_rows = [{
            "file_name": "NO_FASTQ_GENERATED",
            "absolute_path": str(FASTQ_DIR),
            "size_bytes": 0,
            "sha256": "",
            "upload_policy": "no_full_fastq_available",
            "status": "remote_partial_extraction_timeout_no_fastq"
        }]
    write_tsv(
        OUT / "prjna932556_pilot_file_inventory.tsv",
        inventory_rows,
        ["file_name", "absolute_path", "size_bytes", "sha256", "upload_policy", "status"],
    )

    metrics_rows = [{
        "file_name": "NO_FASTQ_GENERATED",
        "read_number_or_role": "not_validated",
        "number_of_records": 0,
        "min_read_length": "",
        "median_read_length": "",
        "max_read_length": "",
        "all_files_same_number_of_spots": "not_evaluable",
        "common_spot_identifier_preserved_across_files": "not_evaluable",
        "metric_status": "not_computed_no_fastq_generated"
    }]
    write_tsv(
        OUT / "prjna932556_pilot_read_metrics.tsv",
        metrics_rows,
        list(metrics_rows[0].keys()),
    )

    header_rows = [{
        "file_name": "NO_FASTQ_GENERATED",
        "header_index": "",
        "header": "",
        "spot_identifier": "",
        "header_status": "not_observed_no_fastq_generated"
    }]
    write_tsv(
        OUT / "prjna932556_pilot_header_audit.tsv",
        header_rows,
        list(header_rows[0].keys()),
    )

    role_rows = [{
        "run": "SRR23490337",
        "sample": "R3",
        "requested_spots": 100000,
        "observed_fastq_file_count": len(fastq_files),
        "observed_read_lengths": "not_observed",
        "barcode_role": "not_validated",
        "umi_role": "not_validated",
        "cdna_role": "not_validated",
        "chemistry": "not_validated",
        "role_assessment_status": "not_validated_no_fastq_generated",
        "decision": "remote_partial_extraction_failed_no_reads_observed"
    }]
    write_tsv(
        OUT / "prjna932556_pilot_role_assessment.tsv",
        role_rows,
        list(role_rows[0].keys()),
    )

    note = OUT / "NO_FASTQ_GENERATED.txt"
    note.write_text(
        "No FASTQ records were generated.\n"
        "Attempt 1 failed because fastq-dump 3.4.1 does not accept --include-technical.\n"
        "Attempt 2 was terminated after prolonged remote partial extraction with no FASTQ, stdout, or stderr output.\n",
        encoding="utf-8",
    )
    with tarfile.open(OUT / "prjna932556_pilot_first100_records.tar.gz", "w:gz") as tar:
        tar.add(note, arcname="NO_FASTQ_GENERATED.txt")
    note.unlink()

    status = """# 20260711 PRJNA932556 Windows SRA Toolkit Pilot

## Status

- SRA Toolkit 3.4.1 for Windows was installed under `I:/codex-config/tools/sratoolkit/sratoolkit.3.4.1-win64`.
- The downloaded archive MD5 matched the official checksum recorded in the request audit.
- The exact requested `fastq-dump` command failed because SRA Toolkit 3.4.1 does not support `--include-technical`.
- A documented compatibility command without `--include-technical` was attempted; it was terminated after prolonged remote partial extraction with no FASTQ, stdout, or stderr output.
- No complete `SRR23490337` run was downloaded, no Cell Ranger step was run, and no count matrix was generated.

## Decision

`remote_partial_extraction_failed_no_reads_observed`

Barcode, UMI, cDNA, and chemistry roles remain `not_validated` because no FASTQ records were observed.
"""
    (OUT / "STATUS.md").write_text(status, encoding="utf-8")

    sha_rows = []
    for path in sorted(p for p in OUT.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt"):
        sha_rows.append({"sha256": sha256(path), "file_name": path.name})
    write_tsv(OUT / "SHA256SUMS.txt", sha_rows, ["sha256", "file_name"])

    root_status = EXCHANGE / "STATUS.md"
    root_status.write_text("""# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-11

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

No open PC-side request is pending in this exchange directory.

The most recent request, `20260711_prjna932556_windows_sra_pilot_request.md`, is complete as a failed bounded pilot audit. SRA Toolkit 3.4.1 was installed and checked. The exact requested `fastq-dump` command failed because `--include-technical` is not supported by this version. A documented compatibility command was attempted and then stopped after prolonged remote partial extraction with no FASTQ output. No complete run, quantification, or count matrix was generated.

## Completed Request Files

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

- `uploads/20260711_prjna932556_windows_sra_pilot`
- `uploads/20260711_gse235863_timepoint_prjna932556`
- `uploads/20260710_gse189926_r_only_rerun`
- `uploads/20260710_gse189926_annotation_qc_refinement`
- `uploads/20260709_metadata_download`
- `uploads/20260709_small_file_download`
- `uploads/20260709_matrix_download`
- `uploads/20260709_preprocess_priority_response`
- `uploads/20260709_response_object_construction`

Large RDS/h5ad objects and any full FASTQ outputs remain PC-local and are recorded in the relevant upload inventories when generated.
""", encoding="utf-8")


if __name__ == "__main__":
    main()
