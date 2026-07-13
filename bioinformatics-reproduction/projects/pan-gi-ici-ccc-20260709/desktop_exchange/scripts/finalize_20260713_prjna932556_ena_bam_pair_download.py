import csv
import hashlib
from datetime import datetime
from pathlib import Path


PROJECT = Path("I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709")
EXCHANGE = PROJECT / "desktop_exchange"
OUT = EXCHANGE / "uploads/20260713_prjna932556_ena_bam_pair_download"
LOCAL = EXCHANGE / "local_only/20260713_prjna932556_ena_bam_pair_download"
LOGDIR = EXCHANGE / "logs/20260713_prjna932556_ena_bam_pair_download"
REPORT = OUT / "ena_srr23490337_file_report.tsv"


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
    LOCAL.mkdir(parents=True, exist_ok=True)
    with open(REPORT, newline="", encoding="utf-8") as fh:
        row = list(csv.DictReader(fh, delimiter="\t"))[0]

    bai_path = LOCAL / "R3_possorted_genome_bam.bam.1.bai"
    bai_bytes = bai_path.stat().st_size if bai_path.exists() else 0
    bai_md5 = md5(bai_path) if bai_path.exists() else ""
    bai_sha = sha256(bai_path) if bai_path.exists() else ""

    # BITS stores the incomplete BAM in BIT*.tmp until completion. Its visible
    # size is sparse/full target size, so the audited partial amount is the
    # BytesTransferred value captured from the suspended job.
    partial_bam_transferred = 13917920681
    expected_bam_bytes = 21494429186
    expected_bai_bytes = 10390008
    expected_bam_md5 = "056a5c996212a64af708d48c5b585e48"
    expected_bai_md5 = "f391d287338dbf9d3856630338957bb8"

    audit_rows = [
        {
            "run_accession": "SRR23490337",
            "file_role": "BAM",
            "ena_ftp_url": "https://ftp.sra.ebi.ac.uk/vol1/run/SRR234/SRR23490337/R3_possorted_genome_bam.bam.1",
            "download_method": "BITS Start-BitsTransfer asynchronous",
            "download_status": "incomplete_suspended_remote_connection_closed_repeatedly",
            "start_time": "2026-07-13T18:52:29+08:00",
            "end_time": "2026-07-13T21:09:59+08:00",
            "exit_code": "",
            "expected_bytes": expected_bam_bytes,
            "observed_local_bytes": partial_bam_transferred,
            "byte_match_result": "False",
            "expected_md5": expected_bam_md5,
            "observed_md5": "",
            "md5_match_result": "not_computed_incomplete",
            "observed_sha256": "",
            "absolute_pc_local_path": str(LOCAL / "R3_possorted_genome_bam.bam.1"),
            "notes": "BITS job 0d4845b2-530e-461c-841a-6a7c9066e1ba suspended after repeated RemoteFile transient errors: connection closed prematurely",
        },
        {
            "run_accession": "SRR23490337",
            "file_role": "BAI",
            "ena_ftp_url": "https://ftp.sra.ebi.ac.uk/vol1/run/SRR234/SRR23490337/R3_possorted_genome_bam.bam.1.bai",
            "download_method": "BITS Start-BitsTransfer synchronous",
            "download_status": "complete_verified",
            "start_time": "",
            "end_time": "",
            "exit_code": "0",
            "expected_bytes": expected_bai_bytes,
            "observed_local_bytes": bai_bytes,
            "byte_match_result": str(bai_bytes == expected_bai_bytes),
            "expected_md5": expected_bai_md5,
            "observed_md5": bai_md5,
            "md5_match_result": str(bai_md5 == expected_bai_md5),
            "observed_sha256": bai_sha,
            "absolute_pc_local_path": str(bai_path),
            "notes": "BAI completed after BAM BITS job was suspended",
        },
    ]
    write_tsv(OUT / "ena_bam_pair_download_audit.tsv", audit_rows, list(audit_rows[0].keys()))

    inventory_rows = [
        {
            "run_accession": "SRR23490337",
            "ena_sample_accession": "SAMN33196641",
            "ena_experiment_accession": "SRX19384059",
            "study_accession": "PRJNA932556",
            "file_role": "BAM",
            "reported_submitted_file_name": "R3_possorted_genome_bam.bam.1",
            "reported_submitted_bytes": expected_bam_bytes,
            "observed_local_bytes": partial_bam_transferred,
            "reported_md5": expected_bam_md5,
            "observed_md5": "",
            "md5_match_result": "not_computed_incomplete",
            "observed_sha256": "",
            "absolute_pc_local_path": str(LOCAL / "R3_possorted_genome_bam.bam.1"),
            "download_status": "incomplete_suspended",
        },
        {
            "run_accession": "SRR23490337",
            "ena_sample_accession": "SAMN33196641",
            "ena_experiment_accession": "SRX19384059",
            "study_accession": "PRJNA932556",
            "file_role": "BAI",
            "reported_submitted_file_name": "R3_possorted_genome_bam.bam.1.bai",
            "reported_submitted_bytes": expected_bai_bytes,
            "observed_local_bytes": bai_bytes,
            "reported_md5": expected_bai_md5,
            "observed_md5": bai_md5,
            "md5_match_result": str(bai_md5 == expected_bai_md5),
            "observed_sha256": bai_sha,
            "absolute_pc_local_path": str(bai_path),
            "download_status": "complete_verified",
        },
    ]
    write_tsv(OUT / "ena_bam_pair_file_inventory.tsv", inventory_rows, list(inventory_rows[0].keys()))

    linux_audit = [
        "Runtime audit was inherited from the immediately preceding ENA BAM acquisition task on the same PC unless rerun logs are present.",
        "",
        "No WSL, Docker, Cell Ranger, bamtofastq, Seurat, or quantification installation was performed in this task.",
    ]
    write_tsv(
        OUT / "linux_runtime_audit.txt",
        [{"text": "\n".join(linux_audit)}],
        ["text"],
    )

    status = f"""# 20260713 PRJNA932556 ENA BAM/BAI Pair Download

## Status

Not successful. The ENA row gate passed, and the BAI file completed and matched MD5, but the BAM download did not complete.

## Gate

- run_accession: `{row['run_accession']}`
- sample_accession: `{row['sample_accession']}`
- experiment_accession: `{row['experiment_accession']}`
- study_accession: `{row['study_accession']}`
- submitted_format: `{row['submitted_format']}`
- gate_status: `passed`

## Download

- BAM expected bytes: `{expected_bam_bytes}`
- BAM observed transferred bytes before suspension: `{partial_bam_transferred}`
- BAM MD5: `not_computed_incomplete`
- BAI expected bytes: `{expected_bai_bytes}`
- BAI observed bytes: `{bai_bytes}`
- BAI MD5 match: `{bai_md5 == expected_bai_md5}`

## Decision

`incomplete_bam_download_not_successful`

The BAM/BAI success definition is not met. The BITS BAM transfer was suspended after repeated remote connection closures.
"""
    (OUT / "STATUS.md").write_text(status, encoding="utf-8")

    root_status = """# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-13

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

No open PC-side request is pending in this exchange directory.

The most recent request, `20260713_prjna932556_ena_bam_pair_download_request.md`, was completed as an incomplete download audit, not as a successful BAM acquisition. The ENA gate passed and the BAI file was downloaded and MD5-verified. The BAM BITS download was suspended after repeated remote connection closures at 13,917,920,681 transferred bytes out of 21,494,429,186 expected bytes. No BAM MD5 was computed and no BAM/BAI success condition was claimed.

## Completed Request Files

- `20260713_prjna932556_ena_bam_pair_download_request.md`
- `20260713_prjna932556_ena_bam_pair_download_manifest.tsv`
- `20260711_prjna932556_ena_bam_acquisition_request.md`
- `20260711_prjna932556_ena_bam_acquisition_manifest.tsv`
- `20260711_prjna932556_windows_sra_pilot_request.md`
- `20260711_prjna932556_windows_sra_pilot_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_manifest.tsv`
- `20260711_gse235863_timepoint_prjna932556_request.md`

## Completed Upload Directories

- `uploads/20260713_prjna932556_ena_bam_pair_download`
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
"""
    (EXCHANGE / "STATUS.md").write_text(root_status, encoding="utf-8")

    sha_rows = []
    for path in sorted(p for p in OUT.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt"):
        sha_rows.append({"sha256": sha256(path), "file_name": path.name})
    write_tsv(OUT / "SHA256SUMS.txt", sha_rows, ["sha256", "file_name"])
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
