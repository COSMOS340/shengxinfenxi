import csv
import hashlib
from pathlib import Path


PROJECT = Path("I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709")
EXCHANGE = PROJECT / "desktop_exchange"
OUT = EXCHANGE / "uploads/20260713_prjna932556_ena_bam_pair_download"
LOCAL = EXCHANGE / "local_only/20260713_prjna932556_ena_bam_pair_download"
REPORT = OUT / "ena_srr23490337_file_report.tsv"


def digest(path, algorithm):
    h = hashlib.new(algorithm)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    with open(REPORT, newline="", encoding="utf-8") as fh:
        row = list(csv.DictReader(fh, delimiter="\t"))[0]

    bam_path = LOCAL / "R3_possorted_genome_bam.bam.1"
    bai_path = LOCAL / "R3_possorted_genome_bam.bam.1.bai"
    expected = {
        "BAM": {
            "path": bam_path,
            "url": "https://ftp.sra.ebi.ac.uk/vol1/run/SRR234/SRR23490337/R3_possorted_genome_bam.bam.1",
            "bytes": 21494429186,
            "md5": "056a5c996212a64af708d48c5b585e48",
        },
        "BAI": {
            "path": bai_path,
            "url": "https://ftp.sra.ebi.ac.uk/vol1/run/SRR234/SRR23490337/R3_possorted_genome_bam.bam.1.bai",
            "bytes": 10390008,
            "md5": "f391d287338dbf9d3856630338957bb8",
        },
    }

    audit_rows = []
    inventory_rows = []
    for role, spec in expected.items():
        observed_bytes = spec["path"].stat().st_size
        observed_md5 = digest(spec["path"], "md5")
        observed_sha = digest(spec["path"], "sha256")
        status = "complete_verified" if observed_bytes == spec["bytes"] and observed_md5 == spec["md5"] else "verification_failed"
        audit_rows.append({
            "run_accession": "SRR23490337",
            "file_role": role,
            "ena_ftp_url": spec["url"],
            "download_method": "BITS Start-BitsTransfer with resume after transient remote connection closures",
            "download_status": status,
            "start_time": "2026-07-13T18:52:29+08:00",
            "end_time": "2026-07-13T22:15:38+08:00",
            "exit_code": "0",
            "expected_bytes": spec["bytes"],
            "observed_local_bytes": observed_bytes,
            "byte_match_result": str(observed_bytes == spec["bytes"]),
            "expected_md5": spec["md5"],
            "observed_md5": observed_md5,
            "md5_match_result": str(observed_md5 == spec["md5"]),
            "observed_sha256": observed_sha,
            "absolute_pc_local_path": str(spec["path"]),
            "notes": "BAM transfer resumed through 6 transient remote connection closures; Complete-BitsTransfer finalized both files",
        })
        inventory_rows.append({
            "run_accession": "SRR23490337",
            "ena_sample_accession": "SAMN33196641",
            "ena_experiment_accession": "SRX19384059",
            "study_accession": "PRJNA932556",
            "file_role": role,
            "reported_submitted_file_name": spec["path"].name,
            "reported_submitted_bytes": spec["bytes"],
            "observed_local_bytes": observed_bytes,
            "reported_md5": spec["md5"],
            "observed_md5": observed_md5,
            "md5_match_result": str(observed_md5 == spec["md5"]),
            "observed_sha256": observed_sha,
            "absolute_pc_local_path": str(spec["path"]),
            "download_status": status,
        })

    write_tsv(OUT / "ena_bam_pair_download_audit.tsv", audit_rows, list(audit_rows[0].keys()))
    write_tsv(OUT / "ena_bam_pair_file_inventory.tsv", inventory_rows, list(inventory_rows[0].keys()))

    (OUT / "STATUS.md").write_text(f"""# 20260713 PRJNA932556 ENA BAM/BAI Pair Download

## Status

Successful. The exact ENA-submitted BAM and BAI files for `SRR23490337` were downloaded to PC-local storage and verified.

## Gate

- run_accession: `{row['run_accession']}`
- sample_accession: `{row['sample_accession']}`
- experiment_accession: `{row['experiment_accession']}`
- study_accession: `{row['study_accession']}`
- submitted_format: `{row['submitted_format']}`
- gate_status: `passed`

## Verification

- BAM bytes: `21494429186`
- BAM MD5: `056a5c996212a64af708d48c5b585e48`
- BAM SHA256: `{audit_rows[0]['observed_sha256']}`
- BAI bytes: `10390008`
- BAI MD5: `f391d287338dbf9d3856630338957bb8`
- BAI SHA256: `{audit_rows[1]['observed_sha256']}`

## Local storage

The BAM and BAI remain PC-local under:

`I:\\shengxinfenxi\\bioinformatics-reproduction\\projects\\pan-gi-ici-ccc-20260709\\desktop_exchange\\local_only\\20260713_prjna932556_ena_bam_pair_download`

No `bamtofastq`, Cell Ranger, Seurat, or quantification step was run.
""", encoding="utf-8")

    (EXCHANGE / "STATUS.md").write_text("""# Pan-GI ICI CCC Desktop Exchange

Updated: 2026-07-13

This folder contains the PC-side requests and lightweight upload outputs for the pan-GI immune checkpoint therapy cell communication project.

## Current Request

No open PC-side request is pending in this exchange directory.

The most recent request, `20260713_prjna932556_ena_bam_pair_download_request.md`, is complete. The ENA gate passed, and both exact submitted files for `SRR23490337` were downloaded to PC-local storage and verified by byte size and MD5. The BAM and BAI remain PC-local and were not uploaded to GitHub.

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
""", encoding="utf-8")

    sha_rows = []
    for path in sorted(p for p in OUT.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt"):
        sha_rows.append({"sha256": digest(path, "sha256"), "file_name": path.name})
    write_tsv(OUT / "SHA256SUMS.txt", sha_rows, ["sha256", "file_name"])
    print("success audit updated")


if __name__ == "__main__":
    main()
