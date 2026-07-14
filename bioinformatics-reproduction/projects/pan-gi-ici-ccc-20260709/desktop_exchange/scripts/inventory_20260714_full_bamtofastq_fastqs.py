from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path
import re
import shutil
import sys
import tarfile


output_dir = Path(sys.argv[1])
audit_dir = Path(sys.argv[2])
inventory_path = audit_dir / "full_bamtofastq_fastq_inventory.tsv"
sha256_path = audit_dir / "full_bamtofastq_local_fastq_sha256.tsv"
pairing_path = audit_dir / "full_bamtofastq_pairing_audit.tsv"
readlength_path = audit_dir / "full_bamtofastq_readlength_audit.tsv"
sample_root = audit_dir / "first100_records"
archive_path = audit_dir / "full_bamtofastq_first100_records.tar.gz"
sample_root.mkdir(parents=True, exist_ok=True)

fastq_pattern = re.compile(
    r"^(?P<prefix>.+_L(?P<lane>\d{3}))_(?P<role>[RI][12])_(?P<chunk>\d{3})\.fastq\.gz$"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


fastqs = sorted(output_dir.rglob("*.fastq.gz"))
if not fastqs:
    shutil.rmtree(sample_root)
    raise SystemExit("No gzip FASTQ files were produced")

role_counts = {"R1": 0, "R2": 0, "I1": 0, "I2": 0}
pair_keys: dict[str, set[str]] = {"R1": set(), "R2": set()}
total_bytes = 0
all_gzip = True
parsed_files = 0

with (
    inventory_path.open("w", encoding="utf-8", newline="") as inventory_handle,
    sha256_path.open("w", encoding="utf-8", newline="") as sha256_handle,
    readlength_path.open("w", encoding="utf-8", newline="") as readlength_handle,
):
    inventory_writer = csv.writer(inventory_handle, delimiter="\t", lineterminator="\n")
    sha256_writer = csv.writer(sha256_handle, delimiter="\t", lineterminator="\n")
    readlength_writer = csv.writer(readlength_handle, delimiter="\t", lineterminator="\n")
    inventory_writer.writerow(
        [
            "relative_fastq_path",
            "file_size_bytes",
            "gzip_status",
            "first_read_header",
            "records_inspected",
            "read_length_distribution_first100",
        ]
    )
    sha256_writer.writerow(["relative_fastq_path", "file_size_bytes", "sha256"])
    readlength_writer.writerow(
        [
            "relative_fastq_path",
            "read_role",
            "lane",
            "chunk",
            "records_inspected",
            "read_length_distribution_first100",
        ]
    )

    for fastq in fastqs:
        relative_path = fastq.relative_to(output_dir)
        size = fastq.stat().st_size
        total_bytes += size
        with fastq.open("rb") as binary_handle:
            gzip_magic = binary_handle.read(2) == b"\x1f\x8b"
        all_gzip = all_gzip and gzip_magic

        lines: list[str] = []
        with gzip.open(fastq, "rt", encoding="utf-8", errors="strict", newline="") as input_handle:
            for _ in range(400):
                line = input_handle.readline()
                if not line:
                    break
                lines.append(line)
        complete_line_count = len(lines) - (len(lines) % 4)
        lines = lines[:complete_line_count]
        records_inspected = complete_line_count // 4
        if records_inspected == 0:
            raise SystemExit(f"No complete FASTQ records in {relative_path.as_posix()}")

        lengths: dict[int, int] = {}
        for sequence_line in lines[1::4]:
            length = len(sequence_line.rstrip("\r\n"))
            lengths[length] = lengths.get(length, 0) + 1
        length_distribution = ";".join(
            f"{length}:{lengths[length]}" for length in sorted(lengths)
        )
        first_header = lines[0].rstrip("\r\n").replace("\t", " ")

        match = fastq_pattern.match(fastq.name)
        role = "UNKNOWN"
        lane = ""
        chunk = ""
        if match:
            parsed_files += 1
            role = match.group("role")
            lane = match.group("lane")
            chunk = match.group("chunk")
            role_counts[role] += 1
            if role in pair_keys:
                pair_keys[role].add(
                    f"{relative_path.parent.as_posix()}|{match.group('prefix')}|{chunk}"
                )

        sample_path = sample_root / relative_path
        sample_path = sample_path.with_name(sample_path.name + ".first100.fastq.gz")
        sample_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(sample_path, "wt", encoding="utf-8", newline="") as sample_handle:
            sample_handle.writelines(lines)

        digest = sha256_file(fastq)
        inventory_writer.writerow(
            [
                relative_path.as_posix(),
                size,
                "gzip" if gzip_magic else "not_gzip",
                first_header,
                records_inspected,
                length_distribution,
            ]
        )
        sha256_writer.writerow([relative_path.as_posix(), size, digest])
        readlength_writer.writerow(
            [
                relative_path.as_posix(),
                role,
                lane,
                chunk,
                records_inspected,
                length_distribution,
            ]
        )

missing_r2 = sorted(pair_keys["R1"] - pair_keys["R2"])
missing_r1 = sorted(pair_keys["R2"] - pair_keys["R1"])
pairing_passed = (
    parsed_files == len(fastqs)
    and role_counts["R1"] > 0
    and role_counts["R1"] == role_counts["R2"]
    and not missing_r1
    and not missing_r2
)

with pairing_path.open("w", encoding="utf-8", newline="") as pairing_handle:
    writer = csv.writer(pairing_handle, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "r1_files",
            "r2_files",
            "i1_files",
            "i2_files",
            "all_r1_chunks_have_matching_r2",
            "missing_r1_chunks",
            "missing_r2_chunks",
            "parsed_fastq_files",
            "total_gzip_fastq_count",
            "total_fastq_bytes",
            "all_files_gzip",
            "pairing_audit_passed",
        ]
    )
    writer.writerow(
        [
            role_counts["R1"],
            role_counts["R2"],
            role_counts["I1"],
            role_counts["I2"],
            str(not missing_r1 and not missing_r2),
            len(missing_r1),
            len(missing_r2),
            parsed_files,
            len(fastqs),
            total_bytes,
            str(all_gzip),
            str(pairing_passed and all_gzip),
        ]
    )

with tarfile.open(archive_path, "w:gz") as archive:
    for sample_path in sorted(sample_root.rglob("*")):
        if sample_path.is_file():
            archive.add(sample_path, arcname=sample_path.relative_to(sample_root).as_posix())

if not all_gzip:
    raise SystemExit("At least one FASTQ file did not have the gzip magic header")
if not pairing_passed:
    raise SystemExit(
        f"R1/R2 pairing audit failed: missing_r1={len(missing_r1)} missing_r2={len(missing_r2)}"
    )

