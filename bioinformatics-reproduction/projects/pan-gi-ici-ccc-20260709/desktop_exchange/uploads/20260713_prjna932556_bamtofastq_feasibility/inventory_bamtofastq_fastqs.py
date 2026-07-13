from __future__ import annotations

import csv
import gzip
from pathlib import Path
import shutil
import sys
import tarfile


output_dir = Path(sys.argv[1])
audit_dir = Path(sys.argv[2])
inventory_path = audit_dir / "bamtofastq_fastq_inventory.tsv"
sample_root = audit_dir / "first100_records"
archive_path = audit_dir / "bamtofastq_first100_records.tar.gz"
sample_root.mkdir(parents=True, exist_ok=True)


def is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_text(path: Path, compressed: bool):
    if compressed:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("rt", encoding="utf-8", errors="replace", newline="")


fastqs = sorted(
    path
    for path in output_dir.rglob("*")
    if path.is_file() and ("fastq" in path.name.lower() or path.suffix.lower() in {".fq", ".gz"})
)

with inventory_path.open("w", encoding="utf-8", newline="") as inventory_handle:
    writer = csv.writer(inventory_handle, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "relative_fastq_path",
            "file_size_bytes",
            "gzip_status",
            "first_read_header",
            "records_inspected",
            "read_length_distribution_first100",
        ]
    )

    for fastq in fastqs:
        relative_path = fastq.relative_to(output_dir)
        compressed = is_gzip(fastq)
        lines: list[str] = []
        with open_text(fastq, compressed) as input_handle:
            for _ in range(400):
                line = input_handle.readline()
                if not line:
                    break
                lines.append(line)

        complete_line_count = len(lines) - (len(lines) % 4)
        lines = lines[:complete_line_count]
        records_inspected = complete_line_count // 4
        lengths: dict[int, int] = {}
        for sequence_line in lines[1::4]:
            length = len(sequence_line.rstrip("\r\n"))
            lengths[length] = lengths.get(length, 0) + 1
        length_distribution = ";".join(f"{length}:{lengths[length]}" for length in sorted(lengths))
        first_header = lines[0].rstrip("\r\n").replace("\t", " ") if lines else ""

        sample_path = sample_root / relative_path
        sample_path = sample_path.with_name(sample_path.name + ".first100.fastq.gz")
        sample_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(sample_path, "wt", encoding="utf-8", newline="") as sample_handle:
            sample_handle.writelines(lines)

        writer.writerow(
            [
                relative_path.as_posix(),
                fastq.stat().st_size,
                "gzip" if compressed else "plain",
                first_header,
                records_inspected,
                length_distribution,
            ]
        )

if fastqs:
    with tarfile.open(archive_path, "w:gz") as archive:
        for sample_path in sorted(sample_root.rglob("*")):
            if sample_path.is_file():
                archive.add(sample_path, arcname=sample_path.relative_to(sample_root).as_posix())
else:
    shutil.rmtree(sample_root)
    raise SystemExit("No FASTQ files were produced")
