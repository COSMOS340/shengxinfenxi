#!/usr/bin/env python3
from __future__ import annotations

import base64
import csv
import datetime as dt
import hashlib
from pathlib import Path


CHUNK_CHARS = 800_000

FIGURES = [
    {
        "file_id": "supp_fig_s2_full_cell_composition",
        "priority": "required",
        "source_path": r"I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_supp_fig_s2_full_cell_composition\pc_supp_fig_s2_full_cell_composition_results\outputs\figures\supp_fig_s2_full_cell_composition.png",
        "target_path": "figures/supp_fig_s2_full_cell_composition.png",
        "description": "Supplementary Fig. S2 full six-dataset PC figure",
    },
    {
        "file_id": "comparison_original_supp_fig_s2_vs_pc_full",
        "priority": "required",
        "source_path": r"I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_supp_fig_s2_full_cell_composition\pc_supp_fig_s2_full_cell_composition_results\outputs\figures\comparison_original_supp_fig_s2_vs_pc_full.png",
        "target_path": "figures/comparison_original_supp_fig_s2_vs_pc_full.png",
        "description": "PC-side S2 comparison image, if present in the PC result package",
    },
    {
        "file_id": "supp_fig_s4_s5_spatial_featureplots",
        "priority": "required",
        "source_path": r"I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_supp_fig_s4_s5_spatial_featureplots\pc_supp_fig_s4_s5_spatial_featureplots_results\figures\supp_fig_s4_s5_spatial_featureplots.png",
        "target_path": "figures/supp_fig_s4_s5_spatial_featureplots.png",
        "description": "Supplementary Fig. S4/S5 PC spatial feature plots",
    },
    {
        "file_id": "supp_fig_s6_spatial_signatures",
        "priority": "required",
        "source_path": r"I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_supp_fig_s6_spatial_signatures\pc_supp_fig_s6_spatial_signatures_results\figures\supp_fig_s6_spatial_signatures.png",
        "target_path": "figures/supp_fig_s6_spatial_signatures.png",
        "description": "Supplementary Fig. S6 PC spatial signature and correlation figure",
    },
    {
        "file_id": "fig5e_full_seurat_proportions",
        "priority": "optional",
        "source_path": r"I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-cancer-mono-mac\desktop_exchange\uploads\20260704_fig5e_gse120575_full_seurat\pc_fig5e_full_seurat_results\outputs\figures\fig5e_full_seurat_proportions.png",
        "target_path": "figures/fig5e_full_seurat_proportions.png",
        "description": "Optional Fig. 5E all-zero audit PNG from the completed PC route",
    },
]


REASSEMBLE_SCRIPT = r'''#!/usr/bin/env python3
from __future__ import annotations

import base64
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "upload_manifest.tsv"
CHUNK_DIR = ROOT / "files" / "base64_chunks"
OUT_DIR = ROOT / "reassembled"


def sha256_bytes(data: bytes) -> str:
    h = hashlib.sha256()
    h.update(data)
    return h.hexdigest()


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["status"] != "packaged":
                continue
            file_id = row["file_id"]
            chunk_count = int(row["chunk_count"])
            pieces = []
            for idx in range(1, chunk_count + 1):
                chunk_path = CHUNK_DIR / f"{file_id}.part{idx:04d}.b64.txt"
                pieces.append(chunk_path.read_text(encoding="ascii").strip())
            data = base64.b64decode("".join(pieces).encode("ascii"))
            observed = sha256_bytes(data)
            expected = row["sha256"]
            if observed != expected:
                raise SystemExit(f"SHA256 mismatch for {file_id}: {observed} != {expected}")
            out_path = OUT_DIR / row["target_path"]
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_bytes(data)
            print(f"reassembled {out_path} ({len(data)} bytes)")


if __name__ == "__main__":
    main()
'''


def find_project_root(start: Path) -> Path:
    for path in [start, *start.parents]:
        if path.name == "pan-cancer-mono-mac" and (path / "desktop_exchange").is_dir():
            return path
    raise RuntimeError("Could not find pan-cancer-mono-mac project root from script path")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    script_path = Path(__file__).resolve()
    project_root = find_project_root(script_path)
    out_dir = project_root / "desktop_exchange" / "uploads" / "20260705_pc_binary_figure_upload" / "pc_binary_figure_text_handoff"
    chunk_dir = out_dir / "files" / "base64_chunks"
    direct_dir = out_dir / "figures_direct"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    direct_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []

    for spec in FIGURES:
        source = Path(spec["source_path"])
        row = {
            "file_id": spec["file_id"],
            "priority": spec["priority"],
            "source_path": str(source),
            "target_path": spec["target_path"],
            "description": spec["description"],
            "status": "missing",
            "bytes": 0,
            "sha256": "",
            "base64_chars": 0,
            "chunk_count": 0,
            "chunk_prefix": f"files/base64_chunks/{spec['file_id']}.part",
        }
        if not source.exists():
            manifest_rows.append(row)
            missing_rows.append({
                "file_id": spec["file_id"],
                "priority": spec["priority"],
                "source_path": str(source),
                "target_path": spec["target_path"],
                "description": spec["description"],
                "reason": "source_path_not_found",
            })
            continue

        data = source.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        b64 = base64.b64encode(data).decode("ascii")
        chunks = [b64[i:i + CHUNK_CHARS] for i in range(0, len(b64), CHUNK_CHARS)]
        for idx, chunk in enumerate(chunks, start=1):
            chunk_path = chunk_dir / f"{spec['file_id']}.part{idx:04d}.b64.txt"
            chunk_path.write_text(chunk + "\n", encoding="ascii")

        direct_path = direct_dir / Path(spec["target_path"]).name
        direct_path.write_bytes(data)

        row.update({
            "status": "packaged",
            "bytes": len(data),
            "sha256": digest,
            "base64_chars": len(b64),
            "chunk_count": len(chunks),
        })
        manifest_rows.append(row)

    write_tsv(
        out_dir / "upload_manifest.tsv",
        manifest_rows,
        [
            "file_id",
            "priority",
            "source_path",
            "target_path",
            "description",
            "status",
            "bytes",
            "sha256",
            "base64_chars",
            "chunk_count",
            "chunk_prefix",
        ],
    )
    write_tsv(
        out_dir / "MISSING_FILES.tsv",
        missing_rows,
        ["file_id", "priority", "source_path", "target_path", "description", "reason"],
    )

    required_total = sum(1 for item in FIGURES if item["priority"] == "required")
    required_found = sum(1 for row in manifest_rows if row["priority"] == "required" and row["status"] == "packaged")
    optional_found = sum(1 for row in manifest_rows if row["priority"] == "optional" and row["status"] == "packaged")
    now = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    status = "complete" if required_found == required_total else "partial"
    status_text = f"""# PC binary figure text handoff

Status: {status}
Created: {now}

required_found: {required_found}
required_total: {required_total}
optional_found: {optional_found}

This package converts PNG figures into base64 text chunks so the GitHub connector can upload them as text files.

After this directory is uploaded to GitHub, run:

```bash
python3 reassemble_binary_figures.py
```

from this handoff directory to reconstruct the PNG files under `reassembled/`.
"""
    (out_dir / "STATUS.md").write_text(status_text, encoding="utf-8")
    (out_dir / "reassemble_binary_figures.py").write_text(REASSEMBLE_SCRIPT, encoding="utf-8")

    print(status_text)
    print(f"Output directory: {out_dir}")


if __name__ == "__main__":
    main()
