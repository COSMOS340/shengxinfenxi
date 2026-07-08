# PC Request: Core GEO Download For CRLM CAF/ECM Project

Created: 2026-07-08

## Project

`/Volumes/research/crlm_caf_ecm_spatial_niche_20260707`

## Reason

Local download from `ftp.ncbi.nlm.nih.gov` is too slow in the Codex desktop session. `GSE225857_RAW.tar` stayed around 40 KB/s, which would take many hours for the three core files. Local partial file was stopped intentionally and must not be treated as complete input.

## Files Already Present In Project

These two files have been imported into the project from the earlier local CRLM work folder and have matching GEO listing sizes:

| File | Project Path | Size Bytes | SHA256 |
|---|---:|---:|---|
| `GSE178318_barcodes.tsv.gz` | `/Volumes/research/crlm_caf_ecm_spatial_niche_20260707/03_data_raw/single_cell/GSE178318/GSE178318_barcodes.tsv.gz` | 592450 | `8c65e63a380ee85439ecba2ff81a1af232cd169f300c38d5a74e9fcbedfdc8c1` |
| `GSE178318_genes.tsv.gz` | `/Volumes/research/crlm_caf_ecm_spatial_niche_20260707/03_data_raw/single_cell/GSE178318/GSE178318_genes.tsv.gz` | 264505 | `e13e2618e6b6ca963394aa0c8afbbd901524afe74f4b0bd8b1dc846fdd57bebc` |

## Files To Download

Download the following official GEO supplementary files. Use resumable download and write to `.part` first; rename only after download completes.

| Dataset | URL | Target Path | Expected Size Bytes |
|---|---|---|---:|
| `GSE225857` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225857/suppl/GSE225857_RAW.tar` | `/Volumes/research/crlm_caf_ecm_spatial_niche_20260707/03_data_raw/geo/GSE225857_RAW.tar` | 636508160 |
| `GSE178318` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/suppl/GSE178318_matrix.mtx.gz` | `/Volumes/research/crlm_caf_ecm_spatial_niche_20260707/03_data_raw/single_cell/GSE178318/GSE178318_matrix.mtx.gz` | 545942023 |
| `GSE245552` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245552/suppl/GSE245552_RAW.tar` | `/Volumes/research/crlm_caf_ecm_spatial_niche_20260707/03_data_raw/geo/GSE245552_RAW.tar` | 1065748480 |

## GitHub PC Run

If running from this GitHub request folder, use:

```bash
python bioinformatics-reproduction/projects/crlm-caf-ecm-spatial-niche/desktop_exchange/requests/20260708_core_geo_download/scripts/pc_download_crlm_caf_core_geo.py
```

The script reads `inputs/core_geo_urls.tsv`, reuses the already downloaded PC raw files from the 20260706 CRLM core GEO audit when present, downloads missing files into `raw_core_geo/`, and writes:

- `outputs/core_geo_download_status.tsv`
- `outputs/core_geo_file_integrity.tsv`

Do not commit `raw_core_geo/`. Commit only the `outputs/` tables and any small status note.

## Suggested Commands

```bash
PROJECT="/Volumes/research/crlm_caf_ecm_spatial_niche_20260707"
mkdir -p "$PROJECT/03_data_raw/geo" "$PROJECT/03_data_raw/single_cell/GSE178318" "$PROJECT/desktop_exchange/uploads/pc_core_geo_download_20260708"

download_one() {
  url="$1"
  output="$2"
  mkdir -p "$(dirname "$output")"
  curl --fail --location --continue-at - --retry 8 --retry-delay 10 --connect-timeout 60 --output "${output}.part" "$url"
  mv "${output}.part" "$output"
}

download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225857/suppl/GSE225857_RAW.tar" "$PROJECT/03_data_raw/geo/GSE225857_RAW.tar"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/suppl/GSE178318_matrix.mtx.gz" "$PROJECT/03_data_raw/single_cell/GSE178318/GSE178318_matrix.mtx.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245552/suppl/GSE245552_RAW.tar" "$PROJECT/03_data_raw/geo/GSE245552_RAW.tar"

{
  printf "file_path\tsize_bytes\tsha256\n"
  for f in \
    "$PROJECT/03_data_raw/geo/GSE225857_RAW.tar" \
    "$PROJECT/03_data_raw/single_cell/GSE178318/GSE178318_matrix.mtx.gz" \
    "$PROJECT/03_data_raw/geo/GSE245552_RAW.tar"
  do
    printf "%s\t%s\t%s\n" "$f" "$(stat -f '%z' "$f")" "$(shasum -a 256 "$f" | awk '{print $1}')"
  done
} > "$PROJECT/desktop_exchange/uploads/pc_core_geo_download_20260708/download_manifest.tsv"
```

## Return Requirements

Please return:

1. The three completed files at the exact target paths above, or tell Codex the exact PC-side path if direct placement is not possible.
2. `download_manifest.tsv` with file path, size in bytes, and SHA256.
3. Any terminal log showing failed retries or renamed `.part` files.

## Codex Acceptance Checks After PC Upload

Codex will run:

```bash
python3 /Volumes/research/crlm_caf_ecm_spatial_niche_20260707/02_scripts/00_setup/audit_core_data_files.py
tar -tf /Volumes/research/crlm_caf_ecm_spatial_niche_20260707/03_data_raw/geo/GSE225857_RAW.tar | head
tar -tf /Volumes/research/crlm_caf_ecm_spatial_niche_20260707/03_data_raw/geo/GSE245552_RAW.tar | head
gzip -t /Volumes/research/crlm_caf_ecm_spatial_niche_20260707/03_data_raw/single_cell/GSE178318/GSE178318_matrix.mtx.gz
```
