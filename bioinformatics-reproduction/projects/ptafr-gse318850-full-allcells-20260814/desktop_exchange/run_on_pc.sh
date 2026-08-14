#!/usr/bin/env bash
set -euo pipefail

project_root="${PTAFR_PROJECT_ROOT:-/mnt/i/ptafr/05_公共数据库复现_图ABCD}"
exchange_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_dir="$project_root/04_data_processed/single_cell/GSE318850/reanalysis_full"
log_dir="$project_root/99_logs"
manifest="$exchange_dir/input_manifest.tsv"

counts_h5="$project_root/04_data_processed/single_cell/GSE318850/GSE318850_RNA_counts_10x.h5"
source_metadata="$project_root/04_data_processed/single_cell/GSE318850/reannotation/GSE318850_full_independent_cell_metadata.rds"
qc_doublet_metadata="$project_root/04_data_processed/single_cell/GSE318850/qc_doublet/GSE318850_qc_doublet_metadata.rds"
analysis_script="$exchange_dir/10_GSE318850_full_harmony_all_cells.R"
resource_script="$exchange_dir/detect_pc_resources.py"

for required_file in \
  "$counts_h5" \
  "$source_metadata" \
  "$qc_doublet_metadata" \
  "$analysis_script" \
  "$resource_script" \
  "$manifest"; do
  if [[ ! -f "$required_file" ]]; then
    printf 'Required file is missing: %s\n' "$required_file" >&2
    exit 1
  fi
done

while IFS=$'\t' read -r role relative_path_under_project_root expected_bytes expected_sha256; do
  if [[ "$role" == "role" ]]; then
    continue
  fi
  input_file="$project_root/$relative_path_under_project_root"
  actual_bytes="$(stat -c '%s' "$input_file")"
  actual_sha256="$(sha256sum "$input_file" | cut -d ' ' -f 1)"
  if [[ "$actual_bytes" != "$expected_bytes" || "$actual_sha256" != "$expected_sha256" ]]; then
    printf 'Input verification failed: %s\n' "$role" >&2
    exit 1
  fi
  printf 'Input verified: %s\n' "$role"
done < "$manifest"

mkdir -p "$output_dir" "$log_dir"
timestamp="$(date +%Y%m%d_%H%M%S)"
resource_json="$output_dir/PC_resource_snapshot_${timestamp}.json"
session_info="$output_dir/PC_R_sessionInfo_${timestamp}.txt"
run_log="$log_dir/GSE318850_full_harmony_all_cells_PC_${timestamp}.log"

analysis_script_win="$(wslpath -w "$analysis_script")"
resource_script_win="$(wslpath -w "$resource_script")"
project_root_win="$(wslpath -w "$project_root")"
resource_json_win="$(wslpath -w "$resource_json")"
session_info_win="$(wslpath -w "$session_info")"
counts_h5_win="$(wslpath -w "$counts_h5")"
source_metadata_win="$(wslpath -w "$source_metadata")"
qc_doublet_metadata_win="$(wslpath -w "$qc_doublet_metadata")"
output_dir_win="$(wslpath -w "$output_dir")"

/mnt/h/python312/python.exe "$resource_script_win" \
  --output "$resource_json_win" \
  --workdir "$project_root_win"

/mnt/h/R/bin/Rscript.exe -e \
  "args <- commandArgs(trailingOnly=TRUE); required <- c('BPCells','data.table','harmony','Matrix','MatrixGenerics','Seurat'); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop('Missing R packages: ', paste(missing, collapse=', ')); writeLines(capture.output(sessionInfo()), args[[1]])" \
  "$session_info_win"

printf 'Starting all-cell GSE318850 analysis. Log: %s\n' "$run_log"
/mnt/h/R/bin/Rscript.exe "$analysis_script_win" \
  "$counts_h5_win" \
  "$source_metadata_win" \
  "$qc_doublet_metadata_win" \
  "$output_dir_win" \
  2>&1 | tee "$run_log"
