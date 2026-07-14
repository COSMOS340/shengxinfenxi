#!/usr/bin/env bash
set +e

audit_dir=/work/audit
mkdir -p "$audit_dir"
echo "$(date --iso-8601=seconds)" > "$audit_dir/full_wrapper_start.txt"
bash /work/tools/run_full.sh > "$audit_dir/full_wrapper.stdout.txt" 2> "$audit_dir/full_wrapper.stderr.txt"
exit_code=$?
echo "$exit_code" > "$audit_dir/full.exit_code"
echo "$(date --iso-8601=seconds)" > "$audit_dir/full_wrapper_end.txt"
exit "$exit_code"
