#!/usr/bin/env bash
set +e

audit_dir=/work/audit
mkdir -p "$audit_dir"
echo "$(date --iso-8601=seconds)" > "$audit_dir/pilot_wrapper_start.txt"
bash /work/tools/run_pilot.sh > "$audit_dir/pilot_wrapper.stdout.txt" 2> "$audit_dir/pilot_wrapper.stderr.txt"
exit_code=$?
echo "$exit_code" > "$audit_dir/pilot.exit_code"
echo "$(date --iso-8601=seconds)" > "$audit_dir/pilot_wrapper_end.txt"
exit "$exit_code"
