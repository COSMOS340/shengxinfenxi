#!/usr/bin/env bash
set -uo pipefail

input_dir=/work/input
output_dir=/work/output
audit_dir=/work/audit
tool_dir=/work/tools
bam_name=R3_possorted_genome_bam.bam.1
bai_name=R3_possorted_genome_bam.bam.1.bai
bam_path="$input_dir/$bam_name"
bai_path="$input_dir/$bai_name"
bamtofastq_path="$tool_dir/bamtofastq"
inventory_script="$tool_dir/inventory_fastq.py"
host_source=http://10.0.2.2:18000
official_bamtofastq_url=https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux
expected_tool_sha256=fdf7db4fe6cf8e13a432e8a59815dad506415349627502fdf77f4be208a198ce
expected_bam_md5=056a5c996212a64af708d48c5b585e48
expected_bai_md5=f391d287338dbf9d3856630338957bb8

mkdir -p "$input_dir" "$audit_dir" "$tool_dir"

download_file() {
    local source_url=$1
    local destination=$2
    local expected_bytes=${3:-}
    if [[ -n "$expected_bytes" && -f "$destination" ]]; then
        local actual_bytes
        actual_bytes=$(stat -c %s "$destination")
        if [[ "$actual_bytes" == "$expected_bytes" ]]; then
            echo "Reusing complete local file: $destination" >&2
            return 0
        fi
        if (( actual_bytes > expected_bytes )); then
            rm -f "$destination"
        fi
    elif [[ -f "$destination" ]]; then
        rm -f "$destination"
    fi
    curl --fail --location --retry 8 --retry-delay 5 --continue-at - \
        --output "$destination" "$source_url"
}

download_file "$host_source/$bam_name" "$bam_path" 21494429186
download_file "$host_source/$bai_name" "$bai_path" 10390008
download_file "$host_source/bamtofastq_linux" "$bamtofastq_path" 18540256
download_file "$host_source/inventory_fastq.py" "$inventory_script"
chmod 700 "$bamtofastq_path"

actual_bam_md5=$(md5sum "$bam_path" | awk '{print $1}')
actual_bai_md5=$(md5sum "$bai_path" | awk '{print $1}')
actual_tool_sha256=$(sha256sum "$bamtofastq_path" | awk '{print $1}')
if [[ "$actual_bam_md5" != "$expected_bam_md5" ]]; then
    echo "BAM MD5 mismatch: $actual_bam_md5" >&2
    exit 11
fi
if [[ "$actual_bai_md5" != "$expected_bai_md5" ]]; then
    echo "BAI MD5 mismatch: $actual_bai_md5" >&2
    exit 12
fi
if [[ "$actual_tool_sha256" != "$expected_tool_sha256" ]]; then
    echo "bamtofastq SHA256 mismatch: $actual_tool_sha256" >&2
    exit 13
fi

{
    echo "LINUX_AUDIT_START: $(date --iso-8601=seconds)"
    echo "OFFICIAL_BAMTOFASTQ_URL: $official_bamtofastq_url"
    echo "LOCAL_TRANSFER_URL: $host_source/bamtofastq_linux"
    echo "BAMTOFASTQ_PATH: $bamtofastq_path"
    echo "BAMTOFASTQ_BYTES: $(stat -c %s "$bamtofastq_path")"
    echo "BAMTOFASTQ_SHA256: $actual_tool_sha256"
    echo
    echo 'COMMAND: uname -a'
    uname -a
    echo
    echo 'COMMAND: pwd'
    pwd
    echo
    echo 'COMMAND: which bash'
    command -v bash || true
    echo
    echo 'COMMAND: which samtools'
    command -v samtools || true
    echo
    echo 'COMMAND: samtools --version'
    samtools --version || true
    echo
    echo 'COMMAND: which bamtofastq'
    echo "$bamtofastq_path"
    echo
    echo 'COMMAND: bamtofastq --version'
    "$bamtofastq_path" --version || true
    echo
    echo 'COMMAND: bamtofastq --help'
    "$bamtofastq_path" --help || true
    echo
    echo "LINUX_AUDIT_END: $(date --iso-8601=seconds)"
} > "$audit_dir/runtime_linux_audit.txt" 2>&1

set +e
samtools quickcheck -v "$bam_path" > "$audit_dir/samtools_quickcheck.stdout.txt" 2> "$audit_dir/samtools_quickcheck.stderr.txt"
quickcheck_exit=$?
set -e
if [[ $quickcheck_exit -ne 0 ]]; then
    echo "samtools quickcheck failed with exit code $quickcheck_exit" >&2
    exit 21
fi

samtools view -H "$bam_path" > "$audit_dir/bam_header.sam"
samtools idxstats "$bam_path" > "$audit_dir/bam_idxstats.tsv"

if awk -F '\t' '$1 == "chrM" { found=1 } END { exit !found }' "$audit_dir/bam_idxstats.tsv"; then
    validated_locus=chrM:1-16569
elif awk -F '\t' '$1 == "MT" { found=1 } END { exit !found }' "$audit_dir/bam_idxstats.tsv"; then
    validated_locus=MT:1-16569
else
    echo 'Neither exact contig chrM nor exact contig MT was found.' >&2
    exit 22
fi

command_text="$bamtofastq_path --locus=$validated_locus --reads-per-fastq=1000000 $bam_path $output_dir"
if [[ -e "$output_dir" ]]; then
    if [[ -d "$output_dir" && -z "$(find "$output_dir" -mindepth 1 -print -quit)" ]]; then
        rmdir "$output_dir"
    else
        echo "Refusing to replace non-empty output path: $output_dir" >&2
        exit 23
    fi
fi
start_time=$(date --iso-8601=seconds)
disk_log="$audit_dir/bamtofastq_disk_footprint.tsv"
printf 'timestamp\tbytes\n' > "$disk_log"

set +e
/usr/bin/time -v "$bamtofastq_path" \
    "--locus=$validated_locus" \
    --reads-per-fastq=1000000 \
    "$bam_path" "$output_dir" \
    > "$audit_dir/bamtofastq.stdout.txt" \
    2> "$audit_dir/bamtofastq.stderr.txt" &
pilot_pid=$!
while kill -0 "$pilot_pid" 2>/dev/null; do
    printf '%s\t%s\n' "$(date --iso-8601=seconds)" "$(du -sb "$output_dir" | awk '{print $1}')" >> "$disk_log"
    sleep 5
done
wait "$pilot_pid"
pilot_exit=$?
set -e

end_time=$(date --iso-8601=seconds)
printf '%s\t%s\n' "$(date --iso-8601=seconds)" "$(du -sb "$output_dir" | awk '{print $1}')" >> "$disk_log"
peak_disk_bytes=$(awk -F '\t' 'NR > 1 && $2 > max { max=$2 } END { print max+0 }' "$disk_log")

{
    echo "COMMAND: $command_text"
    echo "VALIDATED_LOCUS: $validated_locus"
    echo "START_TIME: $start_time"
    echo "END_TIME: $end_time"
    echo "EXIT_CODE: $pilot_exit"
    echo "PEAK_OUTPUT_DISK_BYTES: $peak_disk_bytes"
    echo
    echo 'STDOUT_BEGIN'
    cat "$audit_dir/bamtofastq.stdout.txt"
    echo 'STDOUT_END'
    echo
    echo 'STDERR_AND_TIME_BEGIN'
    cat "$audit_dir/bamtofastq.stderr.txt"
    echo 'STDERR_AND_TIME_END'
} > "$audit_dir/bamtofastq_command_log.txt"

if [[ $pilot_exit -ne 0 ]]; then
    exit "$pilot_exit"
fi

python3 "$inventory_script" "$output_dir" "$audit_dir"
