#!/usr/bin/env python3
import argparse
import csv
import math
import os
import sys
import time
import urllib.error
import urllib.request


def parse_args():
    parser = argparse.ArgumentParser(description="Extract selected uncompressed tar members using HTTP Range requests.")
    parser.add_argument("--url", required=True)
    parser.add_argument("--filelist", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--file-type", default="H5")
    parser.add_argument("--chunk-mb", type=int, default=64)
    parser.add_argument("--retries", type=int, default=20)
    parser.add_argument("--manifest", required=True)
    return parser.parse_args()


def request_range(url, start, end, retries):
    headers = {"Range": f"bytes={start}-{end}", "User-Agent": "codex-range-extractor/1.0"}
    last_error = None
    for attempt in range(1, retries + 1):
        req = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=120) as response:
                status = getattr(response, "status", None)
                if status not in (200, 206):
                    raise RuntimeError(f"unexpected HTTP status {status}")
                data = response.read()
            expected = end - start + 1
            if len(data) != expected:
                raise RuntimeError(f"range length mismatch: observed {len(data)}, expected {expected}")
            return data
        except (urllib.error.URLError, TimeoutError, RuntimeError) as exc:
            last_error = exc
            sleep_seconds = min(60, 2 * attempt)
            print(f"RANGE_RETRY\t{start}\t{end}\tattempt={attempt}\terror={exc}", flush=True)
            time.sleep(sleep_seconds)
    raise RuntimeError(f"range request failed after {retries} attempts: {start}-{end}: {last_error}")


def parse_tar_header(block):
    if block == b"\0" * 512:
        return None
    name = block[0:100].split(b"\0", 1)[0].decode("utf-8", errors="replace")
    size_text = block[124:136].split(b"\0", 1)[0].strip() or b"0"
    size = int(size_text, 8)
    typeflag = block[156:157].decode("ascii", errors="replace") or "0"
    prefix = block[345:500].split(b"\0", 1)[0].decode("utf-8", errors="replace")
    full_name = f"{prefix}/{name}" if prefix else name
    return {"name": full_name, "size": size, "typeflag": typeflag}


def read_targets(filelist, file_type):
    targets = {}
    archive_url = None
    with open(filelist, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["record_type"] == "Archive" and row["file_name"] == "GSE306459_RAW.tar":
                archive_url = row["url"]
            if row["record_type"] == "File" and row["file_type"] == file_type:
                targets[row["file_name"]] = int(row["size_bytes"])
    if not targets:
        raise RuntimeError(f"no targets found for file_type={file_type}")
    return archive_url, targets


def download_member(url, data_offset, size, out_path, retries, chunk_bytes):
    existing = os.path.getsize(out_path) if os.path.exists(out_path) else 0
    if existing == size:
        return "skipped_existing"
    if existing > size:
        os.remove(out_path)
        existing = 0

    mode = "ab" if existing else "wb"
    written = existing
    with open(out_path, mode) as out_handle:
        while written < size:
            chunk_start = data_offset + written
            take = min(chunk_bytes, size - written)
            chunk_end = chunk_start + take - 1
            data = request_range(url, chunk_start, chunk_end, retries)
            out_handle.write(data)
            written += len(data)
            out_handle.flush()
            print(f"DOWNLOAD_PROGRESS\t{os.path.basename(out_path)}\t{written}/{size}", flush=True)

    observed = os.path.getsize(out_path)
    if observed != size:
        raise RuntimeError(f"downloaded size mismatch for {out_path}: observed {observed}, expected {size}")
    return "downloaded"


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    archive_url, targets = read_targets(args.filelist, args.file_type)
    if archive_url and archive_url != args.url:
        print(f"INFO\tarchive_url_from_filelist={archive_url}", flush=True)

    chunk_bytes = args.chunk_mb * 1024 * 1024
    found = {}
    offset = 0

    while len(found) < len(targets):
        block = request_range(args.url, offset, offset + 511, args.retries)
        header = parse_tar_header(block)
        if header is None:
            break

        name = header["name"]
        basename = os.path.basename(name)
        size = header["size"]
        data_offset = offset + 512
        padded_size = int(math.ceil(size / 512.0) * 512) if size else 0

        print(f"TAR_MEMBER\t{offset}\t{name}\t{size}\t{header['typeflag']}", flush=True)
        if basename in targets:
            expected_size = targets[basename]
            if size != expected_size:
                raise RuntimeError(f"tar member size mismatch for {basename}: observed {size}, expected {expected_size}")
            out_path = os.path.join(args.out_dir, basename)
            status = download_member(args.url, data_offset, size, out_path, args.retries, chunk_bytes)
            found[basename] = {
                "tar_name": name,
                "data_offset": data_offset,
                "size_bytes": size,
                "local_path": out_path,
                "status": status,
            }

        offset = data_offset + padded_size

    missing = sorted(set(targets) - set(found))
    if missing:
        raise RuntimeError("missing target members: " + ", ".join(missing))

    with open(args.manifest, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["file_name", "tar_name", "data_offset", "size_bytes", "local_path", "status"],
            delimiter="\t",
        )
        writer.writeheader()
        for file_name in sorted(found):
            row = {"file_name": file_name}
            row.update(found[file_name])
            writer.writerow(row)

    print(f"COMPLETE\tfiles={len(found)}\tmanifest={args.manifest}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR\t{exc}", file=sys.stderr, flush=True)
        sys.exit(1)
