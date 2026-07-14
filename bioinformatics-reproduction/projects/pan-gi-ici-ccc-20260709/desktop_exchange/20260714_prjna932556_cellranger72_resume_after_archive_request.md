# PC request: resume SRR23490337 Cell Ranger 7.2.0 count after archive acquisition

Date: 2026-07-14

Project: pan-GI ICI cell communication atlas

## Objective

Resume the previously authorized one-sample Cell Ranger hard run now that the exact Cell Ranger 7.2.0 Linux archive is present on `F:`. Continue through archive verification, installation, reference download and installation, QEMU setup, and one complete `cellranger count` attempt.

The prior manual-download block is resolved. Do not stop with an archive-download status and do not redownload Cell Ranger.

## Parent request and preserved boundaries

Follow all scientific, resource, storage, logging, and stop boundaries in:

```text
20260714_prjna932556_cellranger72_single_sample_hardrun_request.md
```

This resume request changes only the resolved Cell Ranger archive precondition and the upload directory. It does not authorize another sample, another Cell Ranger version, another reference, another quantifier, or downstream analysis.

The single authorized count attempt remains unused because the prior result confirms that `cellranger count` was not executed.

## Verified Cell Ranger archive

Use exactly this PC-local archive:

```text
F:\cellranger-7.2.0.tar.gz
```

Expected audit values:

- filename: `cellranger-7.2.0.tar.gz`
- bytes: `683925475`
- official MD5: `85e2573e80a6f8656a42ff09460463e6`
- local SHA256: `b092bd4e3ab585ad051a231fbdd8f3f0f5cbcd10f657eeab86bec98cd594502c`
- official source page: `https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions`

Before extraction, independently recompute file bytes, MD5, and SHA256. Stop with `CELLRANGER72_COUNT_FAILED_OTHER` if any value differs. Preserve the file and report the exact observation; do not repair, rename, or replace it in the same task.

## Required execution order

1. Re-audit host resources and current free space on `F:` and `I:`. Do not write task files to `E:`.
2. Recompute and verify the Cell Ranger archive byte count, MD5, and SHA256.
3. Use the Linux archive only inside the QEMU Ubuntu workflow. Store the extracted software, reference, QEMU data disk, temporary files, pipestance, and final results physically on `F:`.
4. Record the extracted executable path, bytes, SHA256, `cellranger --version`, and `cellranger count --help` exit code. Continue only if the observed version is exactly Cell Ranger `7.2.0`.
5. Download and install only `refdata-gex-GRCh38-2020-A` from:

```text
https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

6. Verify the existing `SRR23490337` FASTQ directory and use sample prefix `bamtofastq` exactly as specified in the parent request.
7. Start the QEMU guest with the parent request's resource policy: 6 vCPU and 24 GiB first, with one VM-start retry at 20 GiB only if required. Record acceleration and observed guest resources.
8. Run the exact parent-request `cellranger count` command once. Do not add `--expect-cells`, `--force-cells`, downsampling, trimming, or unverified flags.
9. Monitor chemistry, memory, and free space. Preserve the pipestance and logs after either success or failure.
10. Upload all required lightweight audit files and summaries to the new upload directory below. Keep large data PC-local.
11. Finalize every report file before creating `upload_manifest.tsv` and `SHA256SUMS.txt`. Recompute the checksum of each uploaded file from its final bytes and verify the complete checksum list before committing.

## Inputs that must remain unchanged

- run: `SRR23490337`
- sample: `R3`
- FASTQ directory:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260714_prjna932556_srr23490337_full_bamtofastq\W6B_0_1_HGGWHDSX2
```

- expected gzip FASTQ files: `804`
- expected total FASTQ bytes: `29268262419`
- expected R1/R2 pairs: `402` / `402`
- required Cell Ranger: `7.2.0`
- required reference: `refdata-gex-GRCh38-2020-A`
- expected detected chemistry: `SC3Pv2`

## Required upload directory

Upload lightweight outputs to:

```text
bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260714_prjna932556_cellranger72_resume_after_archive
```

Required files:

- `STATUS.md`
- `archive_reverification.tsv`
- `fastq_reverification_summary.tsv`
- `host_resource_before_run.tsv`
- `qemu_acceleration_audit.txt`
- `guest_resource_allocation.tsv`
- `software_install_audit.tsv`
- `reference_install_audit.tsv`
- `disk_component_sizes.tsv`
- `disk_usage_timeseries.tsv`
- `cellranger_count_command.txt`
- `cellranger_count_runtime.tsv`
- `cellranger_count_exit_code.txt`
- `chemistry_detection.txt`
- `cellranger_count_log_tail.txt`
- `cellranger_failure_logs.tar.gz` if the run fails and the archive remains below 100 MB
- `cellranger_outs_inventory.tsv` if `outs` is created
- `metrics_summary.csv` if generated
- `web_summary.html` if generated
- `upload_manifest.tsv`
- `SHA256SUMS.txt`

Do not upload the Cell Ranger archive, reference archive, extracted reference, FASTQs, BAMs, QEMU images, molecule information, HDF5 matrices, matrix directories, or complete pipestance.

## Final status

Use exactly one final status:

- `CELLRANGER72_COUNT_SUCCESS`
- `CELLRANGER72_COUNT_FAILED_MEMORY`
- `CELLRANGER72_COUNT_FAILED_DISK`
- `CELLRANGER72_COUNT_FAILED_CHEMISTRY`
- `CELLRANGER72_COUNT_FAILED_OTHER`
- `BLOCKED_REFERENCE_DOWNLOAD`
- `BLOCKED_QEMU_SETUP`

`CELLRANGER72_ARCHIVE_DOWNLOADED_TO_F` and `BLOCKED_CELLRANGER72_MANUAL_DOWNLOAD` are not valid final statuses for this resumed task because the archive has already been acquired and verified.

Success requires exit code 0, an observed `outs` directory, `metrics_summary.csv`, `web_summary.html`, a complete output inventory, and complete lightweight checksums.
