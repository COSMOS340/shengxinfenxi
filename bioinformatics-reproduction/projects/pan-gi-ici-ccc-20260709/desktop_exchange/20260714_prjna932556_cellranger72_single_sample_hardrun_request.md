# PC request: SRR23490337 Cell Ranger 7.2.0 single-sample hard run

Date: 2026-07-14

Project: pan-GI ICI cell communication atlas

## Authorization and purpose

The user explicitly authorized one PC-side attempt to run `cellranger count` for `SRR23490337` / sample `R3` with the resources actually available on the PC.

This request supersedes the previous fixed `500000000000`-byte free-space gate for this one attempt. Disk use must be based on measured files and live free space. The request does not relax input integrity, software version, reference version, sample identity, logging, or one-sample boundaries.

Do not process the other five `PRJNA932556` runs.

## Verified inputs

Use the existing PC-local FASTQs:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260714_prjna932556_srr23490337_full_bamtofastq\W6B_0_1_HGGWHDSX2
```

Verified state:

- run: `SRR23490337`
- sample: `R3`
- response: resistant / progression
- gzip FASTQ files: `804`
- R1 files: `402`
- R2 files: `402`
- total FASTQ bytes: `29268262419`
- read pairs: `401831188`
- FASTQ prefix: `bamtofastq`
- R1 length: `28`
- R2 length: `91`
- reported assay: Chromium Single Cell 3' Reagent Kit v2
- full FASTQ SHA256 mismatches: `0`

Before any installation or count work, rerun the existing FASTQ checksum and pairing verification. Stop if any mismatch is observed.

## Exact software and reference

Use only:

- Cell Ranger `7.2.0`
- `refdata-gex-GRCh38-2020-A`
- official reference URL: `https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz`
- official Cell Ranger download page: `https://www.10xgenomics.com/support/software/cell-ranger/downloads`

Do not substitute another Cell Ranger release, another genome build, STARsolo, kallisto, bustools, salmon, or any alternative quantifier.

The Cell Ranger archive requires the official 10x download flow. Do not invent an archive URL, identity, email address, institution, or EULA response. Use a legitimately obtained Cell Ranger 7.2.0 archive if it is already available or can be downloaded through the user's existing browser session. If manual user input is required, upload the audit with status `BLOCKED_CELLRANGER72_MANUAL_DOWNLOAD` and stop.

After installation, record:

- exact archive path, bytes, and SHA256
- exact extracted executable path, bytes, and SHA256
- output of `cellranger --version`
- output and exit code of `cellranger count --help`

Continue only when the observed version is exactly `cellranger cellranger-7.2.0` or the exact equivalent emitted by the executable.

Download the reference from the exact URL above. Record archive bytes and SHA256, extraction path, extracted bytes, and reference metadata files. Do not rename the extracted reference directory.

## PC resource allocation

The previous audit observed:

- host CPU: AMD Ryzen 5 5600, 6 physical cores, 12 logical processors
- host RAM: `34282172416` bytes
- host AVX and AVX2: available
- `E:` free bytes: `417535246336`
- `F:` free bytes: `393195249664`
- existing QEMU guest allocation: 4 vCPU and approximately 4 GB RAM

`E:` is the user's application volume and is not authorized for this task. Do not write Cell Ranger files, reference files, QEMU data disks, temporary files, pipestance files, or results to `E:` even if it has more free space.

Re-audit these values immediately before setup. Use the existing QEMU Ubuntu environment because WSL and Docker were not installed in the previous audit.

For the count attempt:

1. Close nonessential PC applications before starting the guest. Do not terminate unrelated research jobs without checking them first.
2. Query `qemu-system-x86_64.exe -accel help`. Use `whpx` only if the installed QEMU reports it as available and a short boot test succeeds. Otherwise use the existing QEMU execution mode and record that acceleration was unavailable.
3. First try a guest allocation of 6 vCPU and 24 GiB RAM. If the guest cannot start because the host cannot allocate 24 GiB, make one VM-start retry at 20 GiB. Do not make repeated allocation retries.
4. Inside the guest, set `--localcores` to the observed guest-visible CPU count, capped at 6.
5. Set `--localmem` to the integer GiB represented by observed guest `MemTotal`, minus 2 GiB for the guest OS. Record the calculation. Do not claim that this meets the official 64 GB recommendation.
6. Set user open-file limit to at least `16384` and user process limit to at least `384` when the guest permits it. Record the observed limits used for the run.

## Actual disk allocation policy

Do not require 500 GB of free space. Do not preallocate a 500 GB image.

Use measured storage as follows:

1. Keep the verified FASTQs on `I:`. Do not create another Windows-host copy.
2. Use `F:` as the only Windows volume for new files created by this task. Do not select `E:` based on free space.
3. Store the Cell Ranger archive, extracted Cell Ranger directory, reference archive, extracted reference, QEMU data disk, Cell Ranger pipestance, temporary files, and final outputs on `F:`.
4. Use a sparse, dynamically growing QEMU data disk on `F:`. Set its virtual capacity from current `F:` free space while leaving at least `20000000000` bytes free on `F:`. Record the formula, `F:` free bytes, virtual capacity, and initial physical image bytes.
5. If the FASTQs cannot be exposed read-only to the guest with the existing QEMU setup, copy them once into the guest data disk and record the exact copied bytes and SHA256 verification. Do not retain an unnecessary second guest copy after the final result has been verified, unless it is required for Cell Ranger resume.
6. Record the actual byte size of every component after download or extraction: Cell Ranger archive, extracted Cell Ranger, reference archive, extracted reference, guest FASTQ copy if one was required, pipestance, temporary directories, and final `outs`.
7. Sample host and guest free space at least every 10 minutes during the count. Write the measurements to `disk_usage_timeseries.tsv`.
8. Stop the count cleanly if `F:` falls below `20000000000` free bytes or the guest work filesystem falls below `10000000000` free bytes. Preserve the pipestance and logs for diagnosis or resume. Do not delete user data to continue.

This policy is intended to measure the real high-water mark for this sample instead of enforcing the previous conservative free-space threshold.

## Cell Ranger command

Run from the guest work directory so the pipestance is created on the QEMU data disk:

```bash
cellranger count \
  --id=R3_SRR23490337 \
  --transcriptome=<observed_absolute_path_to_refdata-gex-GRCh38-2020-A> \
  --fastqs=<observed_guest_path_to_W6B_0_1_HGGWHDSX2> \
  --sample=bamtofastq \
  --chemistry=auto \
  --localcores=<observed_guest_cores_capped_at_6> \
  --localmem=<observed_guest_MemTotal_GiB_minus_2> \
  --disable-ui
```

Do not add `--expect-cells`, `--force-cells`, downsampling, read trimming, or unverified flags.

Capture the complete invocation, start time, end time, exit code, stdout, stderr, Martian logs, peak guest memory, peak QEMU host memory, and disk high-water mark.

Monitor chemistry detection. If Cell Ranger reports a chemistry other than `SC3Pv2`, stop the run, preserve logs, and report the exact observed chemistry. Do not override chemistry in the same task.

Do not delete or restart the pipestance after a failed count. This request authorizes one count attempt only.

## Required lightweight upload

Upload lightweight outputs to:

```text
bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260714_prjna932556_cellranger72_single_sample_hardrun
```

Required files:

- `STATUS.md`
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

Do not upload full FASTQs, reference archives, Cell Ranger archives, BAM files, molecule information files, raw/filtered matrix directories, HDF5 matrices, QEMU images, or the complete pipestance to GitHub. Keep those PC-local and record exact paths and checksums in the lightweight audit.

## Status definitions

Use exactly one of these final status values:

- `CELLRANGER72_COUNT_SUCCESS`
- `CELLRANGER72_COUNT_FAILED_MEMORY`
- `CELLRANGER72_COUNT_FAILED_DISK`
- `CELLRANGER72_COUNT_FAILED_CHEMISTRY`
- `CELLRANGER72_COUNT_FAILED_OTHER`
- `BLOCKED_CELLRANGER72_MANUAL_DOWNLOAD`
- `BLOCKED_REFERENCE_DOWNLOAD`
- `BLOCKED_QEMU_SETUP`

Success requires exit code 0, an observed `outs` directory, `metrics_summary.csv`, `web_summary.html`, and a complete output inventory. A preflight pass alone is not success.

## Stop boundaries

- Do not process another `PRJNA932556` run.
- Do not delete existing PC data.
- Do not change the verified FASTQ files.
- Do not substitute software or reference versions.
- Do not perform Seurat, CellChat, LIANA, NicheNet, differential expression, or cell annotation in this task.
- Do not push large binary analysis outputs to GitHub.
