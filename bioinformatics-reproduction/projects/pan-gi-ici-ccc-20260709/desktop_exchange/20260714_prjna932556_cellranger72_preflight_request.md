# PC request: PRJNA932556 Cell Ranger 7.2.0 resource preflight

Date: 2026-07-14

## Purpose

Determine whether the PC can safely run one complete `cellranger count` for `SRR23490337` / sample `R3` using the already verified full FASTQs.

This request is an audit only. Do not download large archives, install Cell Ranger, extract the reference, or start `cellranger count`.

## Verified input state

Use this exact PC-local FASTQ root:

`I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260714_prjna932556_srr23490337_full_bamtofastq`

Expected state:

- gzip FASTQ files: `804`
- R1 files: `402`
- R2 files: `402`
- total FASTQ bytes: `29268262419`
- read pairs: `401831188`
- FASTQ sample prefix: `bamtofastq`
- R1 first-100 length: `28`
- R2 first-100 length: `91`
- existing host SHA256 mismatches: `0`

The full per-file verification table is already committed at:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260714_prjna932556_srr23490337_full_bamtofastq/host_fastq_sha256_verification.tsv`

## Exact analysis target

- sample: `R3`
- SRA run: `SRR23490337`
- response group: resistant / progression
- assay reported by the original study: Chromium Single Cell 3' Reagent Kit v2
- Cell Ranger version used by the published reanalysis: `7.2.0`
- reference used by the published reanalysis: `refdata-gex-GRCh38-2020-A`
- official reference URL: `https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz`
- official Cell Ranger download page: `https://www.10xgenomics.com/support/software/cell-ranger/downloads`

Do not substitute another Cell Ranger version or reference package.

## Required checks

### 1. Reverify the existing FASTQ directory

Record the exact directory, file count, R1/R2 count, total bytes, gzip test result, and whether every R1 chunk has one matching R2 chunk.

Do not rename, move, delete, recompress, or regenerate the FASTQs.

### 2. Audit Windows host resources

Record:

- Windows version
- CPU model
- physical core count
- logical processor count
- total physical RAM in bytes
- currently available physical RAM in bytes
- whether the CPU exposes AVX to the host
- every fixed volume with exact total and free bytes

Do not delete any data to increase free space.

### 3. Audit the Linux execution environment

Recheck WSL and Docker status. For the existing QEMU Ubuntu environment, record:

- QEMU version
- guest architecture
- allocated virtual CPU count
- allocated guest RAM in bytes
- guest-visible AVX flag status
- guest free bytes on the proposed Cell Ranger work volume
- `ulimit -n`
- `ulimit -u`

Do not start a large guest disk image or copy the FASTQs into the guest during this audit.

### 4. Audit exact software and reference availability

Check whether Cell Ranger `7.2.0` and `refdata-gex-GRCh38-2020-A` already exist under the known PC tool and project roots. Record exact paths only when observed.

If Cell Ranger is found, run only:

```bash
cellranger --version
cellranger count --help
```

Record the executable byte size and SHA256. Do not run `cellranger count`.

If Cell Ranger is absent, test access to the official download page without accepting the EULA or inventing a download URL. Record `manual_eula_required` when the exact archive URL cannot be obtained non-interactively.

For `refdata-gex-GRCh38-2020-A`, issue only an HTTP metadata request to the exact URL above. Record HTTP status, redirects, and content length when the server provides it. Do not download the archive.

### 5. Apply the resource gate

The PC is ready for a later single-sample count request only if all of the following are verified:

- at least `8` guest-visible CPU cores
- at least `64000000000` bytes of guest RAM
- AVX visible inside the Linux environment
- at least one usable fixed volume with `500000000000` free bytes before software, reference, and count outputs are created
- exact Cell Ranger `7.2.0` acquisition route is available
- exact `refdata-gex-GRCh38-2020-A` source is reachable

The `500000000000`-byte storage threshold is the project safety gate for this 401,831,188-read-pair run. It is not permission to delete existing data.

If any gate fails, report `LOCAL_CELLRANGER_COUNT_NOT_READY` and stop. Do not attempt a reduced-memory, reduced-disk, partial, or alternative quantification run.

If every gate passes, report `READY_FOR_CELLRANGER72_ACQUISITION_AND_COUNT_REQUEST` and stop. A separate request will authorize downloads and execution.

## Planned command audit

Write the following template into `planned_cellranger_count_command.txt`, replacing placeholders only with exact observed absolute paths and verified resource allocations:

```bash
cellranger count \
  --id=R3_SRR23490337 \
  --transcriptome=<absolute_path_to_refdata-gex-GRCh38-2020-A> \
  --fastqs=<absolute_path_to_directory_containing_W6B_0_1_HGGWHDSX2> \
  --sample=bamtofastq \
  --chemistry=auto \
  --localcores=<verified_guest_core_allocation> \
  --localmem=<verified_guest_memory_in_GB> \
  --disable-ui
```

This is a planning record only. Do not execute it. Chemistry auto-detection must later be reported; the count stage must stop for review if the detected chemistry is not `SC3Pv2`.

## Required lightweight upload

Upload only these audit files to:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260714_prjna932556_cellranger72_preflight`

Required files:

- `STATUS.md`
- `fastq_reverification_summary.tsv`
- `host_system_resources.tsv`
- `host_fixed_volume_audit.tsv`
- `linux_runtime_resource_audit.txt`
- `cellranger72_presence_audit.tsv`
- `official_source_access_audit.tsv`
- `planned_cellranger_count_command.txt`
- `preflight_decision.tsv`
- `upload_manifest.tsv`
- `SHA256SUMS.txt`

## Stop boundary

Do not download or install Cell Ranger. Do not download or extract the reference. Do not run Cell Ranger, STARsolo, kallisto, bustools, salmon, Seurat, CellChat, LIANA, NicheNet, or any quantification workflow. Do not process the other five PRJNA932556 runs.

