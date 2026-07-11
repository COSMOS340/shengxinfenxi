# PC Task: PRJNA932556 ENA BAM Acquisition and Integrity Audit

Created: 2026-07-11

## Goal

Resolve the failed SRA streaming pilot by using the archive route documented for this dataset: query the European Nucleotide Archive (ENA), download the single submitted BAM associated with `SRR23490337`, and verify its integrity.

This task covers one run only:

- run accession: `SRR23490337`
- expected study: `PRJNA932556`
- expected sample label from the existing project mapping: `R3`

Do not process the other five runs. Do not run `bamtofastq`, Cell Ranger, Seurat, or any quantification step in this task.

## Why this route is required

The Windows SRA Toolkit pilot did not observe any read. The exact requested command contained an unsupported option, and the compatibility command attempted remote partial extraction for more than 30 minutes without producing a file.

Published analyses of `PRJNA932556` report obtaining the raw data from ENA as BAM files and using 10x Genomics `bamtofastq` before Cell Ranger counting. Therefore the next reproducible step is to retrieve and verify the ENA-submitted BAM, not to repeat remote `fastq-dump` streaming.

## Phase 1: exact ENA file report

Query the ENA Portal API for `SRR23490337` with `result=read_run` and save the unmodified TSV response.

Required fields:

```text
run_accession,sample_accession,experiment_accession,study_accession,instrument_platform,instrument_model,library_layout,library_strategy,library_source,library_selection,read_count,base_count,submitted_format,submitted_ftp,submitted_md5,submitted_bytes,fastq_ftp,fastq_md5,fastq_bytes,sra_ftp,sra_md5,sra_bytes
```

Record:

- exact API URL
- retrieval timestamp and timezone
- HTTP status
- exit code
- stdout and stderr
- the complete returned row without manually changing identifiers or file names

Use the returned fields as the authority. Do not infer a file name, extension, path, checksum, size, sample accession, or experiment accession.

## Phase 2: download gate

Continue only when all conditions below are satisfied by the returned ENA row:

- `run_accession` is exactly `SRR23490337`
- `study_accession` is exactly `PRJNA932556`
- `submitted_ftp` contains exactly one submitted file entry
- `submitted_md5` contains exactly one checksum entry
- `submitted_bytes` contains exactly one byte-size entry
- ENA reports the submitted file format as BAM

If any condition fails, stop and report the exact returned values. Do not select another file or run.

Before download, record available bytes on the target volume. Require available space greater than the reported BAM size plus 20%. If this check fails, stop without deleting existing project data.

## Phase 3: download and checksum

Download only the exact ENA-submitted BAM reported for `SRR23490337` to a PC-local directory on the `I:` drive. Use an ENA-provided download link and preserve the reported file name.

Requirements:

- use a resumable download method
- preserve the incomplete file if network interruption occurs
- record the exact download command or application settings
- record start time, end time, exit code, final byte size, and absolute path
- compute local MD5 after completion
- require the local MD5 to equal `submitted_md5`
- compute local SHA256 for project inventory

Do not upload the BAM to GitHub.

## Phase 4: environment audit for the later Linux step

Run and record:

```powershell
wsl.exe --status
wsl.exe -l -v
docker.exe --version
```

An unavailable command must be recorded exactly as unavailable. Do not install WSL, Docker, Cell Ranger, or `bamtofastq` in this task.

## Required outputs

- `STATUS.md`
- `ena_api_request_log.txt`
- `ena_srr23490337_file_report.tsv`
- `ena_bam_download_audit.tsv`
- `ena_bam_file_inventory.tsv`
- `linux_runtime_audit.txt`
- `SHA256SUMS.txt`

The BAM inventory must contain:

- run accession
- ENA sample accession
- ENA experiment accession
- study accession
- reported submitted file name
- reported submitted bytes
- observed local bytes
- reported MD5
- observed MD5
- MD5 match result
- observed SHA256
- absolute PC-local path
- download status

## Output directory

Upload only the lightweight audit files to:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260711_prjna932556_ena_bam_acquisition`

## Success condition

Success requires one complete PC-local BAM for `SRR23490337` whose byte size and MD5 match the exact ENA report.

The task is not successful if only metadata were retrieved, if the download is incomplete, or if the checksum does not match.

## Completion signal

Commit and push to:

`pan-gi-ici-metadata-request-20260709`

Commit message:

`Acquire and verify SRR23490337 ENA BAM`
