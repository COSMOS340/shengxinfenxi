# PC Task: PRJNA932556 Windows SRA Toolkit Read-Structure Pilot

Created: 2026-07-11

## Goal

Install the official NCBI SRA Toolkit for Windows and inspect a partial public run from `PRJNA932556`.

This task is a read-structure pilot only. Do not download all six runs, do not run Cell Ranger, and do not generate a count matrix.

## Why this request replaces the previous pilot

The previous PC run only checked PATH. No SRA tool was installed, no FASTQ file was generated, and no read was observed. Therefore `blocked_read_structure_not_validated` describes an environment block rather than evidence that the deposited reads are unusable.

## Required tool installation

Download the current official MS Windows 64-bit SRA Toolkit release from the NCBI SRA Toolkit release page.

Requirements:

- install or extract in a user-writable directory
- do not require administrator access
- record download URL, archive name, archive checksum if supplied, install path, and exact tool version
- add the toolkit `bin` directory to PATH for this task or call executables by absolute path

Required executable checks:

- `fastq-dump --version`
- `prefetch --version`
- `vdb-validate --version`

## Pilot run

Use only:

- run: `SRR23490337`
- sample: `R3`
- requested spots: first 100,000

Use `fastq-dump`, not `fasterq-dump`, because the pilot requires a bounded spot range.

Required extraction behavior:

- split reads into separate files
- include technical reads
- preserve read identifiers
- do not clip reads
- stop after 100,000 spots

Run the equivalent of:

```powershell
fastq-dump.exe --split-files --include-technical --readids --dumpbase -X 100000 --outdir <pilot_output_directory> SRR23490337
```

Record the exact command actually executed. Do not silently replace options.

If remote partial extraction fails, report the exact command, exit code, stdout, and stderr. Do not fall back to downloading the complete 8,033 MB run in this task.

## Required audit

For every FASTQ file generated, calculate:

- file name
- file size
- read number or read role
- number of records
- minimum read length
- median read length
- maximum read length
- first 20 headers
- whether all files contain the same number of spots
- whether read names preserve a common spot identifier across files

Create a gzip file containing the first 100 FASTQ records from each generated read file. These are public sequencing records and are included only for read-structure review.

Do not assign barcode, UMI, cDNA, or chemistry roles unless supported by the observed number of read files, read lengths, sequences, headers, and the original paper method. Record unsupported roles as `not_validated`.

## Required outputs

- `STATUS.md`
- `sra_toolkit_install_audit.tsv`
- `prjna932556_pilot_command_log.txt`
- `prjna932556_pilot_file_inventory.tsv`
- `prjna932556_pilot_read_metrics.tsv`
- `prjna932556_pilot_header_audit.tsv`
- `prjna932556_pilot_role_assessment.tsv`
- `prjna932556_pilot_first100_records.tar.gz`
- `SHA256SUMS.txt`

## Output directory

Upload lightweight outputs to:

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260711_prjna932556_windows_sra_pilot`

Do not upload the full 100,000-record FASTQ files if the combined archive is large. Keep them PC-local and list their absolute paths and SHA256 values in the file inventory.

## Stop condition

Stop after the partial read audit. Full processing requires a separate decision because Cell Ranger is Linux-only and the complete SRA-to-FASTQ conversion may require hundreds of gigabytes of working space.

## Completion signal

Commit and push to:

`pan-gi-ici-metadata-request-20260709`

Commit message:

`Add PRJNA932556 Windows SRA read-structure pilot`
