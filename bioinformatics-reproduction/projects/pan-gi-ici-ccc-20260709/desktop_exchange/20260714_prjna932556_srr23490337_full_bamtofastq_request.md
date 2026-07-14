# PC request: full bamtofastq for SRR23490337 only

Date: 2026-07-14

Project: pan-GI ICI cell communication atlas

## Purpose

The bounded `bamtofastq` pilot for `SRR23490337` succeeded. This task performs full `bamtofastq` conversion for this one run only.

Do not process the other PRJNA932556 runs in this task.

Do not run Cell Ranger, Seurat, CellChat, LIANA, NicheNet, or quantification in this task.

Keep all full FASTQ files PC-local. Upload only lightweight audits and first-record summaries.

## Verified input files

Use the existing PC-local files:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260713_prjna932556_ena_bam_pair_download\R3_possorted_genome_bam.bam.1
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260713_prjna932556_ena_bam_pair_download\R3_possorted_genome_bam.bam.1.bai
```

Before conversion, re-check:

- BAM bytes exactly `21494429186`
- BAM MD5 exactly `056a5c996212a64af708d48c5b585e48`
- BAI bytes exactly `10390008`
- BAI MD5 exactly `f391d287338dbf9d3856630338957bb8`

If any check fails, stop and upload failure logs.

## Tool requirements

Use the same official 10x `bamtofastq` binary validated in the feasibility task:

- path on PC: `I:\codex-config\tools\bamtofastq\v1.4.1\bamtofastq_linux`
- official URL: `https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux`
- version: `v1.4.1`
- SHA256: `fdf7db4fe6cf8e13a432e8a59815dad506415349627502fdf77f4be208a198ce`

Record the runtime and tool audit again. If this exact binary cannot be used, stop and upload the audit.

## Disk gate

The bounded `MT:1-16569` pilot wrote 5,386,519,028 bytes from 76,559,661 read pairs. `bam_idxstats.tsv` reported 390,129,881 mapped records across contigs.

Before full conversion, record free bytes on the output volume. Proceed only if available free space is at least `60000000000` bytes.

This 60 GB gate is intentionally conservative for this single-run full conversion.

If free space is below this gate, stop and upload a disk audit. Do not delete existing data unless explicitly instructed by the user.

## Full conversion command

Use a fresh output directory that does not already exist:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260714_prjna932556_srr23490337_full_bamtofastq
```

Run full conversion without `--locus`:

```bash
bamtofastq --reads-per-fastq=1000000 R3_possorted_genome_bam.bam.1 <output_dir>
```

Record:

- start time
- end time
- exit code
- stdout
- stderr
- `/usr/bin/time -v` output if available
- output disk footprint over time if available

## Output audit

After conversion, create:

- `full_bamtofastq_fastq_inventory.tsv`
- `full_bamtofastq_local_fastq_sha256.tsv`
- `full_bamtofastq_pairing_audit.tsv`
- `full_bamtofastq_readlength_audit.tsv`

The inventory must include, for every FASTQ file:

- relative FASTQ path
- file size
- gzip status
- first read header
- number of records inspected
- read length distribution from the first 100 records where possible

The pairing audit must report:

- number of R1 files
- number of R2 files
- number of I1/I2 files if any
- whether every R1 chunk has a matching R2 chunk by lane and chunk number
- total gzip FASTQ count
- total FASTQ bytes

Upload at most a small tarball containing the first 100 records per FASTQ file. Do not upload complete FASTQ files.

## Required upload directory

Upload lightweight outputs to:

```text
bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260714_prjna932556_srr23490337_full_bamtofastq
```

Required files:

- `STATUS.md`
- `input_bam_bai_reverification.tsv`
- `runtime_tool_audit.txt`
- `full_bamtofastq_command_log.txt`
- `full_bamtofastq_disk_audit.tsv`
- `full_bamtofastq_fastq_inventory.tsv`
- `full_bamtofastq_local_fastq_sha256.tsv`
- `full_bamtofastq_pairing_audit.tsv`
- `full_bamtofastq_readlength_audit.tsv`
- `full_bamtofastq_first100_records.tar.gz` if FASTQ files were generated
- `SHA256SUMS.txt`

## Success definition

Mark success only if:

- the BAM/BAI pair was reverified,
- the validated official 10x `bamtofastq` v1.4.1 binary was used,
- disk gate passed before conversion,
- full `bamtofastq` exited with code `0`,
- gzip FASTQ files were produced in PC-local storage,
- R1/R2 pairing audit passed,
- read length audit was uploaded,
- full FASTQ files were not uploaded to GitHub.

If the task stops earlier, mark the precise stop reason in `STATUS.md`.
