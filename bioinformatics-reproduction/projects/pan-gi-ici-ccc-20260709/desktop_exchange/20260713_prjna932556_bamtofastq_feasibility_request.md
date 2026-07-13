# PC request: SRR23490337 bamtofastq feasibility pilot

Date: 2026-07-13

Project: pan-GI ICI cell communication atlas

## Purpose

The ENA-submitted `SRR23490337` BAM/BAI pair has been downloaded and verified. This task tests whether the verified BAM can be converted back to FASTQ with official 10x `bamtofastq`.

Do not run Cell Ranger, Seurat, CellChat, LIANA, NicheNet, or quantification in this task.

Do not run a full BAM-to-FASTQ conversion unless the bounded pilot described below succeeds and this request explicitly allows the bounded output. Keep all FASTQ outputs PC-local.

## Verified input files

Use the existing PC-local files:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260713_prjna932556_ena_bam_pair_download\R3_possorted_genome_bam.bam.1
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260713_prjna932556_ena_bam_pair_download\R3_possorted_genome_bam.bam.1.bai
```

Before using them, re-check:

- BAM bytes exactly `21494429186`
- BAM MD5 exactly `056a5c996212a64af708d48c5b585e48`
- BAI bytes exactly `10390008`
- BAI MD5 exactly `f391d287338dbf9d3856630338957bb8`

If any check fails, stop and upload failure logs.

## Runtime and tool audit

Create readable logs for:

```powershell
chcp 65001
wsl.exe --status
wsl.exe -l -v
docker.exe --version
```

Inside the Linux environment used for conversion, record:

```bash
uname -a
pwd
which bash
which samtools || true
samtools --version || true
which bamtofastq || true
bamtofastq --version || true
bamtofastq --help || true
```

If `bamtofastq` is absent, obtain the official 10x Genomics `bamtofastq` Linux binary only from an official 10x Genomics page or the official `10XGenomics/bamtofastq` GitHub releases page. Record the exact download URL, version string, local path, file size, and SHA256.

If the official binary cannot be obtained reproducibly, stop and upload the audit. Do not use a third-party converter.

## BAM header and contig audit

If `samtools` is available, run:

```bash
samtools quickcheck -v R3_possorted_genome_bam.bam.1
samtools view -H R3_possorted_genome_bam.bam.1 > bam_header.sam
samtools idxstats R3_possorted_genome_bam.bam.1 > bam_idxstats.tsv
```

If `samtools` is not available, record that explicitly and continue only if `bamtofastq --help` and `bamtofastq --version` work.

For the bounded pilot locus, use this exact rule:

- If `bam_idxstats.tsv` contains a contig name exactly equal to `chrM`, use `chrM:1-16569`.
- Else if it contains a contig name exactly equal to `MT`, use `MT:1-16569`.
- Else stop before conversion and upload the audit.

Do not infer any other contig name.

## Bounded bamtofastq pilot

Run a locus-restricted pilot only:

```bash
bamtofastq --locus=<validated_locus> --reads-per-fastq=1000000 R3_possorted_genome_bam.bam.1 <output_dir>
```

Use output directory:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260713_prjna932556_bamtofastq_feasibility
```

Record start time, end time, exit code, stdout, stderr, and peak disk footprint if available.

After the pilot, create a FASTQ inventory table with:

- relative FASTQ path
- file size
- gzip status
- first read header
- number of records inspected
- read length distribution from the first 100 records where possible

Upload at most a small tarball containing the first 100 records per FASTQ file. Do not upload complete FASTQ files.

## Required upload directory

Upload lightweight outputs to:

```text
bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260713_prjna932556_bamtofastq_feasibility
```

Required files:

- `STATUS.md`
- `input_bam_bai_reverification.tsv`
- `runtime_tool_audit.txt`
- `bam_header.sam` if generated
- `bam_idxstats.tsv` if generated
- `bamtofastq_command_log.txt`
- `bamtofastq_fastq_inventory.tsv`
- `bamtofastq_first100_records.tar.gz` if FASTQ files were generated
- `SHA256SUMS.txt`

## Success definition

Mark success only if:

- the BAM/BAI pair was reverified,
- an official 10x `bamtofastq` binary was used,
- the bounded locus was selected by the exact `chrM` or `MT` rule,
- `bamtofastq` exited with code `0`,
- FASTQ files were produced in PC-local storage,
- first-record/read-length audit was uploaded,
- full FASTQ files were not uploaded to GitHub.

If the task stops earlier, mark the precise stop reason in `STATUS.md`.
