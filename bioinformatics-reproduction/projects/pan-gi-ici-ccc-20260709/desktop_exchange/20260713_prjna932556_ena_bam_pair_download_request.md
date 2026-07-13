# PC request: download and verify SRR23490337 ENA BAM/BAI pair

Date: 2026-07-13

Project: pan-GI ICI cell communication atlas

## Purpose

The previous ENA acquisition task stopped correctly because the request required exactly one submitted file, while ENA returned a matched `BAM;BAI` pair for `SRR23490337`.

This task replaces that overly strict gate. Download and verify the exact ENA-submitted BAM and BAI files for `SRR23490337`. Do not run `bamtofastq`, Cell Ranger, Seurat, or any quantification in this task.

## Required ENA query

Query the ENA Portal API read-run file report for exactly `SRR23490337`:

```text
https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR23490337&result=read_run&fields=run_accession,sample_accession,experiment_accession,study_accession,instrument_platform,instrument_model,library_layout,library_strategy,library_source,library_selection,read_count,base_count,submitted_format,submitted_ftp,submitted_md5,submitted_bytes,fastq_ftp,fastq_md5,fastq_bytes,sra_ftp,sra_md5,sra_bytes&format=tsv&download=false
```

Save the raw response exactly as `ena_srr23490337_file_report.tsv`.

## Download gate

Proceed only if all checks below pass exactly:

- `run_accession` is `SRR23490337`
- `study_accession` is `PRJNA932556`
- `sample_accession` is `SAMN33196641`
- `experiment_accession` is `SRX19384059`
- `submitted_format` is exactly `BAM;BAI`
- `submitted_ftp` contains exactly these two semicolon-separated entries in this order:
  - `ftp.sra.ebi.ac.uk/vol1/run/SRR234/SRR23490337/R3_possorted_genome_bam.bam.1`
  - `ftp.sra.ebi.ac.uk/vol1/run/SRR234/SRR23490337/R3_possorted_genome_bam.bam.1.bai`
- `submitted_md5` contains exactly these two semicolon-separated entries in this order:
  - `056a5c996212a64af708d48c5b585e48`
  - `f391d287338dbf9d3856630338957bb8`
- `submitted_bytes` contains exactly these two semicolon-separated entries in this order:
  - `21494429186`
  - `10390008`

If any check fails, stop and upload a failure audit. Do not choose a different file.

## Download target

Download both files to PC-local storage under:

```text
I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260713_prjna932556_ena_bam_pair_download
```

Preserve ENA filenames exactly:

- `R3_possorted_genome_bam.bam.1`
- `R3_possorted_genome_bam.bam.1.bai`

Before download, check free disk space on `I:`. Required free space is at least:

```text
(21494429186 + 10390008) * 1.20 = 25857455032.8 bytes
```

Record observed free bytes before download.

## Verification

For each downloaded file, record:

- ENA FTP URL
- expected bytes
- observed local bytes
- byte match result
- expected MD5
- observed MD5
- MD5 match result
- observed SHA256
- absolute PC local path

Success requires both files to match expected bytes and expected MD5.

## Runtime audit

Repeat the runtime audit, but force readable text output if possible:

```powershell
chcp 65001
wsl.exe --status
wsl.exe -l -v
docker.exe --version
```

This task does not require installing WSL, Docker, Cell Ranger, or `bamtofastq`.

## Required upload directory

Upload lightweight results to:

```text
bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260713_prjna932556_ena_bam_pair_download
```

Do not upload the BAM or BAI files to GitHub.

Required uploaded files:

- `STATUS.md`
- `ena_api_request_log.txt`
- `ena_srr23490337_file_report.tsv`
- `ena_bam_pair_download_audit.tsv`
- `ena_bam_pair_file_inventory.tsv`
- `linux_runtime_audit.txt`
- `SHA256SUMS.txt`

## Success definition

Mark success only if:

- the exact ENA row was retrieved,
- both exact submitted files were downloaded,
- BAM observed bytes equal `21494429186`,
- BAM observed MD5 equals `056a5c996212a64af708d48c5b585e48`,
- BAI observed bytes equal `10390008`,
- BAI observed MD5 equals `f391d287338dbf9d3856630338957bb8`,
- the BAM and BAI remain PC-local,
- all required lightweight audit files are uploaded.

If successful, stop after upload. The next task will be a separate `bamtofastq` feasibility pilot.
