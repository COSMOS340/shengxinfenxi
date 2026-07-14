# SRR23490337 full bamtofastq status

Status: SUCCESS
Completed: 2026-07-14T19:28:04+08:00

## Request

- Request file: `20260714_prjna932556_srr23490337_full_bamtofastq_request.md`
- Manifest file: `20260714_prjna932556_srr23490337_full_bamtofastq_manifest.tsv`
- Upload directory: `bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260714_prjna932556_srr23490337_full_bamtofastq`
- Full FASTQ location: `I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260714_prjna932556_srr23490337_full_bamtofastq`
- Full FASTQ upload to GitHub: no

## Verified Inputs

- BAM bytes: `21494429186`
- BAM MD5: `056a5c996212a64af708d48c5b585e48`
- BAI bytes: `10390008`
- BAI MD5: `f391d287338dbf9d3856630338957bb8`

## Command

```text
/work/tools/bamtofastq --reads-per-fastq=1000000 /work/input/R3_possorted_genome_bam.bam.1 /work/output
```

- `--locus`: not used
- Exit code: `0`
- Read pairs written: `401831188`
- Wall time: `3:07:01`

## Output Audit

- Gzip FASTQ files: `804`
- R1 files: `402`
- R2 files: `402`
- I1 files: `0`
- I2 files: `0`
- Total FASTQ bytes: `29268262419`
- Pairing audit passed: `True`
- Host copy verification passed: `True`
- Host SHA256 mismatches: `0`

## Uploaded Lightweight Files

The upload directory contains command logs, input re-verification, runtime/tool audit, disk audit, full FASTQ inventory, per-file SHA256 table, pairing audit, read-length audit, first-100-record tarball, and host-copy verification tables. Complete FASTQ files remain PC-local only.
