# SRR23490337 bamtofastq feasibility status

Result: SUCCESS
Completed: 2026-07-13T16:54:00.569334+00:00

## Input verification

- Host BAM bytes and MD5: PASS
- Host BAI bytes and MD5: PASS
- Guest BAM MD5: PASS
- Guest BAI MD5: PASS
- `samtools quickcheck -v`: exit 0 with no reported errors

## Runtime

- Windows WSL feature: disabled
- Docker: not found
- Portable runtime used instead: QEMU 11.0.50 with Ubuntu 24.04 LTS, stored under `I:\`
- `samtools`: 1.19.2
- Official 10x `bamtofastq`: v1.4.1
- Official binary SHA256: `fdf7db4fe6cf8e13a432e8a59815dad506415349627502fdf77f4be208a198ce`

## Bounded pilot

- Exact contig rule selected: `MT:1-16569`
- Command: `bamtofastq --locus=MT:1-16569 --reads-per-fastq=1000000 R3_possorted_genome_bam.bam.1 <output_dir>`
- Final exit code: `0`
- Read pairs observed and written: `76559661`
- FASTQ files: `154`
- FASTQ bytes: `5386519028`
- Gzip files: `154/154`
- First records audited: `15400` total, 100 per FASTQ
- R1 read length in first 100 records: `28`
- R2 read length in first 100 records: `91`
- Peak output footprint: `5386519028` bytes
- Full FASTQ local path: `I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260713_prjna932556_bamtofastq_feasibility`

The first attempt stopped before writing FASTQ because the empty output directory already existed; official `bamtofastq` requires a non-existing output path. The attempt-1 logs are retained. The corrected final attempt exited 0.

Full FASTQ files remain PC-local and are not included in the GitHub upload. No Cell Ranger, Seurat, CellChat, LIANA, NicheNet, or quantification was run.
