# SRR23490337 full bamtofastq status

Result: STOPPED_BEFORE_CONVERSION
Checked: 2026-07-14T13:21:24.1008016+08:00

## Verification

- BAM bytes and MD5: PASS
- BAI bytes and MD5: PASS
- Official 10x `bamtofastq` v1.4.1 SHA256: PASS
- Fresh output path requirement: PASS; the requested output directory did not exist
- Disk gate: FAIL

## Precise stop reason

- Output volume: `I:`
- Observed free bytes: `36503392256`
- Required free bytes: `60000000000`
- Shortfall bytes: `23496607744`

The mandatory 60 GB disk gate failed. The full conversion command was not executed, the portable Linux runtime was not started, the full output directory was not created, and no existing data was deleted.

The planned command remains:

```bash
bamtofastq --reads-per-fastq=1000000 R3_possorted_genome_bam.bam.1 <output_dir>
```

No FASTQ files were generated or uploaded. Cell Ranger, Seurat, CellChat, LIANA, NicheNet, and quantification were not run.
