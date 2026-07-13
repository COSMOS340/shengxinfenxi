# 20260713 PRJNA932556 ENA BAM/BAI Pair Download

## Status

Not successful. The ENA row gate passed, and the BAI file completed and matched MD5, but the BAM download did not complete.

## Gate

- run_accession: `SRR23490337`
- sample_accession: `SAMN33196641`
- experiment_accession: `SRX19384059`
- study_accession: `PRJNA932556`
- submitted_format: `BAM;BAI`
- gate_status: `passed`

## Download

- BAM expected bytes: `21494429186`
- BAM observed transferred bytes before suspension: `13917920681`
- BAM MD5: `not_computed_incomplete`
- BAI expected bytes: `10390008`
- BAI observed bytes: `10390008`
- BAI MD5 match: `True`

## Decision

`incomplete_bam_download_not_successful`

The BAM/BAI success definition is not met. The BITS BAM transfer was suspended after repeated remote connection closures.
