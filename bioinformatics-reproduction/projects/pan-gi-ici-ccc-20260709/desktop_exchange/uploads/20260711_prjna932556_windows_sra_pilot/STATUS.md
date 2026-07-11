# 20260711 PRJNA932556 Windows SRA Toolkit Pilot

## Status

- SRA Toolkit 3.4.1 for Windows was installed under `I:/codex-config/tools/sratoolkit/sratoolkit.3.4.1-win64`.
- The downloaded archive MD5 matched the official checksum recorded in the request audit.
- The exact requested `fastq-dump` command failed because SRA Toolkit 3.4.1 does not support `--include-technical`.
- A documented compatibility command without `--include-technical` was attempted; it was terminated after prolonged remote partial extraction with no FASTQ, stdout, or stderr output.
- No complete `SRR23490337` run was downloaded, no Cell Ranger step was run, and no count matrix was generated.

## Decision

`remote_partial_extraction_failed_no_reads_observed`

Barcode, UMI, cDNA, and chemistry roles remain `not_validated` because no FASTQ records were observed.
