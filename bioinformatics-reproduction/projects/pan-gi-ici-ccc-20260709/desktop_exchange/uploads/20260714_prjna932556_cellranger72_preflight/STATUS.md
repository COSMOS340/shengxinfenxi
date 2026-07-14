# Cell Ranger 7.2.0 preflight status

Status: LOCAL_CELLRANGER_COUNT_NOT_READY
Completed: 2026-07-14T20:13:50+08:00

## Scope

- Request: `20260714_prjna932556_cellranger72_preflight_request.md`
- Sample: `R3`
- Run: `SRR23490337`
- Audit only: no downloads, no installation, no Cell Ranger count, no alternative quantification
- Full FASTQ upload to GitHub: no

## FASTQ Reverification

- FASTQ root: `I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\desktop_exchange\local_only\20260714_prjna932556_srr23490337_full_bamtofastq`
- Gzip FASTQ files: `804`
- R1 files: `402`
- R2 files: `402`
- Total FASTQ bytes: `29268262419`
- Pairing passed: `True`
- Existing host SHA256 mismatches: `0`

## Resource Gate

- Guest-visible CPU cores: `4` / required `8`
- Guest RAM bytes: `4104351744` / required `64000000000`
- Guest AVX visible: `True`
- Maximum fixed-volume free bytes: `417535246336` / required `500000000000`
- Cell Ranger 7.2.0 local executable: `False`
- Cell Ranger archive route: manual EULA flow required
- GRCh38-2020-A source reachable: `True`

## Decision

`LOCAL_CELLRANGER_COUNT_NOT_READY`

Failed criteria: `guest_visible_cpu_cores, guest_ram_bytes, fixed_volume_free_bytes, cellranger_7_2_0_acquisition_route`
