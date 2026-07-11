# 20260711 GSE235863 timepoint rebuild and PRJNA932556 pilot

## Status

- GSE235863 CD45: completed timepoint-preserving sample and patient-level composition rebuild.
- GSE235863 CD8: completed timepoint-preserving doublet sensitivity tables after excluding unresolved tissue/timepoint samples.
- PRJNA932556: stopped after pilot with `blocked_read_structure_not_validated`; no full count matrices were generated.

## Key audit decisions

- `GSM7510911` raw GEO fields were preserved. Its title indicates `P27`, while `characteristics_ch1::patient` is `P18`; `P27-pre-P` was excluded from analytical mapping.
- Pre/post and blood/liver tumor strata were not collapsed.
- SRA metadata alone was not used to infer 10x chemistry or force single-cell quantification.

## Output directory

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260711_gse235863_timepoint_prjna932556`
