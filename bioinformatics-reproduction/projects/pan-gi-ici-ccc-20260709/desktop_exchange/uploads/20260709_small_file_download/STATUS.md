# Pan-GI ICI CCC small file download

Completed: 2026-07-11 refresh of manual rows

GitHub handoff directory: bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_small_file_download
Local input directory: I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\02_input

Summary:
- Manifest rows: 14
- Direct URL downloads: 9 downloaded, 0 failed
- Manual rows: 1 manual_done, 4 verified_unavailable_public_sources; no unresolved manual rows remain
- No RDS, h5ad, FASTQ, BAM/CRAM, or spatial image files were downloaded.

Manual resolution notes:
- STAD_PRJEB25780_CIDE_metadata: verified_unavailable_public_sources; CIDE pages are reachable, but no reproducible STAD-PRJEB25780 sample metadata export endpoint/static table was found after refresh.
- GSE205506_response_mapping: verified_unavailable_public_sources; GEO sample metadata retained, but exact pCR/non-pCR response labels were not publicly accessible from GEO/PubMed/web/Cell supplement paths checked on this PC.
- GSE236581_response_mapping: verified_unavailable_public_sources; GEO sample metadata retained, but exact CR/PR/SD labels were not publicly accessible from GEO/SRA/PubMed/web/Cell supplement paths checked on this PC.
- PRJNA932556_sample_response_mapping: manual_done; SRA runinfo maps sample names to accessions, and article-derived snippets identify S1-S3 as sensitive/disease remission and R1-R3 as resistant/progression.
- PRJNA932556_processed_matrix_links: verified_unavailable_public_sources; public metadata/search paths expose sequencing runs but no processed matrix endpoint.

Evidence is recorded in manual_source_audit.tsv. Unknown response labels were not inferred.
