# GSE189926 R-only rerun status

Status: complete
Language: R

- Raw input: 22 gzip text matrices
- Cells imported: 121452
- QC-pass cells: 89587
- Feature union: 51757
- Unintegrated routes: unfiltered and QC-pass Seurat log-normalized PCA/UMAP
- Sample-aware route: Harmony by sample_accession
- Clustering resolutions audited: 0.3, 0.5, 0.8
- Annotation resolution: 0.5
- Marker export: raw and technical-gene-excluded top 50 per cluster
- Large RDS objects: PC-local only and recorded in object_inventory.tsv
- GSE235863 CD45/CD8 summaries: rebuilt in R from existing uploaded count tables
- Downsampling: none

Stop conditions respected: CellChat, LIANA, NicheNet, differential communication, final statistical testing, and manuscript figures were not run.
