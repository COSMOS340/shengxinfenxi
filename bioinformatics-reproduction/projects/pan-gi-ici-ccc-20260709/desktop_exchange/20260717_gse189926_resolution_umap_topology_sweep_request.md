# PC request: GSE189926 resolution and UMAP topology sweep

Status: `REQUESTED_GSE189926_RESOLUTION_UMAP_SWEEP`

Date: 2026-07-17

## Purpose

The current GSE189926 Harmony UMAP shows long, thin, filament-like structures and weak visual separation in several regions. Determine whether this is driven by cluster resolution, UMAP parameters, residual QC structure, sample dominance, or the current Harmony representation.

Cluster resolution does not change UMAP coordinates. This request therefore contains two separate controlled experiments:

1. a cluster-resolution sweep on one fixed Harmony graph and embedding;
2. a UMAP-parameter sweep on one fixed Harmony reduction.

Do not select a result only because it looks more compact. The accepted result must preserve marker-supported biology, avoid sample-dominated islands, and reduce unsupported filament-like geometry.

## Exact existing object

Load the existing R object recorded in the accepted inventory:

`I:/shengxinfenxi/bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/04_objects/20260710_gse189926_r_only_rerun/gse189926_r_qc_pass_annotated.rds`

Expected inventory values:

- format: Seurat RDS
- cells: 89,587
- file size: `3,739,942,614` bytes
- SHA-256: `bf112efde4b1816348c8939ab2127f00ae1e959f4132adda17a74d8673052abb`

Verify the actual path, size, SHA-256, cell count, assays, layers, reductions, graphs, and metadata field names before analysis. Stop and report the exact mismatch if any value or required object component differs. Do not infer a replacement path or field name.

Do not overwrite the existing RDS. Save all new reductions, cluster fields, and output objects under a new directory.

## Language and scope

- Use R only for object analysis, clustering, UMAP, metrics, markers, and figures.
- Do not rerun response testing, communication analysis, or downstream manuscript statistics.
- Do not change the existing QC-pass cell set in the primary sweep.
- QC metrics may be used to diagnose filament-like regions. Any alternative QC exclusion must be exported only as a separate sensitivity result and must not replace the primary object in this request.

## 1. Confirm the baseline

Confirm that the existing baseline has:

- Harmony reduction derived by `sample_accession`;
- Harmony dimensions 1:30 available;
- `harmony_snn` graph present or reproducibly rebuilt from Harmony dimensions 1:30 with `k.param = 15`;
- existing UMAP parameters `n.neighbors = 15`, `min.dist = 0.3`, `seed.use = 340`;
- existing resolution fields 0.3, 0.5, and 0.8.

Export an exact baseline audit instead of assuming these components exist.

## 2. Fixed-graph cluster-resolution sweep

Using the same `harmony_snn` graph and random seed 340, run:

```r
resolution_values <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
```

For every resolution export:

- number of clusters;
- cell count per cluster;
- clusters below 50, 100, and 200 cells;
- sample and patient contribution per cluster;
- maximum sample fraction and maximum patient fraction per cluster;
- each cluster's own Top10 and Top50 markers from RNA expression;
- broad-marker positive and conflicting evidence;
- adjusted Rand index between adjacent resolutions if the required R package is available; otherwise report that the metric was not computed and do not substitute an unverified implementation;
- cluster-to-cluster transition table between adjacent resolutions.

Draw all seven resolution labels on the same baseline UMAP coordinates in a readable side-by-side figure. This figure evaluates cluster partitioning only and must not be described as seven different UMAPs.

## 3. Fixed-reduction UMAP parameter sweep

Using the same Harmony reduction, same cell set, and seed 340, calculate these distinct UMAP reductions:

| route | Harmony dimensions | n.neighbors | min.dist |
|---|---:|---:|---:|
| baseline | 1:30 | 15 | 0.3 |
| neighbor30 | 1:30 | 30 | 0.3 |
| neighbor50 | 1:30 | 50 | 0.3 |
| compact30 | 1:30 | 30 | 0.1 |
| diffuse30 | 1:30 | 30 | 0.5 |
| diffuse50 | 1:30 | 50 | 0.5 |
| dims20_neighbor30 | 1:20 | 30 | 0.3 |

Record every actual `RunUMAP` argument, method, metric, package version, thread count, and reduction name. Do not overwrite `umap.harmony`.

For each route draw separate panels colored by:

- the exact existing `broad_label` field;
- the exact existing `cluster_harmony_r0_5` field for coordinate comparison;
- sample accession;
- patient;
- raw outcome;
- timepoint;
- `nFeature_RNA`;
- `nCount_RNA`;
- mitochondrial percentage;
- ribosomal percentage;
- hemoglobin percentage.

Use identical color scales and plotting order across routes so geometry is directly comparable.

## 4. Quantify filament-like geometry

Do not rely only on visual judgment. For each UMAP route, on a reproducible stratified sample of at most 20,000 cells:

1. calculate 30 nearest neighbors in the two-dimensional UMAP coordinates;
2. calculate the two eigenvalues of each cell's local neighbor covariance matrix;
3. define local linearity as `lambda1 / max(lambda2, .Machine$double.eps)`;
4. report median, 90th percentile, and 99th percentile local linearity overall and by resolution 0.5 cluster;
5. report the fraction of sampled cells above fixed local-linearity values 10, 20, and 50;
6. map local linearity back onto each UMAP for visual review.

This metric is a topology diagnostic, not a biological quality score. Do not remove cells based on it.

Also export for each route:

- mean and median same-sample 30-neighbor fraction;
- mean and median same-patient 30-neighbor fraction;
- number of sample-pure connected regions identified by the existing diagnostic method, if that method is already implemented;
- cluster centroid separation and within-cluster two-dimensional dispersion for the fixed resolution 0.5 labels.

## 5. Diagnose the long arms

Identify cells in the visually and quantitatively most filament-like regions without deleting them. Export:

- cell IDs;
- exact `cluster_harmony_r0_5` value;
- exact `broad_label` and `refined_label` values;
- sample and patient;
- response and timepoint;
- QC metrics;
- Top50 markers for the affected cluster;
- sample and patient contribution;
- scDblFinder or other doublet field only if that exact field already exists in the R object.

Determine whether each long arm is associated with a supported lineage trajectory, one sample or patient, a QC metric, mixed broad markers, or unsupported low-signal cells. Use evidence labels in a review table; do not rename or delete cells automatically.

## 6. Color and rendering requirements

Retain the approved bright palette from the first R/Seurat run, but increase point visibility compared with the current six-panel image.

- broad T cells: `#0072B2`
- broad CD8 T/NK: `#56B4E9`
- plasma: `#CC79A7`
- inflammatory myeloid: `#D55E00`
- macrophage: `#009E73`
- monocyte: `#006D2C`
- mast: `#7A5195`
- epithelial: `#B79F00`
- endothelial: `#17BECF`
- stromal/perivascular: `#8C564B`
- unresolved: `#666666`

Use a white background, alpha at least 0.85 for cell-type panels, and a point size that remains visible after figure reduction. Randomize plotting order with seed 340. Export both PNG and PDF. Inspect every final comparison sheet for clipping, text overlap, legend overflow, overly faint points, and colors that cannot be distinguished.

## 7. Selection rules

Do not choose cluster resolution from UMAP appearance. Rank resolutions by:

1. marker-supported broad and refined separation;
2. absence of unsupported small mixed clusters;
3. stability across adjacent resolutions;
4. sample and patient representation;
5. suitability for later compartment-specific reclustering.

Rank UMAP routes separately by:

1. reduced unsupported local linearity;
2. acceptable sample mixing without erasing patient biology;
3. preservation of known broad lineages;
4. absence of QC-dominated or one-sample arms;
5. readable publication geometry.

The selected resolution and selected UMAP route may be different decisions. Export the full evidence table and request Mac review before changing the formal GSE189926 object.

## 8. Required outputs

Upload to:

`uploads/20260717_gse189926_resolution_umap_topology_sweep`

At minimum include:

- `STATUS.md`
- `MAC_REVIEW_REQUEST.md`
- R scripts, logs, and `sessionInfo()`
- `gse189926_sweep_input_object_audit.tsv`
- `gse189926_resolution_sweep_summary.tsv`
- `gse189926_resolution_cluster_counts.tsv`
- `gse189926_resolution_sample_patient_dominance.tsv`
- `gse189926_resolution_transition_tables.tsv`
- `gse189926_resolution_top10_markers.tsv`
- `gse189926_resolution_top50_markers.tsv`
- `gse189926_umap_route_parameters.tsv`
- `gse189926_umap_route_mixing_metrics.tsv`
- `gse189926_umap_local_linearity_summary.tsv`
- `gse189926_umap_local_linearity_by_cluster.tsv`
- `gse189926_filament_region_cell_audit.tsv.gz`
- `gse189926_filament_region_review.tsv`
- `gse189926_resolution_selection_evidence.tsv`
- `gse189926_umap_selection_evidence.tsv`
- fixed-coordinate seven-resolution comparison as PNG and PDF
- broad-label seven-route UMAP comparison as PNG and PDF
- fixed-cluster seven-route UMAP comparison as PNG and PDF
- sample, patient, QC-metric, and local-linearity comparison sheets as PNG and PDF
- `output_manifest.tsv` with relative paths, sizes, SHA-256 values, and local-only status

Save the sweep object locally without overwriting the prior object, and list its exact path, size, and SHA-256 in the manifest.

## Completion gate

The request is complete only when the seven resolutions and seven UMAP routes are both present, resolution and UMAP decisions are kept separate, the long arms have cell-level evidence, and all comparison figures pass visual inspection.
