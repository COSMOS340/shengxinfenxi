# GSE205506 desktop completion status

Status: `COMPLETED_WITH_REVIEW_REQUEST`

Request: `20260716_stop_prjna932556_and_build_gse205506_request.md`

Verified request blob: `ff3a9abc0fcf9499fb6dff72a026b44afb646645`

## 1. PRJNA932556 stop and deletion

- Six exact paths were established from the existing request/log/process evidence before deletion.
- `134161587134` bytes were present and `134161587134` bytes were deleted.
- Post-deletion audit: 6/6 paths absent, 0 bytes remaining, 0 audit errors.
- Cell Ranger software, references, GSE205506, and unrelated projects were retained.
- I: free-space increase: `34677944320` bytes. F: free-space increase: `84182069248` bytes.

## 2. GSE205506 inputs

- GEO archive: `F:/pan-gi-ici-ccc-20260709/GSE205506/raw/GSE205506_RAW.tar`
- Actual GEO URL: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205506/suppl/GSE205506_RAW.tar`
- Archive size: `1585182720` bytes
- Archive SHA-256: `07A1D08CA5D90CD9E7D8336DDEA506273B3E6FACE605F04796CBAD4C7F82B7F3`
- Extracted structure: 120 files, 40 exact matrix/barcode/feature triplets.
- Table S1: `F:/pan-gi-ici-ccc-20260709/GSE205506/supplement/mmc2.xlsx`, 18,791 bytes, SHA-256 `67649F6FBA1C27773FF50EED566382FDE04D32D0CD5E7A17007009854C736D25`.
- Author marker workbook: `F:/pan-gi-ici-ccc-20260709/GSE205506/supplement/mmc3.xlsx`, 1,302,551 bytes, SHA-256 `28AA8BC71620AE3B667D0243891B631DA3B8EB1576108C88463B2F66979EAC80`.
- The publisher supplement link entered a redirect loop. The user supplied the original `mmc2.xlsx`; no download URL was inferred or fabricated.

## 3. Exact response mapping

- Table S1 title and table structure were asserted before parsing.
- 19/19 patients were parsed: 15 pCR and 4 non-pCR.
- Exact non-pCR subjects: `P12`, `P18`, `P26`, `P31`.
- GEO join: 40/40 samples and 19/19 subjects matched.
- Unmatched samples: 0; duplicate mapping subjects: 0; many-to-many joins: 0; blank response fields: 0.

## 4. R/Seurat object

- The formal object was built entirely in R from the 40 author-processed matrices.
- Final retained singlets: 238,934 cells and 33,538 features.
- Table S1 groups: 192,323 pCR cells and 46,611 non-pCR cells.
- Timepoints: 61,442 pre-treatment and 177,492 post-treatment cells.
- QC thresholds were derived per sample, with sensitivity, doublet, and sample-retention audits exported.
- Only filtered matrices were available; ambient RNA contamination was audited, but SoupX correction was not run without unfiltered droplets.
- Normalization: LogNormalize, scale factor 10,000; 3,000 variable genes; PCA 50 PCs.
- Harmony was run by exact `geo_accession` on PCs 1:30; uncorrected PCA was retained for audit.
- Mean same-sample 30-neighbor fraction changed from `0.2262533` (PCA) to `0.1146883` (Harmony).
- Clustering resolutions 0.2/0.4/0.6/0.8/1.0 produced 16/22/28/35/40 clusters.
- Working resolution: 0.6, with 28 clusters.
- Full-cell R Wilcoxon/AUC markers were exported as exactly 280 Top10 and 1,400 Top50 rows.

## 5. Annotation review and author marker comparison

- `mmc3.xlsx` contains 45 author clusters across five compartment-specific sheets: T/I/NK 19, B 7, Myeloid 10, Endothelial 6, Fibroblast 3.
- Current broad annotations were reviewed against canonical markers and the author Top50 lists.
- Strong same-compartment matches include B cluster 1 -> Bn (10 shared markers), Myeloid cluster 9 -> FOLR2+Mac (10), Plasma cluster 11 -> pB IgG (30), Fibroblast cluster 16 -> CXCL12+Fibro (10), Pericyte/SMC cluster 17 -> Myofibroblast (26), and Lymphatic endothelial cluster 22 -> LYVE1+LEC (18).
- Author overlap alone was not used to force subtype labels. For example, current CD4 T cluster 3 overlaps the author CD8+MAIT list by 14 genes, demonstrating that whole-dataset clustering and compartment-specific reclustering are not directly interchangeable.
- `mmc3.xlsx` has no epithelial sheet; epithelial clusters were not assigned an author subtype from this workbook.

## 6. QC issue requiring review

- Nine clusters (2, 6, 12, 15, 19, 21, 23, 24, 27) contain at least five mitochondrial genes in their Top10 markers.
- These clusters contain 53,708 cells (22.478% of the retained object).
- Median percent mitochondrial reads include 92.06% in cluster 2, 96.76% in cluster 12, and 99.82% in cluster 27.
- These cells remain in the object and are explicitly flagged; they were not silently removed.
- See `MAC_REVIEW_REQUEST.md` for the requested decision before subtype-level analysis.

## 7. Local-only final object and visual QA

- Final reviewed object: `F:/pan-gi-ici-ccc-20260709/GSE205506/r_analysis/gse205506_seurat_final_reviewed.rds`
- Size: `13973709204` bytes
- SHA-256: `3916592E5CCE0B9F490D304847BDBBA947E0C9500AFDCDA59F7D699893D3CE49`
- Large matrices and RDS files remain on the PC and are not uploaded to GitHub.
- Reviewed PNG/PDF figures were visually inspected for clipping, overlap, canvas size, legend overflow, and color contrast.
- Final palettes use blue/orange for pre/post, green/magenta for pCR/non-pCR, and distinct high-contrast broad cell-type colors.

