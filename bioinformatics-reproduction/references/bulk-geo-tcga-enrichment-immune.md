# Bulk GEO/TCGA, Enrichment, Immune Infiltration, and Subtyping

Use this reference when a reproduction task depends on bulk transcriptomics from GEO, TCGA, or both. It covers expression processing, identifier mapping, sample metadata checks, batch correction, differential expression, pathway analysis, immune scoring, and sample subtyping.

## Minimum Route

1. Lock the dataset design.
   - Record accession, organism, tissue, platform, sample count, phenotype definition, and analysis role.
   - Assign each dataset to discovery, training, validation, external validation, or mechanistic support before feature selection.

2. Build sample metadata.
   - Create a table with `Sample`, `Group`, `Dataset`, `Platform`, `Batch`, `Include`, and `ExcludeReason`.
   - Verify that expression matrix column names and metadata sample names match exactly.
   - Stop the analysis when any expression column has no metadata row.

3. Process identifiers.
   - GEO microarray: map probes to gene symbols using the exact GPL or platform annotation file.
   - TCGA: map Ensembl IDs to gene symbols using `rowData()` from `GDCprepare()` or a recorded GTF release.
   - Write mapping failures and duplicate-gene handling to a QC table.

4. Normalize and merge.
   - Check expression scale before log transformation.
   - Merge datasets only after each dataset has gene symbols and numeric expression.
   - Save matrices before and after batch correction.
   - Export PCA and boxplots before and after correction.

5. Run differential analysis.
   - Use a declared contrast such as `Disease - Control` or `Tumor - Normal`.
   - Write full differential results and filtered differential genes.
   - Record log fold-change threshold, multiple-testing method, and adjusted P-value threshold.

6. Add pathway or immune layers only after the expression matrix and grouping are locked.
   - ORA: GO/KEGG with identifier conversion reports.
   - GSEA: ranked full differential table, not only filtered genes.
   - GSVA/ssGSEA: cached GMT files and current GSVA API.
   - Immune infiltration: reference or gene set source, score matrix, and group comparison.

7. Add subtyping only when it answers a stated biological question.
   - Save consensus matrices, CDF or delta area evidence, cluster sizes, labels, and downstream phenotype comparisons.

## GEO Microarray Guardrails

### Required Files

Keep these files in the project scaffold:

| File | Purpose |
|---|---|
| Series matrix or expression table | Probe-level expression |
| GPL or platform annotation | Probe-to-gene mapping |
| sample metadata | Phenotype, tissue, platform, batch, inclusion |
| gene-symbol matrix | Processed expression |
| mapping QC table | Mapping losses and duplicate-gene decisions |

### Probe Mapping

1. Read the platform annotation header and use the exact gene-symbol field.
2. Record how multi-gene probe annotations are handled.
3. Aggregate duplicated gene symbols with one declared rule, such as mean expression or maximum-variance probe.
4. Avoid probe-pattern assumptions when the platform file provides the mapping.
5. Report input probe count, mapped probe count, gene count, and duplicate-gene count.

### Expression Scale

Use the source documentation and matrix range to decide transformation:

- If the matrix is raw intensity or values are large, transform with `log2(x + 1)` or the platform-specific rule.
- If the matrix is already log2-transformed, do not transform again.
- If multiple datasets are merged, perform scale checks per dataset before merging.

## TCGA Guardrails

### Preferred Download Pattern

Use `TCGAbiolinks` when the task can use GDC directly:

```r
query <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
se <- GDCprepare(query)
expr <- assay(se, "unstranded")
gene_info <- rowData(se)
sample_info <- colData(se)
```

Record `project`, `data.category`, `data.type`, `workflow.type`, assay name, and GDC query date.

### Manual GDC Pattern

When using downloaded GDC files:

1. Keep the manifest and JSON metadata.
2. Parse sample IDs from metadata, not only from file names.
3. Build a gene-by-sample matrix from the same quantification column for every file.
4. Remove Ensembl version suffixes only after storing the raw ID.
5. Map Ensembl IDs to symbols using a recorded GTF or GDC row annotation.
6. Export unmapped IDs and duplicated symbols.

### Sample Type Code

Extract the two-digit TCGA sample type code from the barcode. For a tumor-normal comparison:

| Code | Use |
|---|---|
| `01` | Primary tumor |
| `11` | Solid tissue normal |

Do not treat every code beginning with `0` as tumor. Keep non-`01` and non-`11` samples in an exclusion table unless the study design names a different sample type.

## Batch Correction

Batch correction is valid only after phenotype labels and dataset identities are fixed.

Minimum outputs:

| Output | Purpose |
|---|---|
| `merged_before_batch_removal` | Audit raw merged expression |
| `merged_after_batch_removal` | Matrix for downstream bulk modules |
| `sample_batch_info` | Sample-to-batch mapping |
| PCA before and after | Batch structure review |
| Boxplots before and after | Distribution review |

Do not mix before-correction and after-correction matrices across downstream steps.

## Differential Expression

Use `limma` for log-scale microarray or transformed expression:

```r
group <- factor(metadata$Group, levels = c("Control", "Disease"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr, design)
contrast <- makeContrasts(Disease - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))
deg <- topTable(fit2, adjust.method = "fdr", number = Inf)
```

Record:

1. Contrast direction and meaning of positive logFC.
2. Multiple-testing method.
3. Filter thresholds.
4. Number of tested genes.
5. Number of upregulated and downregulated genes.
6. Full table path and filtered table path.

## Enrichment Analysis

### GO and KEGG ORA

Use ORA for a locked gene set, such as differential genes or intersection genes:

1. Convert symbols to Entrez IDs.
2. Export converted and failed genes.
3. Run GO and KEGG with declared organism.
4. Filter by P-value and adjusted P-value.
5. Plot only from the result table used for reporting.

### GSEA

Use GSEA with a full ranked gene list:

1. Start from the complete differential table.
2. Remove duplicated gene symbols with a declared rule.
3. Rank by logFC, signal-to-noise, or another declared statistic.
4. Use cached GMT files or recorded MSigDB releases.
5. Report NES, adjusted P-value, core enrichment, and gene set version.

## GSVA and ssGSEA

Use current GSVA constructors:

```r
library(GSVA)
library(GSEABase)

gene_sets <- getGmt("immune_or_pathway.gmt", geneIdType = SymbolIdentifier())
param <- ssgseaParam(exprData = expr, geneSets = gene_sets)
scores <- gsva(param, BPPARAM = BiocParallel::SerialParam())
```

For GSVA:

```r
param <- gsvaParam(expr, gene_sets, minSize = 5, maxSize = 500, kcdf = "Gaussian")
scores <- gsva(param)
```

Guardrails:

1. Use gene symbols that match the GMT.
2. Cache GMT files in the project.
3. Normalize scores only with a documented rule.
4. Compare scores only after sample groups are verified.
5. Correct for multiple tests across pathways or immune cell types.

## Immune Infiltration

Use immune scoring as transcriptomic evidence. It does not prove direct cell counts unless paired with orthogonal validation.

Supported patterns:

| Method | Required input | Output |
|---|---|---|
| ssGSEA immune gene sets | Expression matrix and immune GMT | Immune score matrix |
| CIBERSORT-style SVR | Expression matrix and reference signature | Cell fraction table with quality metrics |

For CIBERSORT-style workflows:

1. Intersect reference genes and expression genes.
2. Decide whether quantile normalization is required by platform.
3. Keep P-value, correlation, and RMSE columns.
4. Filter low-quality samples only after writing the full result.
5. Run training and validation datasets separately.

## Consensus Subtyping

Use subtyping when it will be interpreted biologically or clinically.

Minimum records:

1. Input feature matrix and feature source.
2. `maxK`, final K, clustering algorithm, distance metric, resampling rate, repetitions, and seed.
3. Consensus matrices for each K.
4. CDF or delta area summary.
5. Cluster sizes.
6. Sample-to-subtype label table.
7. Downstream pathway, immune, clinical, or survival comparison.

Do not fix K from a tutorial setting alone. K must be supported by clustering evidence and downstream interpretability.

## Generated-Data Guardrail

Tutorial scripts often include random toy data for plotting demonstrations. Remove these blocks from reproduction scripts or disable them with an explicit stop. Generated toy values must not enter reported figures, tables, feature lists, or supplementary files.

## Minimal Project Scaffold

```text
project/
  metadata/
    sample_metadata.tsv
    mapping_qc.tsv
    sample_exclusions.tsv
  data_raw/
    geo/
    tcga/
    gene_sets/
  data_processed/
    expression/
    batch_correction/
  results/
    differential_expression/
    enrichment/
    immune/
    subtype/
  figures/
  logs/
    sessionInfo.txt
    run_parameters.yml
```

## Report Checklist

Before reporting the analysis as complete:

1. State exact datasets and accessions.
2. State exact group definitions and sample exclusions.
3. State identifier mapping source and mapping loss.
4. State expression scale and normalization.
5. State contrast direction and thresholds.
6. State gene set or immune reference source.
7. State validation design.
8. Confirm figure visual inspection.
