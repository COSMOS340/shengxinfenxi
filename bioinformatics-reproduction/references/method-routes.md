# Method Routes

Use this reference to select the smallest defensible analysis route for a bioinformatics reproduction task.

## Bulk GEO or TCGA Expression

Use when the central data type is bulk transcriptomics.

For detailed guardrails on GEO/TCGA expression processing, identifier mapping, batch correction, enrichment, immune scoring, and subtyping, read `bulk-geo-tcga-enrichment-immune.md`.

Minimum route:

1. Download expression matrix and platform annotation.
2. Map probes or Ensembl IDs to gene symbols using the platform file or authoritative annotation.
3. Aggregate duplicated gene symbols with a declared rule.
4. Verify sample labels from metadata.
5. Normalize or transform according to platform and source.
6. Run differential analysis with a fixed contrast and multiple-testing correction.
7. Export expression matrix, metadata, differential table, and QC plots.

QC focus:

1. Label source.
2. Batch and platform effects.
3. Gene ID mapping.
4. Filtering threshold.
5. Contrast direction.

## WGCNA

Use WGCNA to identify co-expression modules associated with verified traits.

For WGCNA routes that feed LASSO, random forest, SVM-RFE, ANN, SHAP, or large multi-model diagnostic workflows, read `wgcna-machine-learning.md` before coding.

Minimum route:

1. Use a cleaned expression matrix with samples as rows and genes as columns after transposition.
2. Remove low-quality samples and genes with WGCNA checks.
3. Cluster samples and compare outliers against metadata.
4. Select soft-thresholding power using scale-free topology and mean connectivity.
5. Build adjacency, TOM, dynamic tree modules, and merged modules.
6. Correlate module eigengenes with traits.
7. Export module gene lists and module-trait heatmaps.

QC focus:

1. Trait labels must come from metadata.
2. Sample count must support module analysis.
3. The selected power and module merge threshold must be recorded.
4. Module selection must be made before downstream model evaluation.

## Machine Learning Feature Selection

Use machine learning when the goal is diagnostic, classification, or feature prioritization.

For detailed guardrails on WGCNA plus machine learning diagnostic routes, read `wgcna-machine-learning.md`.

Supported route levels:

1. Three-method route: LASSO, random forest, and SVM-RFE.
2. Four-method route: add logistic regression, gradient boosting, or another declared model.
3. Nine- or twelve-model route: benchmark a fixed panel of algorithms using the same resampling design.
4. Large multi-model route: use algorithm combinations only when training, validation, and ranking rules are fixed before evaluation.

Minimum route:

1. Split discovery/training and validation before feature selection.
2. Fit feature selectors only on the training data.
3. Record seeds, folds, tuning grid, and final feature list.
4. Evaluate the locked feature set in validation data.
5. Report AUC with confidence interval and class direction.

QC focus:

1. Do not select genes using all datasets before validation.
2. Do not tune the model on external validation.
3. Do not add random noise or generated values to improve figures.
4. Keep gene order, sample order, and class labels checked at every matrix merge.

## ANN Diagnostic Model

Use an artificial neural network after a small feature list has been locked by upstream analysis.

Minimum route:

1. Use locked feature genes as input nodes.
2. Train only on the training set.
3. Record architecture, activation, seed, and training settings.
4. Export prediction scores.
5. Evaluate ROC in training and independent validation data.

QC focus:

1. ANN does not replace external validation.
2. Feature selection must precede ANN training and remain fixed.
3. Score direction must be verified before ROC reporting.

## MR, SMR, and Colocalization

Use genetic evidence when the analysis claims causal direction or gene-level support.

Minimum route:

1. Define exposure and outcome resources.
2. Harmonize alleles and remove weak instruments.
3. Run MR with sensitivity tests.
4. Run colocalization or SMR/HEIDI when gene-level evidence is needed.
5. Correct for multiple tests.

QC focus:

1. Instrument strength.
2. Ancestry and sample overlap.
3. Allele harmonization.
4. Horizontal pleiotropy.
5. HEIDI and colocalization interpretation.

## Single-Cell RNA-Seq

Use single-cell analysis to localize genes, define disease cell states, or support mechanisms.

For detailed guardrails on Seurat QC, SingleR annotation, target-gene localization, Monocle3 trajectory, spatial transcriptomics, CellChat communication, and scTenifoldKnk virtual perturbation, read `single-cell-spatial-trajectory.md`.

Minimum route:

1. Load count matrices with sample metadata.
2. Apply cell and gene QC thresholds.
3. Normalize, select variable features, scale, reduce dimensions, cluster, and annotate.
4. Validate annotation with marker genes and reference mapping when appropriate.
5. Test expression of bulk-derived genes in relevant cell types.

QC focus:

1. Sample identity and batch.
2. Mitochondrial fraction.
3. Doublets.
4. Marker specificity.
5. Overinterpretation of cell-type abundance without study design support.

## Spatial Transcriptomics, Trajectory, Communication, and Virtual Perturbation

Use these modules only when they answer a defined biological question beyond bulk differential expression.

For the detailed workflow, read `single-cell-spatial-trajectory.md`.

Minimum route:

1. Verify single-cell or spatial accession, platform, sample metadata, and tissue context.
2. Build a Seurat object with explicit QC and retained-cell counts.
3. Annotate cell types with reference mapping and marker support.
4. Use target-gene plots to localize locked genes from upstream analysis.
5. Use Monocle3 trajectory only with a declared root-state rule.
6. Use spatial analysis only after image-coordinate alignment is checked.
7. Use communication analysis only when ligand-receptor inference is actually run.
8. Use scTenifoldKnk outputs as virtual perturbation evidence.

QC focus:

1. Sample identity and batch effects.
2. Cell and spot filtering thresholds.
3. Doublet and contamination handling.
4. Cell-type annotation support.
5. Pseudotime root selection.
6. Spatial image alignment.
7. Communication model completeness.
8. Simulation wording for perturbation results.

## Immune Infiltration and Enrichment

Use immune and pathway scoring as interpretation layers.

For detailed guardrails on GO/KEGG, GSEA, GSVA, ssGSEA, CIBERSORT-style workflows, and ConsensusClusterPlus subtyping, read `bulk-geo-tcga-enrichment-immune.md`.

Minimum route:

1. Use normalized expression matrix with gene symbols.
2. Choose one scoring method and record the gene sets.
3. Compare scores between verified groups.
4. Correct for multiple testing when many cells or pathways are tested.
5. Correlate scores with locked genes only after validation design is clear.

QC focus:

1. Gene-set source.
2. Multiple-testing control.
3. Batch effects in score comparisons.
4. Claim wording: scores indicate transcriptomic signatures, not direct cell counts.

## Network Toxicology, Docking, and Molecular Dynamics

Use for compound mechanism studies.

For detailed guardrails on compound target integration, AutoDock Vina, GROMACS, free-energy landscapes, and AOP narratives, read `network-toxicology-docking-md.md`.

Minimum route:

1. Collect compound targets and disease genes from named databases.
2. Merge and deduplicate identifiers with documented mapping.
3. Intersect with expression-based disease evidence.
4. Run PPI, enrichment, docking, and simulation in a declared order.
5. Preserve docking parameters, protein structure source, ligand source, box settings, and simulation settings.

QC focus:

1. Database version and query date.
2. Identifier mapping.
3. PPI confidence threshold.
4. Docking reproducibility.
5. Simulation length and force-field settings.
