# Scientific QC

Use this reference before reporting any bioinformatics reproduction result as manuscript-ready.

## Data Source Checks

1. Record accession, database, organism, tissue, platform, and download date.
2. Verify sample groups from metadata, not sample names alone.
3. Check whether the dataset contains paired samples, repeated measures, treatment time points, or mixed tissues.
4. Confirm genome build, probe annotation, gene symbol mapping, and duplicated-gene handling.
5. Keep raw files unchanged and write processed files separately.

## Expression Matrix Checks

1. Confirm genes and samples are in the expected orientation.
2. Confirm row and column names stay unique after mapping.
3. Confirm metadata row order matches expression columns.
4. Inspect library size, density, PCA, and batch structure.
5. Record normalization, filtering, and transformation rules.

## Differential Analysis Checks

1. State the exact contrast direction.
2. Use multiple-testing correction for genome-wide testing.
3. Record log fold-change and adjusted P value thresholds.
4. Export the full result table, not only significant genes.
5. Avoid generated values in volcano, heatmap, or enrichment plots.

## WGCNA Checks

1. Use metadata-derived traits.
2. Record outlier removal decisions.
3. Record soft power, minimum module size, merge threshold, and network type.
4. Export module membership and gene significance tables.
5. Lock selected modules before feature selection.

## Machine Learning Checks

1. Split training and validation before feature selection.
2. Fit preprocessing parameters on training data only.
3. Keep feature-selection results from training data only.
4. Record seeds, folds, tuning grids, selected features, and final model.
5. Report validation AUC, confidence intervals, and class direction.
6. Use external validation for claims about generalization.

## MR and Colocalization Checks

1. Record exposure and outcome dataset identifiers.
2. Check ancestry, sample overlap, and phenotype definition.
3. Filter weak instruments and record the F statistic rule.
4. Run sensitivity analyses for heterogeneity, pleiotropy, and leave-one-out effects.
5. Interpret SMR and colocalization as genetic support, not direct experimental proof.

## Single-Cell Checks

1. Record per-sample cell counts before and after QC.
2. Check mitochondrial fraction, detected genes, detected UMIs, and doublet strategy.
3. Annotate clusters using markers and reference sources.
4. Avoid claiming cell proportion differences without a design that supports that inference.
5. Validate highlighted genes in the relevant cell types and datasets.

## Network Toxicology and Structural Checks

1. Record database names, versions, and query dates.
2. Separate database-derived target lists from expression-derived disease evidence.
3. Record PPI confidence thresholds and hub ranking methods.
4. Record protein structure source, ligand source, docking box, exhaustiveness, and scoring function.
5. Record molecular dynamics force field, solvent model, ions, temperature, pressure, timestep, and simulation length.

## Claim Boundary Checks

Before writing a conclusion, mark every claim as one of:

1. Association in transcriptomic data.
2. Diagnostic classification in a validation dataset.
3. Genetic support from MR, SMR, or colocalization.
4. Cell-type localization from single-cell or spatial data.
5. Structural support from docking or simulation.
6. Experimental validation from wet-lab data.

Use only the claim level that the evidence supports.
