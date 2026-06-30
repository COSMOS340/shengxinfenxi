# Figure Narrative Templates

Use this reference to turn a bioinformatics analysis route into manuscript figures. It focuses on what each figure proves, which source table supports it, and where claim boundaries sit.

## Core Principle

Main figures should advance the argument, not list every analysis that was run. A manuscript-facing bioinformatics figure set usually covers six to eight moves:

1. Study design and dataset roles.
2. Data quality and expression preprocessing.
3. Discovery evidence such as DEG, WGCNA, MR, or target intersection.
4. Feature or gene convergence.
5. Model or survival validation.
6. Pathway, immune, or cellular interpretation.
7. Orthogonal support such as single-cell localization, spatial localization, genetic evidence, docking, molecular dynamics, or experiments.
8. Supplementary completeness and robustness checks.

If a layer is not supported by data, remove that figure instead of adding a loosely related analysis.

## Universal Figure Skeleton

| Figure | Purpose | Common panels | Required source |
|---|---|---|---|
| Figure 1 | Study design | Workflow, cohort table, dataset split, method chain | Data manifest and sample metadata |
| Figure 2 | Data quality and discovery | PCA, batch correction, volcano, heatmap | Processed expression and differential table |
| Figure 3 | Convergence | WGCNA modules, Venn or UpSet, PPI, hub ranking | Gene lists and module tables |
| Figure 4 | Model construction | LASSO, random forest, SVM-RFE, model benchmark, risk score | Training expression and feature outputs |
| Figure 5 | Validation | ROC, calibration, DCA, KM, timeROC, external validation | Prediction score table and validation metadata |
| Figure 6 | Interpretation | GO, KEGG, GSEA, GSVA, ssGSEA, immune infiltration, gene-score correlation | Enrichment tables or score matrices |
| Figure 7 | Orthogonal support | UMAP, spatial feature plot, MR, SMR, colocalization, docking, MD | Single-cell, spatial, QTL/GWAS, or structure outputs |
| Supplement | Completeness | Full tables, all ROC, all pathways, mapping loss, seeds, parameter sweeps | Logs and source tables |

## Route Templates

### MR/GWAS plus GEO plus Machine Learning

Use when the manuscript links genetic risk factors to disease expression and diagnostic genes.

| Figure | Panels |
|---|---|
| Figure 1 | GWAS exposure, GWAS outcome, GEO datasets, QTL resources, workflow |
| Figure 2 | MR forest, scatter, heterogeneity, pleiotropy, leave-one-out |
| Figure 3 | Colocalization summary, GEO DEG volcano and heatmap |
| Figure 4 | Immune infiltration and gene-immune correlations |
| Figure 5 | Model comparison, ROC, nomogram, calibration, DCA |
| Figure 6 | eQTL, pQTL, SMR, or gene-level genetic evidence |

Claim boundary: MR and colocalization support genetic evidence. They do not replace functional experiments.

### Single-Cell plus Bulk Diagnostic Model

Use when single-cell data prioritizes disease cell states and bulk data builds a diagnostic model.

| Figure | Panels |
|---|---|
| Figure 1 | Single-cell and bulk datasets, discovery and validation split |
| Figure 2 | Single-cell QC, UMAP, cell-type annotation, marker heatmap |
| Figure 3 | Disease cell-type markers, bulk DEG, intersection |
| Figure 4 | LASSO, random forest, SVM-RFE, feature convergence |
| Figure 5 | ROC, nomogram, calibration, DCA, external validation |
| Figure 6 | ssGSEA, immune scores, subtype analysis, gene-score correlations |
| Figure 7 | Experimental or orthogonal validation when available |

Claim boundary: single-cell data localizes expression. It does not prove diagnostic performance.

### Comorbidity WGCNA plus Theme Gene Set

Use when the topic combines disease, comorbidity, and a named biological theme such as aging, metabolism, ferroptosis, or inflammation.

| Figure | Panels |
|---|---|
| Figure 1 | Disease dataset, comorbidity dataset, theme gene source |
| Figure 2 | Disease DEG and comorbidity WGCNA module-trait heatmap |
| Figure 3 | Theme-gene intersection, enrichment, PPI |
| Figure 4 | Model comparison and ROC |
| Figure 5 | Nomogram, calibration, DCA, external validation |
| Figure 6 | Immune infiltration and gene-immune correlations |

Claim boundary: cross-dataset convergence is not same-patient comorbidity evidence.

### Treatment Response Modeling

Use when response or resistance is the central phenotype.

| Figure | Panels |
|---|---|
| Figure 1 | Treatment cohort, response definition, time points |
| Figure 2 | WGCNA, DEG, and response-associated gene intersection |
| Figure 3 | Hub genes, feature selection, ANN or model training |
| Figure 4 | Clinical score association, ROC, tissue validation |
| Figure 5 | Immune infiltration and single-cell localization |
| Figure 6 | GSEA, regulatory network, genetic or structural extension |

Claim boundary: clinical response groups must be defined before expression analysis.

### TCGA Prognosis plus Single-Cell or Spatial Localization

Use when TCGA identifies prognostic genes and other data localizes them.

| Figure | Panels |
|---|---|
| Figure 1 | TCGA cohort, gene set, single-cell and spatial datasets |
| Figure 2 | TCGA DEG, pathway-gene intersection, enrichment |
| Figure 3 | Cox or LASSO model, risk score, KM, timeROC |
| Figure 4 | Nomogram, calibration, independent prognostic forest |
| Figure 5 | Single-gene expression and survival |
| Figure 6 | Single-cell localization and trajectory |
| Figure 7 | Spatial localization or experimental validation |

Claim boundary: localization data supports where a gene is expressed, not the functional role by itself.

### Network Toxicology plus Docking, Molecular Dynamics, and AOP

Use for compound-disease mechanism studies.

| Figure | Panels |
|---|---|
| Figure 1 | Compound structure, target databases, disease databases, workflow |
| Figure 2 | Target intersection, PPI, GO, KEGG |
| Figure 3 | GEO expression validation or diagnostic model |
| Figure 4 | Immune infiltration and gene-immune correlations |
| Figure 5 | Docking poses and binding energies |
| Figure 6 | MD metrics such as RMSD, RMSF, Rg, SASA, FEL |
| Figure 7 | AOP or integrated mechanism schematic |

Claim boundary: docking and molecular dynamics support structural plausibility. They do not prove clinical effect.

### Multi-Cohort WGCNA plus Model Interpretation

Use when the manuscript uses multiple GEO cohorts, WGCNA, large model panels, and SHAP.

| Figure | Panels |
|---|---|
| Figure 1 | Dataset split and workflow |
| Figure 2 | Batch-corrected multi-cohort expression, PCA, DEG |
| Figure 3 | WGCNA module-trait heatmap and gene intersection |
| Figure 4 | Model benchmark and locked feature set |
| Figure 5 | External ROC and calibration |
| Figure 6 | SHAP summary and dependence plots |
| Figure 7 | Mechanism or structural follow-up |

Claim boundary: SHAP explains the fitted model. It does not prove biological causality.

### Bulk WGCNA plus Immune plus SMR/Colocalization

Use when expression-derived diagnostic genes are strengthened with immune context and genetic evidence.

| Figure | Panels |
|---|---|
| Figure 1 | Discovery, test, external validation, QTL and GWAS resources |
| Figure 2 | DEG and WGCNA |
| Figure 3 | PPI, model ranking, feature selection |
| Figure 4 | Train, test, and external ROC |
| Figure 5 | Immune infiltration and gene-immune correlations |
| Figure 6 | TF or regulatory network |
| Figure 7 | SMR, HEIDI, and colocalization |

Claim boundary: SMR and colocalization require tissue, ancestry, and instrument limitations in the caption or result text.

## Supplementary Figure Pattern

| Supplement type | Contents |
|---|---|
| Data QC | Sample inclusion, PCA before correction, missing values, mapping loss |
| Parameter selection | WGCNA soft power, consensus clustering CDF, LASSO lambda, SVM-RFE feature count |
| Full results | All DEGs, modules, pathways, immune cells, model metrics |
| Sensitivity checks | MR pleiotropy, leave-one-out, ROC direction, alternate thresholds |
| Boundary results | Unmapped genes, unavailable QTL lookup, non-significant score comparisons |
| Figure QA | Exported file format, figure dimensions, visual inspection notes |

## Figure-To-Table Contract

Every figure needs a source table.

| Figure type | Required source table |
|---|---|
| Workflow | Data accession and cohort role table |
| Volcano or heatmap | Differential expression table and sample metadata |
| Venn or UpSet | Input gene lists and intersection table |
| WGCNA heatmap | Module-trait correlation table |
| Feature selection | Feature output and seed log |
| ROC | Prediction score, true class, and score direction |
| Nomogram or calibration | Model formula and validation table |
| Enrichment | Enrichment result and input gene list or ranked list |
| Immune scoring | Score matrix and group metadata |
| Single-cell or spatial | Object metadata, annotation, and marker table |
| MR, SMR, colocalization | Harmonized variant table and result table |
| Docking or MD | Protein source, ligand source, docking log, trajectory metrics |

## Final Checks

1. Remove figures that do not support the thesis.
2. Keep validation figures before interpretation figures when model performance is central.
3. Keep pathway and immune figures after the gene set is locked.
4. Put parameter sweeps and full lists in supplement.
5. Visually inspect every exported figure before reporting completion.
