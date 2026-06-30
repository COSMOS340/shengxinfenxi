# WGCNA and Machine Learning Diagnostic Workflows

Use this reference when a bioinformatics reproduction plan combines WGCNA, differential expression, feature selection, diagnostic modeling, ANN, SHAP, or multi-model benchmarking.

## Route Overview

A defensible route keeps discovery, feature selection, modeling, validation, and interpretation separate:

1. Build a metadata-verified expression matrix.
2. Run differential expression with a fixed contrast and multiple-testing correction.
3. Run WGCNA on a cleaned expression matrix and verified traits.
4. Intersect differential genes, module genes, and any prespecified theme gene set.
5. Split training and validation before feature selection.
6. Run feature selection or model benchmarking on training data only.
7. Lock features and model choice before validation.
8. Evaluate locked features in internal and external validation datasets.
9. Explain the final model with SHAP or importance analysis after validation design is fixed.

## WGCNA Guardrails

Record:

1. Expression matrix source and filtering rule.
2. Sample label source.
3. Sample outlier decisions.
4. Soft-thresholding power.
5. Network type.
6. Minimum module size.
7. Module merge threshold.
8. Module-trait correlation method.
9. Full module gene tables.

Do not use sample-name suffixes as the final source of phenotype labels. Use metadata, clinical annotations, or manually verified sample sheets.

## Feature Selection Routes

### Three-Method Route

Use this for a compact diagnostic paper:

1. LASSO with binomial `glmnet` and cross-validation.
2. Random forest with recorded tree count and importance metric.
3. SVM-RFE with recorded folds and feature ranking rule.
4. Intersect or rank genes using a rule declared before validation.

Checks:

1. Feature selection runs only on training data.
2. Random seeds and folds are saved.
3. The final gene list is locked before validation.
4. Top-N rules are declared before ROC analysis.

### Nine- or Twelve-Model Route

Use this when the paper compares common classifiers.

Common models:

1. RF.
2. SVM.
3. XGBoost.
4. GLM.
5. GBM.
6. KNN.
7. NNET.
8. LASSO or glmnet.
9. Decision tree.
10. Naive Bayes.
11. AdaBoost.
12. Bagging.

Checks:

1. Use the same split, folds, positive class, and performance metric for all models.
2. Save residual plots, ROC curves, feature-importance tables, and the full performance table.
3. Do not select the final model by training AUC alone.

### Large Multi-Model Route

Use this only when the project can report the full model list and validation design.

Required records:

1. The exact model-combination file.
2. The preprocessing and feature-selection methods for every combination.
3. The minimum selected-variable rule.
4. The model object files.
5. Risk score matrix.
6. Class prediction matrix.
7. Selected feature table.
8. AUC matrix across all datasets.
9. Final model selection rule.

Checks:

1. Remove any code that adds random noise to training or validation matrices.
2. Keep external validation outside model selection and ranking.
3. Report the full AUC matrix, not only the highest-scoring model.
4. Keep the final model mode explicit, such as plain model, logistic ensemble, or stacked model.

## ANN Route

Use ANN after the feature list is locked.

Record:

1. Input features.
2. Training data.
3. Hidden-layer setting.
4. Activation or default package behavior.
5. Seed.
6. Prediction-score output.
7. Training ROC and validation ROC.

Checks:

1. Do not use ANN to discover features after validation data have been inspected.
2. Keep prediction score direction consistent with ROC class direction.
3. Report validation performance separately from training performance.

## SHAP Route

Use SHAP to explain the final locked model.

Required records:

1. Model chosen for SHAP.
2. Training/test split.
3. Prediction probability column.
4. SHAP engine.
5. Feature-importance table.
6. Dependence, beeswarm, bar, waterfall, or force plots used.

Checks:

1. SHAP explains a prediction model; it does not establish mechanism.
2. SHAP ranking must not replace external validation.
3. The positive class must match ROC and manuscript wording.

## ROC and Validation

For every ROC result, record:

1. Dataset role.
2. Sample labels.
3. Positive class.
4. Prediction score.
5. AUC and confidence interval.
6. Whether the ROC is single-gene or multigene model performance.

Do not describe a single-gene ROC as a multigene diagnostic model.

## Output Checklist

Create these files for a manuscript-ready workflow:

1. `data_manifest.tsv`.
2. `methods_log.tsv`.
3. Full differential-expression table.
4. WGCNA module gene tables.
5. Feature-selection tables for each method.
6. Model performance matrix.
7. Locked feature list.
8. Validation ROC tables.
9. Figure manifest linking each panel to a script and source table.

## Red Flags

Fix these before reporting results:

1. Labels parsed only from sample names.
2. Feature selection run before validation split.
3. External validation used to select the model.
4. Added random noise in model matrices.
5. AUC reported without class direction.
6. Only the best model shown.
7. Generated figure values that are not derived from real data.
8. Hard-coded local paths inside reusable scripts.
