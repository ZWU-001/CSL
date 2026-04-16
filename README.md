# Conformalized Super Learning for Uncertainty Quantification in Predictions

This repository contains the R code used in the case study in:

> Wu, Z., Leisen, F., Luque-Fernandez, M.A. and Rubio, F.J. (2026). Conformalized Super Learning for Uncertainty Quantification in Predictions. Submitted.

# Description

The data used in this study are drawn from the National Health and Nutrition Examination Survey (NHANES), August 2021–August 2023 cycle. We consider serum creatinine, `LBXSCR` (mg/dL), as the response variable, as it is a key biomarker of kidney function and a central component in the estimation of glomerular filtration rate (eGFR). Participants with kidney conditions are excluded because impaired renal function directly affects creatinine levels. Pregnant participants are also excluded, as creatinine levels are systematically lower during pregnancy (Davison and Hytten, 1974). The analysis includes 18 covariates, comprising 12 continuous and 6 categorical variables, and the resulting sample size is 5027. As a sensitivity analysis, observations with response values greater than 1.5 are removed, yielding a reduced sample size of 4923.

We apply the Conformalized Super Learner (CSL) framework to the log-transformed response, and the prediction intervals will be transformed back to the original scale at the end. This approach combines conformal prediction (CP) with the Super Learner (SL) ensemble through majority vote to construct prediction intervals under both split and full conformal settings. The base learners are linear regression (LM), Least Absolute Shrinkage and Selection Operator (LASSO), generalized additive models (GAM), neural networks (NNET), random forests (RF), and generalized additive models for location, scale and shape (GAMLSS). The non-conformity score is taken as the absolute raw residual for LM, LASSO, GAM, NNET, and RF, and as the absolute quantile residual for GAMLSS since it accounts for heteroscedasticity.

This case study has two main objectives:

1. Assess the performance of Split-CSL and Full-CSL in terms of empirical coverage and average interval width on the testing set.
2. Investigate, for three representative profiles with different serum creatinine levels, how the corresponding prediction intervals and point predictions change as one covariate varies.

# Folders and Files

The repository is organised according to the two main empirical tasks considered in the case study, under both Full-CSL and Split-CSL settings, together with benchmark and sensitivity analyses.

## `Full-CSL-Task1/`
Files for the testing-set performance assessment under the full conformal setting.

- `full_csl_test_performance.R`: main script for evaluating Full-CSL on the testing set.
- `results_full_csl_test_full_data.xlsx`: results for Full-CSL based on the full data set.
- `results_full_csl_test_outlier_removed.xlsx`: results for Full-CSL after removing observations with response values greater than 1.5.
- `plot_full_csl_test_intervals.R`: script for plotting Full-CSL prediction intervals for ordered testing observations.

## `Full-CSL-Task2/`
Files for the covariate-specific analysis under the full conformal setting.

- `full_csl_profiles_continuous.R`: main script for the analysis of continuous covariates across three representative profiles.
- `full_csl_profiles_categorical.R`: main script for the analysis of categorical covariates across three representative profiles, including the corresponding plots.
- `results_full_csl_profiles_full_data.xlsx`: results for the continuous covariate-specific analysis based on the full data set.
- `results_full_csl_profiles_outlier_removed.xlsx`: results for the continuous covariate-specific analysis after removing observations with response values greater than 1.5.
- `plot_full_csl_profiles.R`: script for plotting Full-CSL covariate-specific results for continuous covariates.

## `Split-CSL-Task1/`
Files for the testing-set performance assessment under the split conformal setting.

- `split_csl_test_performance.R`: main script for evaluating Split-CSL on the testing set, including the corresponding plots.

## `Split-CSL-Task2/`
Files for the covariate-specific analysis under the split conformal setting.

- `split_csl_profiles_continuous.R`: script for the analysis of continuous covariates across three representative profiles, including the corresponding plots.
- `split_csl_profiles_categorical.R`: script for the analysis of categorical covariates across three representative profiles, including the corresponding plots.

## Additional files

- `classical_pi_benchmark.R`: benchmark analysis based on classical linear-model prediction intervals.
- `split_rule_sensitivity.R`: sensitivity analysis for alternative training-calibration split ratios in Split-CSL.
- `tree_plots.R`: scripts for variable-importance and tree-based plots used to interpret the dominant random forest learner.



