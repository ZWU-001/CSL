# Conformalized Super Learning for Uncertainty Quantification in Predictions

This repository contains the R code used in the case study in:

> Wu, Z., Leisen, F., Luque-Fernandez, M.A. and Rubio, F.J. (2026). Conformalized Super Learning for Uncertainty Quantification in Predictions. Submitted.

# Description

The data used in this study are drawn from the National Health and Nutrition Examination Survey (NHANES), August 2021–August 2023 cycle. We consider serum creatinine, `LBXSCR` (mg/dL), as the response variable, as it is a key biomarker of kidney function and a central component in the estimation of glomerular filtration rate (eGFR). Participants with kidney conditions are excluded because impaired renal function directly affects creatinine levels. Pregnant participants are also excluded, as creatinine levels are systematically lower during pregnancy (Davison and Hytten, 1974). The analysis includes 18 covariates, comprising 12 continuous and 6 categorical variables, and the resulting sample size is 5027. As a sensitivity analysis, observations with response values greater than 1.5 are removed, yielding a reduced sample size of 4932.

We apply the Conformalized Super Learner (CSL) framework to the log-transformed response. This approach combines conformal prediction (CP) with the Super Learner (SL) ensemble through majority vote to construct prediction intervals under both split and full conformal settings. The base learners are linear regression (LM), Least Absolute Shrinkage and Selection Operator (LASSO), generalized additive models (GAM), neural networks (NNET), random forests (RF), and generalized additive models for location, scale and shape (GAMLSS). The non-conformity score is taken as the absolute raw residual for LM, LASSO, GAM, NNET, and RF, and as the absolute quantile residual for GAMLSS since it accounts for heteroscedasticity.

This case study has two main objectives:

1. Assess the performance of Split-CSL and Full-CSL in terms of empirical coverage and average interval width on the testing set.
2. Investigate, for three representative profiles with different serum creatinine levels, how the corresponding prediction intervals and point predictions change as one covariate varies.

# Folders and Files

The repository is organised according to the two main empirical tasks considered in the case study, under both Full-CSL and Split-CSL settings, together with benchmark and sensitivity analyses.

## `Full-CSL-Task1/`
Code and output files for **Task 1** under the full conformal setting, namely the evaluation of Full-CSL on the testing set in terms of empirical coverage and average interval width.

- `Full-CSL_Task1.R`: main script for Full-CSL testing-set evaluation.
- `Plot_Full_Task1.R`: script for plotting Full-CSL prediction intervals for ordered testing observations.
- `Results_Full_Task1_original.xlsx`: results for the main analysis based on the full data set.
- `Results_Full_Task1_removal.xlsx`: results for the sensitivity analysis after removing observations with `LBXSCR > 1.5`.


