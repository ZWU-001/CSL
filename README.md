# CSL: Conformalized Super Learning for Uncertainty Quantification in Predictions

## Overview

This repository contains the R code and results used to reproduce the case
study in:

> Wu, Z., Leisen, F., Luque-Fernandez, M.A. and Rubio, F.J. (2026).
> Conformalized Super Learning for Uncertainty Quantification in Predictions.
> *Submitted.*

The **Conformalized Super Learner (CSL)** framework combines
**conformal prediction (CP)** with a **Super Learner (SL)** ensemble and the majority vote to
construct prediction intervals.

## Requirements

The following R packages are required:

```r
install.packages(c("SuperLearner", "glmnet", "mgcv", "nnet",
                   "randomForest", "gamlss", "conformalInference",
                   "ggplot2", "readxl", "writexl"))
```

> Verify exact dependencies against the `library()` calls at the top of each
> script, as the list above may be incomplete.

## Repository structure

```
CSL/
├── Full-CSL-Task1/
│   ├── full_csl_test_performance.R              # Full-CSL evaluation on test set
│   ├── plot_full_csl_test_intervals.R           # Plots of prediction intervals
│   ├── results_full_csl_test_full_data.xlsx     # Results: full dataset
│   └── results_full_csl_test_outlier_removed.xlsx  # Results: outliers removed (> 1.5)
│
├── Full-CSL-Task2/
│   ├── full_csl_profiles_continuous.R           # Continuous covariate profiles
│   ├── full_csl_profiles_categorical.R          # Categorical covariate profiles + plots
│   ├── plot_full_csl_profiles.R                 # Plots for continuous covariate results
│   ├── results_full_csl_profiles_full_data.xlsx # Results: full dataset
│   └── results_full_csl_profiles_outlier_removed.xlsx  # Results: outliers removed
│
├── Split-CSL-Task1/
│   └── split_csl_test_performance.R             # Split-CSL evaluation + plots
│
├── Split-CSL-Task2/
│   ├── split_csl_profiles_continuous.R          # Continuous covariate profiles + plots
│   └── split_csl_profiles_categorical.R         # Categorical covariate profiles + plots
│
├── classical_pi_benchmark.R                     # Benchmark: classical LM prediction intervals
├── split_rule_sensitivity.R                     # Sensitivity: alternative train/calibration splits
└── tree_plots.R                                 # Variable importance and tree-based plots
```

## Case study summary

| Task | Setting | Description |
|---|---|---|
| Task 1 | Full-CSL and Split-CSL | Empirical coverage and average interval width on the test set |
| Task 2 | Full-CSL and Split-CSL | Covariate-specific prediction intervals for three representative profiles |
| Benchmark | — | Classical linear model prediction intervals for comparison |
| Sensitivity | Split-CSL | Alternative training–calibration split ratios |

## Data

The analysis uses data from the **National Health and Nutrition Examination
Survey (NHANES), August 2021–August 2023 cycle**, publicly available at
[https://www.cdc.gov/nchs/nhanes](https://www.cdc.gov/nchs/nhanes). The data
are not bundled in this repository and should be downloaded directly from the
NHANES website. Participants with kidney conditions or pregnancy are excluded,
and the response (serum creatinine) is log-transformed prior to modelling.

## Citation

If you use this code, please cite:

```bibtex
@article{wu:2026,
  author  = {Wu, Z. and Leisen, F. and Luque-Fernandez, M.A. and Rubio, F.J.},
  title   = {Conformalized Super Learning for Uncertainty Quantification
             in Predictions},
  journal = {Submitted},
  year    = {2026}
}
```

