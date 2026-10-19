#Code_PLS_test

## Overview 
[A Two-Step Projection-Based Goodness-of-Fit Test for Ultra-High Dimensional Sparse Regressions]( https://arxiv.org/abs/2412.10721).
This repository contains the R source code to reproduce the Monte Carlo simulations and real data analyses presented in the manuscript and its Supplementary Material. This paper proposes a two-step strategy that first constructs multiple test statistics based on projected predictors from distinct projections to mitigate the Curse of Dimensionality problem. The second step employs p-value combination methods (e.g., the minimum p-value and the Fisher combination) to form the final test statistics. This projection-based approach significantly mitigates the dimensionality problem, enabling the tests to detect local alternatives converging to the null at the rate as if the predictor were univariate.

## Overview of R Scripts and Correspondence to Manuscript Tables
| Script Name | Description | Corresponding Table/Figure |
| :--- | :--- | :--- |
| `simulation_PLS_linear_fisher_min.R` | Main simulation study for high-dimensional **linear regression** models using the Fisher combination/minimum p-value approach. | [ ] |
| `simulation_PLS_logit_fisher_min.R` | Main simulation study for high-dimensional **logistic regression** models using the Fisher combination/minimum p-value approach. | [ ] |
| `realdata_linear_data_test.R` | Real data applications on the linear regression models (e.g., Communities and Crime dataset(communities.csv)). | [ ] |
| `realdata_logistic_data_test.R` | Real data applications on the logistic regression models (e.g., OHSU AML dataset). | [ ] |
| `renv.lock` | R package dependency lock file used to ensure the reproducibility of the computational environment. | - |



