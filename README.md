## Code_PLS_test

## Overview 
[A Two-Step Projection-Based Goodness-of-Fit Test for Ultra-High Dimensional Sparse Regressions]( https://arxiv.org/abs/2412.10721).
This repository contains the R source code to reproduce the Monte Carlo simulations and real data analyses presented in the manuscript and its Supplementary Material. This paper proposes a two-step strategy that first constructs multiple test statistics based on projected predictors from distinct projections to mitigate the Curse of Dimensionality problem. The second step employs p-value combination methods (e.g., the minimum p-value and the Fisher combination) to form the final test statistics. This projection-based approach significantly mitigates the dimensionality problem, enabling the tests to detect local alternatives converging to the null at the rate as if the predictor were univariate.
## Overview of R Scripts and Correspondence to Manuscript Tables
| Script Name | Description | Corresponding Table/Figure | Runtime |
| :--- | :--- | :--- | :--- |
| `simulation_PLS_linear_fisher_min.R` | Main simulation study for high-dimensional **linear regression** models using the Fisher combination/minimum p-value approach. | Table 1, Tables S1-S31 | 552 h |
| `simulation_PLS_logit_fisher_min.R` | Main simulation study for high-dimensional **logistic regression** models using the Fisher combination/minimum p-value approach. | Table 2, Tables 32-38 | 334 h |
| `realdata_linear_data_test.R` | Real data applications on the linear regression models (e.g., Communities and Crime dataset (CommViolPredUnnormalizedData.xlsx)). | Figure 1 | 2.88min  |
| `realdata_logistic_data_test.R` | Real data applications on the logistic regression models (e.g., musks dataset (clean1.csv)). | - | 3.2min |
| `renv.lock` | R package dependency lock file used to ensure the reproducibility of the computational environment. | - | - |
##  Computational Time and Quick Test Settings

Please note that running the full simulation grids for the main tables can be highly computationally intensive (e.g., the linear regression simulation takes approximately 552 hours). To facilitate easier navigation and faster verification, we provide the following recommendations for a **"quick test"**:

1. **Reduce Monte Carlo Replications**: You can significantly reduce the runtime by setting a smaller number of Monte Carlo replications (e.g., setting `s = 100` or lower in the simulation scripts) to quickly test the pipeline and ensure the computational environment is correctly set up.
2. **Parallel Computing**: For full-scale simulations, we highly recommend utilizing parallel computing. You can increase the number of CPU cores used (e.g., via the `parallel` or `future` packages in R) to speed up the execution.
3. **Intermediate Results**: Since the full simulation grid takes a long time to run, the intermediate result CSVs used to generate the tables in the manuscript have been included in the repository.

### Data Description and Preprocessing

`CommViolPredUnnormalizedData.xlsx` is a cleaned version of the raw dataset, which is from  the UCI Machine Learning Repository(‘https://archive.ics.uci.edu/dataset/211/communities+and+crime+unnormalized’). The first five columns of the raw data cannot be used directly as covariates and were therefore removed. Columns 104–120 (inclusive) contain a large number of missing values; using them directly as covariates would be problematic, so these columns were also removed. Accordingly, the only difference between `CommViolPredUnnormalizedData.xlsx` and the raw dataset is that these columns have been deleted. Column indices refer to the original dataset before any columns were removed. Subsequent data processing steps can be found in `realdata_linear_test.R`.

`clean1.csv` is a cleaned version of the raw dataset, which is obtained from the UCI Machine Learning Repository (`https://archive.ics.uci.edu/dataset/74/musk+version+1`). The first two columns of the raw dataset are `molecule_name` and `conformation_name`, which are not predictors and were therefore removed. The resulting dataset contains 168 features and 1 response variable. Accordingly, the only difference between `clean1.csv` and the raw dataset is the removal of these two columns. Subsequent data processing solely involves data standardization. Details can be found in `realdata_logit_test.R`.
