# Penalised spline estimation of covariate-specific time-dependent ROC curves: Supplementary Materials and Software in R
The repository contains the `R`-code needed to replicate the simulation study outlined in Section 3 of the paper

   *Penalised spline estimation of covariate-specific time-dependent ROC curves*

by Maria Xose Rodriguez Alvarez and Vanda Inacio.

Available at: XXX

## Instructions
### Step 1: Download and install the `R`-package `CondTimeROC`
``` r
# From GitHub:
# install.packages("devtools")
devtools::install_github("mxrodriguezUVIGO/CondTimeROC")
```
### Step 2: Download the `R`-code files
 * `r-code-functions-simulations-cTimeROC.R`: Contains the `R` functions required to run the simulation study.
 * `r-code-simulations-m1.R`: Replicates the simulation study for Scenario I.
 * `r-code-simulations-m2.R`: Replicates the simulation study for Scenario II.
 * `r-code-simulations-m3.R`: Replicates the simulation study for Scenario III.

To run the simulation study, open any of the last three files and execute the code. **Note**: The simulations may be time-consuming.
