
# DR_WCLS_LASSO

<!-- badges: start -->

<!-- badges: end -->

## Introduction

Micro-randomized trials (MRTs) are designed to evaluate the
effectiveness of mobile health (mHealth) interventions delivered via
smartphones. In practice, the assumptions required for MRTs are often
difficult to satisfy: randomization probabilities can be uncertain,
observations are frequently incomplete, and prespecifying features from
high-dimensional contexts for linear working models is also challenging.
To address these issues, the **doubly robust weighted centered least
squares (DR-WCLS)** framework provides a flexible procedure for variable
selection and inference. The methods incorporates supervised learning
algorithms and enables valid inference on time-varying causal effects in
longitudinal settings.

## Installation

You can install the development version of MRTpostInfLASSO via github:

    require("remotes")
    remotes::install_github("WHD-Lab/DR_WCLS_LASSO")

Then run:

``` r
library(devtools)
install()
```

## Set-up

To configure a Python virtual environment in R, please run the following
code:

``` r
library(MRTpostInfLASSO)

# Configure virtual environment
venv_info = venv_config()
venv = venv_info$hash
print(venv)
# [1] "a9c268bc"

# Do the python deps load?
library(reticulate)
np = import("numpy", convert = FALSE)
lasso_mod = import("selectinf.randomized.lasso", convert = FALSE)$lasso

# Simple test
run_simple_lasso_test_snigdha_fixed(venv)
```

Or, if already configured, use

``` r
library(MRTpostInfLASSO)
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'MRTpostInfLASSO'
#> Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
#> loading 'MRTpostInfLASSO'
#> Warning: replacing previous import 'dplyr::slice' by 'xgboost::slice' when
#> loading 'MRTpostInfLASSO'
library(reticulate)
use_virtualenv("a9c268bc", required = TRUE)
np = import("numpy", convert = FALSE)
lassopy = import("selectinf.randomized.lasso", convert = FALSE)$lasso
selected_targets = import("selectinf.base", convert = FALSE)$selected_targets  
const = lassopy$gaussian
exact_grid_inference = import("selectinf.randomized.exact_reference", convert = FALSE)$exact_grid_inference  
```

## Detailed tutorial

A detailed tutorial and parameter explanation can be found here
<a href="https://whd-lab.github.io/DR_WCLS_LASSO/index.html" target="_blank" rel="noopener noreferrer">
https://whd-lab.github.io/DR_WCLS_LASSO/index.html </a>.
