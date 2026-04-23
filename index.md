# DR_WCLS_LASSO

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

``` r
require("remotes")
#> Loading required package: remotes
remotes::install_github("WHD-Lab/DR_WCLS_LASSO")
#> Using github PAT from envvar GITHUB_PAT. Use `gitcreds::gitcreds_set()` and unset GITHUB_PAT in .Renviron (or elsewhere) if you want to use the more secure git credential store instead.
#> Downloading GitHub repo WHD-Lab/DR_WCLS_LASSO@HEAD
#> cli         (3.6.5    -> 3.6.6   ) [CRAN]
#> cpp11       (0.5.2    -> 0.5.4   ) [CRAN]
#> vctrs       (0.6.5    -> 0.7.3   ) [CRAN]
#> tibble      (3.3.0    -> 3.3.1   ) [CRAN]
#> purrr       (1.2.0    -> 1.2.2   ) [CRAN]
#> magrittr    (2.0.4    -> 2.0.5   ) [CRAN]
#> lifecycle   (1.0.4    -> 1.0.5   ) [CRAN]
#> glue        (1.8.0    -> 1.8.1   ) [CRAN]
#> dplyr       (1.1.4    -> 1.2.1   ) [CRAN]
#> openssl     (2.3.4    -> 2.4.0   ) [CRAN]
#> curl        (7.0.0    -> 7.1.0   ) [CRAN]
#> xml2        (1.5.0    -> 1.5.2   ) [CRAN]
#> selectr     (0.4-2    -> 0.5-1   ) [CRAN]
#> httr        (1.4.7    -> 1.4.8   ) [CRAN]
#> rappdirs    (0.3.3    -> 0.3.4   ) [CRAN]
#> ps          (1.9.1    -> 1.9.3   ) [CRAN]
#> base64enc   (0.1-3    -> 0.1-6   ) [CRAN]
#> tinytex     (0.57     -> 0.59    ) [CRAN]
#> htmltools   (0.5.8.1  -> 0.5.9   ) [CRAN]
#> bslib       (0.9.0    -> 0.10.0  ) [CRAN]
#> yaml        (2.3.10   -> 2.3.12  ) [CRAN]
#> xfun        (0.54     -> 0.57    ) [CRAN]
#> highr       (0.11     -> 0.12    ) [CRAN]
#> processx    (3.8.6    -> 3.9.0   ) [CRAN]
#> rstudioapi  (0.17.1   -> 0.18.0  ) [CRAN]
#> rmarkdown   (2.30     -> 2.31    ) [CRAN]
#> knitr       (1.50     -> 1.51    ) [CRAN]
#> fs          (1.6.6    -> 2.1.0   ) [CRAN]
#> bit64       (4.6.0-1  -> 4.8.0   ) [CRAN]
#> vroom       (1.6.6    -> 1.7.1   ) [CRAN]
#> textshaping (1.0.4    -> 1.0.5   ) [CRAN]
#> systemfonts (1.3.1    -> 1.3.2   ) [CRAN]
#> backports   (1.5.0    -> 1.5.1   ) [CRAN]
#> tidyr       (1.3.1    -> 1.3.2   ) [CRAN]
#> broom       (1.0.10   -> 1.0.12  ) [CRAN]
#> timechange  (0.3.0    -> 0.4.0   ) [CRAN]
#> readr       (2.1.6    -> 2.2.0   ) [CRAN]
#> uuid        (1.2-1    -> 1.2-2   ) [CRAN]
#> gargle      (1.6.0    -> 1.6.1   ) [CRAN]
#> viridisLite (0.4.2    -> 0.4.3   ) [CRAN]
#> S7          (0.2.1    -> 0.2.2   ) [CRAN]
#> isoband     (0.2.7    -> 0.3.0   ) [CRAN]
#> data.table  (1.17.8   -> 1.18.2.1) [CRAN]
#> DBI         (1.2.3    -> 1.3.0   ) [CRAN]
#> blob        (1.2.4    -> 1.3.0   ) [CRAN]
#> Rcpp        (1.1.0    -> 1.1.1-1 ) [CRAN]
#> ragg        (1.5.0    -> 1.5.2   ) [CRAN]
#> lubridate   (1.9.4    -> 1.9.5   ) [CRAN]
#> ggplot2     (4.0.2    -> 4.0.3   ) [CRAN]
#> dtplyr      (1.3.2    -> 1.3.3   ) [CRAN]
#> dbplyr      (2.5.1    -> 2.5.2   ) [CRAN]
#> png         (0.1-8    -> 0.1-9   ) [CRAN]
#> zoo         (1.8-14   -> 1.8-15  ) [CRAN]
#> xgboost     (1.7.11.1 -> 3.2.1.1 ) [CRAN]
#> ranger      (0.17.0   -> 0.18.0  ) [CRAN]
#> reticulate  (1.44.1   -> 1.46.0  ) [CRAN]
#> Installing 56 packages: cli, cpp11, vctrs, tibble, purrr, magrittr, lifecycle, glue, dplyr, openssl, curl, xml2, selectr, httr, rappdirs, ps, base64enc, tinytex, htmltools, bslib, yaml, xfun, highr, processx, rstudioapi, rmarkdown, knitr, fs, bit64, vroom, textshaping, systemfonts, backports, tidyr, broom, timechange, readr, uuid, gargle, viridisLite, S7, isoband, data.table, DBI, blob, Rcpp, ragg, lubridate, ggplot2, dtplyr, dbplyr, png, zoo, xgboost, ranger, reticulate
#> Installing packages into '/private/var/folders/06/fdt8md_s73l4lbmj_4r6ydzm0000gn/T/RtmpmVn654/temp_libpathb085611eb9e8'
#> (as 'lib' is unspecified)
#> 
#> The downloaded binary packages are in
#>  /var/folders/06/fdt8md_s73l4lbmj_4r6ydzm0000gn/T//Rtmpxqmrtf/downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/06/fdt8md_s73l4lbmj_4r6ydzm0000gn/T/Rtmpxqmrtf/remotes135985d0961e3/WHD-Lab-DR_WCLS_LASSO-b42e4ba/DESCRIPTION’ ... OK
#> * preparing ‘MRTpostInfLASSO’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> Omitted ‘LazyData’ from DESCRIPTION
#> * building ‘MRTpostInfLASSO_1.0.0.tar.gz’
#> Installing package into '/private/var/folders/06/fdt8md_s73l4lbmj_4r6ydzm0000gn/T/RtmpmVn654/temp_libpathb085611eb9e8'
#> (as 'lib' is unspecified)
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
```

Or, if already configured, use

``` r
library(MRTpostInfLASSO)
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'MRTpostInfLASSO'
#> Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
#> loading 'MRTpostInfLASSO'
library(reticulate)
use_virtualenv("a9c268bc", required = TRUE)
```

## Virtual Environment Testing

``` r
# Do the python deps load?
library(reticulate)
np = import("numpy", convert = FALSE)
lasso_mod = import("selectinf.randomized.lasso", convert = FALSE)$lasso

# Simple test
run_simple_lasso_test_snigdha_fixed(venv)
```

## Detailed tutorial

A detailed tutorial and parameter explanation can be found here
[https://whd-lab.github.io/DR_WCLS_LASSO/index.html](https://whd-lab.github.io/DR_WCLS_LASSO/index.html)
.
