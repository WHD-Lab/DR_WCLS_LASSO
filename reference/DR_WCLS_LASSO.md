# DR_WCLS_LASSO

Fit a basic doubly robust weighted centered least squares (DR-WCLS)
model with variable selection.

## Usage

``` r
DR_WCLS_LASSO(
  data,
  fold,
  ID,
  time,
  Ht,
  St,
  At,
  prob,
  outcome,
  method_pseu,
  lam = NULL,
  noise_scale = NULL,
  splitrat = 0.8,
  virtualenv_path = "",
  beta = NULL,
  level = 0.9,
  core_num = NULL,
  max_tol = 10^{
-3
 },
  varSelect_program = "Python",
  standardize_x = TRUE,
  standardize_y = TRUE
)
```

## Arguments

- data:

  A data frame with one row per decision time (raw data; no
  pseudo-outcomes).

- fold:

  Number of folds used when estimating nuisance functions for the
  pseudo-outcome.

- ID:

  Column name of the participant identifier.

- time:

  Column name of the time-in-study (decision point).

- Ht:

  A vector specifying history features.

- St:

  A vector specifying moderator features.

- At:

  Column name of the treatment indicator.

- prob:

  Column name of the design probability \\p_t(A_t=1 \mid H_t)\\.

- outcome:

  Column name of the outcome.

- method_pseu:

  ML method used to estimate nuisance functions for the pseudo-outcome.
  One of `"CVLASSO"`, `"RandomForest"`, `"GradientBoosting"`.

- lam:

  Penalty value for randomized LASSO; if `NULL`, a default is used.
  Default is \\\sqrt{2n\* logp} \rho sd(y)\\ where \\\rho\\ is the split
  rate and \\n\\ is the number of rows.

- noise_scale:

  Gaussian noise added to the objective. Default is
  \\\sqrt{\frac{1-\rho}{\rho}\\ n}\\\mathrm{sd}(y)\\ where \\\rho\\ is
  the split rate and \\n\\ is the number of rows.

- splitrat:

  Data splitting rate \\\rho\\; used only if `noise_scale` or `lam` is
  `NULL`.

- virtualenv_path:

  Path to a Python virtual environment (for `reticulate`) when
  `varSelect_program = "Python"`.

- beta:

  True coefficients (for simulation use only).

- level:

  Confidence level (e.g., `0.90` for a 90% interval).

- core_num:

  Number of cores to use for parallel computation when compute
  pseudo-outcome.

- max_tol:

  Maximum tolerance for the pivot error. Default \\10^{-3}\\.

- varSelect_program:

  `"Python"` (requires a valid `virtualenv_path`) or `"R"`.

- standardize_x:

  Logical flag for design matrix standardization, prior to the model
  selection.

- standardize_y:

  Logical flag for outcome standardization, prior to the model
  selection.

## Value

- E:

  Selected variables.

- GEE_est:

  GEE estimates without adjusting for selection events.

- lowCI:

  lower bound of confidence interval.

- upperCI:

  upper bound of confidence interval.

- prop_low:

  the exact quantile for the lower bound of confidence interval.

- prop_up:

  the exact quantile for the upper bound of confidence interval.

- p_value:

  P-values for selected variables.

- post_true:

  condition on the selection events, the true values for parameter
  \\\beta_E\\ if simulation is conducted and true \\\beta\\ values are
  provided; Otherwise, this value will not be present

- true_signal:

  logical value indicating whether the selected parameter is one of the
  true signals if simulation is conducted and true \\\beta\\ values are
  provided; Otherwise, this value will not be present

## Details

The function generates a pseudo-outcome, performs variable selection,
then conducts post-selective inference to obtain valid confidence
intervals adjusted for data dependent model selection.

## Examples

``` r
  sim_data = generate_dataset(N = 1000, T = 40, P = 50, sigma_residual = 1.5, sigma_randint = 1.5, main_rand = 3, rho = 0.7,
  beta_logit = c(-1, 1.6 * rep(1/50, 50)), model = ~ state1 + state2 + state3 + state4,
  beta = matrix(c(-1, 1.7, 1.5, -1.3, -1),ncol = 1),
  theta1 = 0.8)
  Ht = unlist(lapply(1:50, FUN = function(X) paste0("state",X)))
  St = unlist(lapply(1:25, FUN = function(X) paste0("state",X)))

  UI_return = DR_WCLS_LASSO(data = sim_data,
  fold = 5, ID = "id",
  time = "decision_point",
  Ht = Ht, St = St, At = "action",
  prob = "prob", outcome = "outcome",
  method_pseu = "CVLASSO", lam = NULL, noise_scale = NULL, splitrat = 0.7,
  varSelect_program = "R", standardize_x = F, standardize_y = F)
#> Loading required package: parallel
#> [1] "remove 0 lines of data due to NA produced for yDR"
#> [1] "The current lambda value is: 386.489363835013"
#> [1] "select predictors: (Intercept)" "select predictors: state1"     
#> [3] "select predictors: state2"      "select predictors: state3"     
#> [5] "select predictors: state4"     
#> [1] FALSE
#>  [1] "state5"  "state6"  "state7"  "state8"  "state9"  "state10" "state11"
#>  [8] "state12" "state13" "state14" "state15" "state16" "state17" "state18"
#> [15] "state19" "state20" "state21" "state22" "state23" "state24" "state25"
#> Loading required package: dplyr
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union



```
