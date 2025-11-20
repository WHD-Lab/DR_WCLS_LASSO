# pesudo_outcome_generator_rf_v2

Cross-fitting is applied when generate the target pseudo-outcome. This
function uses Random Forest to conduct model training on folds then
estimates \\E\[Y\_{t+1}\|H_t, A_t\], E\[A_t\|H_t\], E\[A_t\|S_t\]\\ for
the reserved fold. Then the function returns a dataset with column named
"yDR" for DR-WCLS pseudo outcome.

## Usage

``` r
pseudo_outcome_generator_rf_v2(
  fold,
  ID,
  data,
  Ht,
  St,
  At,
  prob,
  outcome,
  core_num = NULL
)
```

## Arguments

- fold:

  number of folds to split when do corss-fitting

- ID:

  the name of column where participants' ID are stored

- data:

  dataset name

- Ht:

  a vector that contains column names of control variables

- St:

  a vector that contains column names of moderator variables; St should
  be a subset of Ht

- At:

  column names of treatment (At)

- prob:

  column names of \\p_t(A_t = 1\|H_t)\\, the experiment design treatment
  probability

- outcome:

  column names of outcome variable

- core_num:

  number of cores will be used for calculation

## Value

This function returns a dataset with pseudo outcome. It learns
appropriate working models with Random Forest and generates pseudo
outcome using the DR-WCLS.

## Examples

``` r
sim_data = generate_dataset(N = 1000, T = 40, P = 50, sigma_residual = 1.5, sigma_randint = 1.5, main_rand = 3, rho = 0.7,
  beta_logit = c(-1, 1.6 * rep(1/50, 50)), model = ~ state1 + state2 + state3 + state4,
  beta = matrix(c(-1, 1.7, 1.5, -1.3, -1),ncol = 1),
  theta1 = 0.8)

Ht = unlist(lapply(1:50, FUN = function(X) paste0("state",X)))
St = unlist(lapply(1:25, FUN = function(X) paste0("state",X)))

pseudo_outcome_generator_rf_v2(fold = 5,ID = "id", sim_data, Ht, St, "action", "outcome",core_num = 5)
#> Error in pseudo_outcome_generator_rf_v2(fold = 5, ID = "id", sim_data,     Ht, St, "action", "outcome", core_num = 5): argument "outcome" is missing, with no default
```
