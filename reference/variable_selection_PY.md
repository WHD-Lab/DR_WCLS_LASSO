# variable_selection_PY

This function interacts with Python virtual environment to conduct
randomized LASSO variable selection. The intercept is not penalized.

## Usage

``` r
variable_selection_PY(
  data,
  ID,
  moderator_formula,
  lam = NULL,
  noise_scale = NULL,
  splitrat = 0.8,
  virtualenv_path = "",
  beta = NULL
)
```

## Arguments

- data:

  the output of pseudo_outcomecal function

- ID:

  the name of column where participants' ID are stored

- moderator_formula:

  determines the formula for the f(St)T\*beta function

- lam:

  the value of penalty term of randomized LASSO. If it is not provided,
  the default value will be used. Default is \\\sqrt{2n\* logp} \rho
  sd(y)\\ where \\\rho\\ is the split rate and \\n\\ is the number of
  rows.

- noise_scale:

  Standard deviation of Gaussian noise added to objective. Default is
  \\\sqrt{\frac{1-\rho}{\rho}\\ n}\\\mathrm{sd}(y)\\ where \\\rho\\ is
  the split rate and \\n\\ is the number of rows. The random noises,
  \\\omega\\, are iid drawn from normal distribution with standard
  deviation noise_scale

- splitrat:

  this value is corresponding to the data splitting rate. Details can
  read "Exact Selective Inference with Randomization" page 15 equation
  (10). This value will be used only when user doesn't provide the
  `noise_scale` or `lam`.

- virtualenv_path:

  Python virtual environment path

- beta:

  True coefficients (for simulation use only).

## Value

- formula:

  the raw model before selection

- E:

  Selected variables

- NE:

  Not selected variables

- n:

  the number of participants in the study

- perturb:

  Random noise added for each variables. This value is scaled by -2 for
  later inference purpose.

- ori_lam:

  Regularization parameter (lambda) used in the penalization.

- lam:

  Scale orginal regularization parameter by -2 for later inference
  purpose.

- Z:

  Subgradient of unselected variables.

- OMEGA:

  The scaled variance for the added random noise. This value is scaled
  by 4.

- sign_soln:

  Signs of the estimated coefficients for selected variables

- soln:

  Estimated coefficients of selected variables from the penalized
  randomized regression.

- postbeta:

  True beta projected on space spinned by post-selection predictors

- nonzero:

  Boolean vector indicating which variables were selected.

## Examples

``` r
sim_data = generate_dataset(N = 1000, T = 40, P = 50, sigma_residual = 1.5, sigma_randint = 1.5, main_rand = 3, rho = 0.7,
  beta_logit = c(-1, 1.6 * rep(1/50, 50)), model = ~ state1 + state2 + state3 + state4,
  beta = matrix(c(-1, 1.7, 1.5, -1.3, -1),ncol = 1),
  theta1 = 0.8)

Ht = unlist(lapply(1:50, FUN = function(X) paste0("state",X)))
St = unlist(lapply(1:25, FUN = function(X) paste0("state",X)))

pseudo_outcome_CVlasso = pseudo_outcome_generator_CVlasso(fold = 5,ID = "id", sim_data, Ht,
  St, "action", "outcome",core_num = 5)
#> Error in pseudo_outcome_generator_CVlasso(fold = 5, ID = "id", sim_data,     Ht, St, "action", "outcome", core_num = 5): argument "outcome" is missing, with no default

my_formula = as.formula(paste("yDR ~ ", paste(St, collapse = " + ")))

var_selection_python = variable_selection_PY(data = pseudo_outcome_CVlasso, ID = "id",
  my_formula, lam = NULL, noise_scale = NULL, splitrat = 0.7,
  virtualenv_path = "fakepath/wcls", beta = matrix(c(-1, 1.7, 1.5, -1.3, -1, rep(0,21)), ncol = 1))
#> Error: object 'pseudo_outcome_CVlasso' not found
```
