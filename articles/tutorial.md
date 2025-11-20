# MRTpostInfLASSO

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

This vignette introduces the **MRTpostInfLASSO** package. Its core
function, **DR_WCLS_LASSO**, allows users to perform variable selection,
estimate time-varying causal effects and make valid inferences
conditional on the selected variables.

Individual-level data of an MRT can be summarized as
$\left\{ O_{1},A_{1},O_{2},A_{2},\cdots,O_{T},A_{T},O_{T + 1} \right\}$
where $T$ is the total decision times, $O_{t}$ is the information
collected between $t - 1$ and $t$, and $A_{t}$ is the treatment provided
at time $t$. Here we consider treatment
$A_{t} \in \left\{ 0,1 \right\}$. Treatment options are intended to
influence a proximal outcome $Y_{t + 1} \in O_{t + 1}$.

Denote history
$H_{t} = \left\{ O_{1},A_{1},O_{2},A_{2},\cdots,A_{t - 1},O_{t} \right\}$
and randomized probabilities
$\mathbf{p} = \left\{ p_{t}\left( A_{t} \mid H_{t} \right) \right\}_{t = 1}^{T}$.
The DR-WCLS criterion is given by

\$\$ \mathbb{P}\_n \Bigl\[ \sum\_{t=1}^{T} \tilde{\sigma}^2_t(S_t)\\
\Bigl( \frac{W_t(A_t-\tilde{p}\_t(1 \mid S_t))
(Y\_{t+1}-g_t(H_t,A_t))}{\tilde{\sigma}^2_t(S_t)}\\
+\beta(t;H_t)-f_t(S_t)^T \beta\Bigr)\\ f_t(S_t) \Bigr\] = 0 \$\$

where
$\beta\left( t;H_{t} \right):=g_{t}\left( H_{t},1 \right) - g_{t}\left( H_{t},0 \right)$
is the causal excursion effect under the fully observed history $H_{t}$,
and
${\widetilde{\sigma}}_{t}^{2}\left( S_{t} \right):={\widetilde{p}}_{t}\left( 1 \mid S_{t} \right)\left( 1 - {\widetilde{p}}_{t}\left( 1 \mid S_{t} \right) \right)$.

The ${\widehat{\beta}}_{n}^{(DR)}$ is a consistent estimator of the true
$\beta$ if either the randomization probability
$p_{t}\left( A_{t} \mid H_{t} \right)$ or the conditional expectation
$g_{t}\left( H_{t},A_{t} \right)$ is correctly specified.

The **DR_WCLS** algorithm is as follows:

Step I: Randomly split the $n$ individuals into $K$ equal folds
$\left\{ I_{k} \right\}_{k = 1}^{K}$ assuming $n$ is a multiple of $K$.
Let $I_{k}^{\complement}$ denote the complement of fold k.

Step II: For each fold $k$, use data from $I_{k}^{\complement}$ to
estimate the nuisance functions
${\widehat{g}}_{t}^{(k)}\left( H_{t},A_{t} \right)$,${\widehat{p}}_{t}^{(k)}\left( 1 \mid H_{t} \right)$,
${\widehat{\widetilde{p}}}_{t}^{(k)}\left( 1 \mid S_{t} \right)$, and
compute the weight
${\widehat{W}}_{t}^{(k)} = {\widehat{\widetilde{p}}}_{t}^{(k)}\left( 1 \mid S_{t} \right)/{\widehat{p}}_{t}^{(k)}\left( 1 \mid H_{t} \right)$.

Step III: For each $j \in I_{k}$ and time $t$, construct the
pseudo-outcome ${\widetilde{Y}}_{t + 1}^{(DR)}$ as follows, then regress
it on $f_{t}\left( S_{t} \right)^{T}\beta$ using weights
${\widetilde{p}}_{t}^{(k)}\left( 1 \mid S_{t} \right)\left( 1 - {\widetilde{p}}_{t}^{(k)}\left( 1 \mid S_{t} \right) \right)$.

\$\$ \tilde{Y}^{(DR)}\_{t+1,j} :=
\frac{\hat{W}\_{t,j}^{(k)}(A\_{t,j}-\hat{\tilde{p}}\_t^{(k)}(1 \mid
S\_{t,j}))
(Y\_{t+1,j}-\hat{g}\_t^{(k)}(H\_{t,j},A\_{t,j}))}{\hat{\tilde{p}}\_t^{(k)}(1
\mid S\_{t,j})(1-\hat{\tilde{p}}\_t^{(k)}(1 \mid S\_{t,j}))} \\ + \Bigl(
\hat{g}\_t^{(k)}(H\_{t,j},1) - \hat{g}\_t^{(k)}(H\_{t,j},0) \Bigr) \$\$

To conduct variable selection, **DR_WCLS_LASSO** solves the problem

\$\$ \min\_{\beta} \frac{1}{n} \sum\_{i=1}^{n}\sum\_{t=1}^{T} \Bigl\[
\hat{\tilde{p}}^{(k)}\_{t}(1 \mid S_t)\\
\bigl(1-\hat{\tilde{p}}^{(k)}\_{t}(1 \mid S_t)\bigr)\\
\bigl(\tilde{Y}^{(DR)}\_{t+1,i}-f_t(S_t)^{\top}\beta\bigr)^2 \Bigr\] +
\lambda \lVert \beta \rVert\_{1} - w^{\top}\beta, \$\$

where $\lambda$ is the LASSO regularization parameter, $\omega$ is the
noise vector.

After the variable selection procedure, we conduct post-selection
inference using DR-WCLS conditional on the selected variables. This
provides estimates for the selected variables along with their
corresponding confidence intervals.

\$\$ \min\_{\beta} \frac{1}{n} \sum\_{i=1}^{n}\sum\_{t=1}^{T} \Bigl\[
\hat{\tilde{p}}^{(k)}\_{t}(1 \mid S_t)\\
\bigl(1-\hat{\tilde{p}}^{(k)}\_{t}(1 \mid S_t)\bigr)\\
\bigl(\tilde{Y}^{(DR)}\_{t+1,i}-f_t(S_t)^{\top}\beta_E\bigr)^2 \Bigr\]
\$\$

## Installation

The package can be installed from our GitHub repository.

``` r
# Install From GitHub
# remotes::install_github("WHD-Lab/DR_WCLS_LASSO")
```

### Loading the Package

``` r
library(MRTpostInfLASSO)
```

## Real Data Example

### HeartSteps

We illustrate the functions in the MRTpostInfLASSO package using the
`HeartSteps` dataset from the MRTAnalysis package. We demonstrate how to
generate the pseudo-outcome, perform variable selection, and conduct
valid post-selection inference.

HeartSteps is a mobile health intervention designed to encourage
physical activity by delivering tailored suggestions. The dataset comes
from a 6-week micro-randomized trial involving 37 participants.
Participants were randomized at 5 decision points per day and received
2-5 notificaitons daily. The dataset contains 7,770 records, with 5
observations per day for each participants. Each record includes the
decision point, whether a notification was sent, participant
availability for walking, and 30-minute step counts before and after the
decision point.

To begin, we load the `HeartSteps` data using the following code. A
summary of `data_mimicHeartSteps` is as follows:

``` r
# Load HeartSteps Data
library(MRTAnalysis)
data(data_mimicHeartSteps)
head(data_mimicHeartSteps)
#>   userid decision_point day_in_study logstep_30min logstep_30min_lag1
#> 1      1              1            0     2.3902011          0.0000000
#> 2      1              2            0    -0.6931472          2.3902011
#> 3      1              3            0     2.4646823         -0.6931472
#> 4      1              4            0     0.1206936          2.4646823
#> 5      1              5            0     0.8322060          0.1206936
#> 6      1              6            1     1.8450452          0.8322060
#>   logstep_pre30min is_at_home_or_work intervention rand_prob avail
#> 1       -0.6931472                  1            0       0.6     0
#> 2        2.1962380                  1            0       0.6     1
#> 3        4.5894007                  1            1       0.6     1
#> 4        3.1791124                  1            1       0.6     1
#> 5        3.2945170                  0            0       0.6     0
#> 6        4.6658254                  1            0       0.6     0
```

We first specify the variable names for the participant ID,
history$H_{t}$, moderator$S_{t}$, treatment$A_{t}$, proximal
outcome$Y_{t}$, and randomization probability$p_{t}$.

``` r
set.seed(100)
ID = 'userid'
Ht = c('logstep_30min_lag1','logstep_pre30min','is_at_home_or_work', 'day_in_study')
St = c('logstep_30min_lag1','logstep_pre30min','is_at_home_or_work', 'day_in_study')
At = 'intervention'
outcome = 'logstep_30min'
prob = 'rand_prob'
```

#### Generating Pseudo-outcome

To generate the pseudo-outcome, we provide three methods for estimating
the nuisance functions: LASSO, random forest and gradient boosting.

\$\$ \tilde{Y}^{(DR)}\_{t+1,j} :=
\frac{\hat{W}\_{t,j}^{(k)}(A\_{t,j}-\hat{\tilde{p}}\_t^{(k)}(1 \mid
S\_{t,j}))
(Y\_{t+1,j}-\hat{g}\_t^{(k)}(H\_{t,j},A\_{t,j}))}{\hat{\tilde{p}}\_t^{(k)}(1
\mid S\_{t,j})(1-\hat{\tilde{p}}\_t^{(k)}(1 \mid S\_{t,j}))} \\ + \Bigl(
\hat{g}\_t^{(k)}(H\_{t,j},1) - \hat{g}\_t^{(k)}(H\_{t,j},0) \Bigr) \$\$

We illustrate their use use via the functions
`pseudo_outcome_generator_CVlasso`, `pseudo_outcome_generator_rf_v2`,
and `pseudo_outcome_generator_gbm`.

``` r
pseudo_outcome_CVlasso = pseudo_outcome_generator_CVlasso(fold = 5,ID = ID,
                                                     data = data_mimicHeartSteps, 
                                                     Ht = Ht, St = St, At = At, 
                                                     prob = prob, outcome = outcome,
                                                     core_num = 1)
#> Loading required package: parallel

pseudo_outcome_RF = pseudo_outcome_generator_rf_v2(fold = 5,ID = ID,
                                                   data = data_mimicHeartSteps, 
                                                   Ht = Ht, St = St, At = At, 
                                                   prob = prob, outcome = outcome,
                                                   core_num = 1)


pseudo_outcome_GBM = pseudo_outcome_generator_gbm(fold = 5,ID = ID,
                                                  data = data_mimicHeartSteps, 
                                                  Ht = Ht, St = St, At = At, 
                                                  prob = prob, outcome = outcome,
                                                  core_num = 1)
```

#### Variable Selection

To perform variable selection, **DR_WCLS_LASSO** solves the problem

\$\$ \min\_{\beta} \frac{1}{n} \sum\_{i=1}^{n}\sum\_{t=1}^{T} \Bigl\[
\hat{\tilde{p}}^{(k)}\_{t}(1 \mid S_t)\\
\bigl(1-\hat{\tilde{p}}^{(k)}\_{t}(1 \mid S_t)\bigr)\\
\bigl(\tilde{Y}^{(DR)}\_{t+1,i}-f_t(S_t)^{\top}\beta\bigr)^2 \Bigr\] +
\lambda \lVert \beta \rVert\_{1} - w^{\top}\beta \$\$

Variable selection is performed using the
`variable_selection_PY_penal_int` or `FISTA_backtracking` function. We
first define `my_formula`, specifying the outcome and candidate
variables. Below, we demonstrate the use of both functions.

Python version

``` r
my_formula = as.formula(paste("yDR ~", paste(c("logstep_30min_lag1", "logstep_pre30min",
                                      "is_at_home_or_work", "day_in_study"),
                                       collapse = " + ")))
# 
# set.seed(100)
# var_selection_python = variable_selection_PY_penal_int(data = pseudo_outcome_CVlasso,ID,
#                                                        my_formula,
#                                                        lam = NULL, noise_scale = NULL,
#                                                        splitrat = 0.8,
#                                                        # virtualenv_path = "fakepath/wcls",
#                                                        ridge_term = 0, beta = NULL)
# 
# cat('The selected variable list:',var_selection_python$E)
```

R version

``` r
set.seed(100)
var_selection_R = FISTA_backtracking(data = pseudo_outcome_CVlasso, ID, my_formula,
                  lam = NULL, noise_scale = NULL,splitrat = 0.8,
                  max_ite = 10^(5), tol = 10^(-4), beta = NULL)
var_selection_R$E
#> [1] "(Intercept)"  "day_in_study"

cat('The selected variable list:',var_selection_R$E)
#> The selected variable list: (Intercept) day_in_study
```

#### Post-selection Inference

After variable selection, we conduct post-selection inference using
DR-WCLS, conditioning on the selected variables. This yields GEE
coefficient estimates and corresponding confidence intervals.

\$\$ \min\_{\beta} \frac{1}{n} \sum\_{i=1}^{n}\sum\_{t=1}^{T} \Bigl\[
\hat{\tilde{p}}^{(k)}\_{t}(1 \mid S_t)\\
\bigl(1-\hat{\tilde{p}}^{(k)}\_{t}(1 \mid S_t)\bigr)\\
\bigl(\tilde{Y}^{(DR)}\_{t+1,i}-f_t(S_t)^{\top}\beta_E\bigr)^2 \Bigr\]
\$\$

`DR_WCLS_LASSO` is used for post-selection inference. The input of the
function are the data for inference, number of folds, and variable names
for ID, time, $H_{t}$, $S_{t}$, $A_{t}$, $Y_{t}$ and randomization
probability. `method_pseu` specifies the method used for pseudo-outcome
generation (e.g., “CVLASSO”, “RandomForest”, “GradientBoosting”).
`varSelect_program` indicates the variable selection programming
language.

``` r
# set.seed(123)

# reticulate::py_run_string("
# import random
# import numpy as np
# 
# random.seed(123)
# np.random.seed(123)
# ")

# UI_return_python = DR_WCLS_LASSO(data = data_mimicHeartSteps,
#                           fold = 5, ID = ID,
#                           time = "decision_point",
#                           Ht = Ht, St = St, At = At,
#                           prob = prob, outcome = outcome,
#                           # virtualenv_path = 'fakepath/wcls',
#                           method_pseu = "CVLASSO",
#                           varSelect_program = "Python",
#                           standardize_x = F, standardize_y = F)
# 
# UI_return_python
```

``` r
set.seed(100)

UI_return_R = DR_WCLS_LASSO(data = data_mimicHeartSteps,
                          fold = 5, ID = ID,
                          time = "decision_point",
                          Ht = Ht, St = St, At = At,
                          prob = prob, outcome = outcome,
                          method_pseu = "CVLASSO", 
                          varSelect_program = "R",
                          standardize_x = F, standardize_y = F)
#> [1] "remove 0 lines of data due to NA produced for yDR"
#> [1] "The current lambda value is: 130.868168720783"
#> [1] "select predictors: (Intercept)"  "select predictors: day_in_study"
#> [1] FALSE
#> [1] "logstep_30min_lag1" "logstep_pre30min"   "is_at_home_or_work"
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
UI_return_R
#>              E     GEE_est      lowCI     upperCI   prop_low   prop_up
#> 1  (Intercept)  0.56714755  0.3711479  0.72711034 0.05046974 0.9499734
#> 2 day_in_study -0.02083236 -0.0286465 -0.01228228 0.04950436 0.9503883
#>         pvalue
#> 1 1.279054e-06
#> 2 1.674468e-04
```

#### A Comparison of Using Randomized LASSO and Weighted Centered Least Squares

Select variables using randomized LASSO

``` r
library(selectiveInference)
#> Loading required package: glmnet
#> Loading required package: Matrix
#> Loaded glmnet 4.1-10
#> Loading required package: intervals
#> 
#> Attaching package: 'intervals'
#> The following object is masked from 'package:Matrix':
#> 
#>     expand
#> Loading required package: survival
#> Loading required package: adaptMCMC
#> Loading required package: coda
#> Loading required package: MASS
#> 
#> Attaching package: 'MASS'
#> The following object is masked from 'package:dplyr':
#> 
#>     select
res_randomizedLASSO = randomizedLasso(X = as.matrix(data_mimicHeartSteps[,St]), 
                y = as.matrix(data_mimicHeartSteps[,outcome]), 
                lam = 100000, 
                family="gaussian",
                noise_scale=NULL, 
                ridge_term=NULL, 
                max_iter=100,       
                kkt_tol=1.e-4,      
                parameter_tol=1.e-8,
                objective_tol=1.e-8,
                objective_stop=FALSE,
                kkt_stop=TRUE,
                parameter_stop=TRUE)

St[res_randomizedLASSO$active_set]
#> [1] "day_in_study"
```

Make inference using weighted centered least squares

``` r

data_mimicHeartSteps$intervention = as.numeric(data_mimicHeartSteps$intervention)
data_mimicHeartSteps$logstep_30min = as.numeric(data_mimicHeartSteps$logstep_30min)
data_mimicHeartSteps$logstep_30min_lag1 = as.numeric(data_mimicHeartSteps$logstep_30min_lag1)
data_mimicHeartSteps$logstep_pre30min = as.numeric(data_mimicHeartSteps$logstep_pre30min)
data_mimicHeartSteps$is_at_home_or_work = as.numeric(data_mimicHeartSteps$is_at_home_or_work)
data_mimicHeartSteps$day_in_study = as.numeric(data_mimicHeartSteps$day_in_study)

wcls_fit = wcls(
  data = data_mimicHeartSteps,
  id = 'userid', 
  outcome = 'logstep_30min', 
  treatment = 'intervention', 
  rand_prob = 0.5, 
  moderator_formula= ~ day_in_study,
  control_formula = ~logstep_30min_lag1 + logstep_pre30min + day_in_study,
  availability = NULL,
  numerator_prob = NULL,
  verbose = TRUE
)
#> availability = NULL: defaulting availability to always available.
#> Constant randomization probability 0.5 is used.
#> Constant numerator probability 0.5 is used.

wcls_res = summary(wcls_fit)
wcls_res$causal_excursion_effect
#>                Estimate     95% LCL    95% UCL      StdErr Hotelling df1 df2
#> (Intercept)   0.5722731  0.37187918  0.7726669 0.098255727  33.92273   1  31
#> day_in_study -0.0209880 -0.03004471 -0.0119313 0.004440621  22.33854   1  31
#>                   p-value
#> (Intercept)  2.024579e-06
#> day_in_study 4.698739e-05

# UI_return_python
UI_return_R
#>              E     GEE_est      lowCI     upperCI   prop_low   prop_up
#> 1  (Intercept)  0.56714755  0.3711479  0.72711034 0.05046974 0.9499734
#> 2 day_in_study -0.02083236 -0.0286465 -0.01228228 0.04950436 0.9503883
#>         pvalue
#> 1 1.279054e-06
#> 2 1.674468e-04
```

![](tutorial_files/figure-html/unnamed-chunk-13-1.png)![](tutorial_files/figure-html/unnamed-chunk-13-2.png)

#### Arguments in DR_WCLS_LASSO

The pseudo outcome generation method can be specified using the
`method_pseu` argument. Currently, three options are supported: CVLASSO,
RandomForest, and GradientBoosting. The default is CVLASSO.

``` r
set.seed(100)

UI_return_method_pseu = DR_WCLS_LASSO(data = data_mimicHeartSteps,
                                 fold = 5, ID = ID,
                                 time = "decision_point",
                                 Ht = Ht, St = St, At = At,
                                 prob = prob, outcome = outcome,
                                 # virtualenv_path = 'fakepath/wcls',
                                 method_pseu = "CVLASSO",
                                 varSelect_program = "R",
                                 standardize_x = F, standardize_y = F)
#> [1] "remove 0 lines of data due to NA produced for yDR"
#> [1] "The current lambda value is: 130.888825329133"
#> [1] "select predictors: (Intercept)"  "select predictors: day_in_study"
#> [1] FALSE
#> [1] "logstep_30min_lag1" "logstep_pre30min"   "is_at_home_or_work"

UI_return_method_pseu
#>              E     GEE_est       lowCI    upperCI   prop_low   prop_up
#> 1  (Intercept)  0.56830179  0.37300973  0.7275392 0.05095176 0.9494623
#> 2 day_in_study -0.02088079 -0.02867089 -0.0123520 0.04972344 0.9500787
#>         pvalue
#> 1 1.168007e-06
#> 2 1.578669e-04
```

The LASSO penalty can be adjusted by setting ‘lam’ in the
`DR_WCLS_LASSO` function.

``` r
set.seed(100)

UI_return_lambda = DR_WCLS_LASSO(data = data_mimicHeartSteps,
                                 fold = 5, ID = ID,
                                 time = "decision_point",
                                 Ht = Ht, St = St, At = At,
                                 prob = prob, outcome = outcome,
                                 # virtualenv_path = 'fakepath/wcls',
                                 method_pseu = "CVLASSO",
                                 varSelect_program = "R",
                                 lam = 100,
                                 standardize_x = F, standardize_y = F)
#> [1] "remove 0 lines of data due to NA produced for yDR"
#> [1] "The current lambda value is: 100"
#> [1] "select predictors: (Intercept)"  "select predictors: day_in_study"
#> [1] FALSE
#> [1] "logstep_30min_lag1" "logstep_pre30min"   "is_at_home_or_work"

UI_return_lambda
#>              E     GEE_est       lowCI     upperCI   prop_low   prop_up
#> 1  (Intercept)  0.56744699  0.37069497  0.72732422 0.05022313 0.9502763
#> 2 day_in_study -0.02084919 -0.02868708 -0.01233107 0.04958886 0.9509224
#>         pvalue
#> 1 1.021045e-06
#> 2 1.333082e-04
```

The data split rate in Step 1 of the DR_WCLS algorithm can be set using
‘splitrat’ in the `DR_WCLS_LASSO` function.

``` r
set.seed(100)
UI_return_splitrat = DR_WCLS_LASSO(data = data_mimicHeartSteps, 
                                   fold = 5, ID = ID, 
                                   time = "decision_point", 
                                   Ht = Ht, St = St, At = At, 
                                   prob = prob, outcome = outcome,
                                   method_pseu = "CVLASSO", 
                                   varSelect_program = "R",
                                   splitrat = 0.8,
                                   standardize_x = F, standardize_y = F)
#> [1] "remove 0 lines of data due to NA produced for yDR"
#> [1] "The current lambda value is: 130.874387390531"
#> [1] "select predictors: (Intercept)"  "select predictors: day_in_study"
#> [1] FALSE
#> [1] "logstep_30min_lag1" "logstep_pre30min"   "is_at_home_or_work"
UI_return_splitrat
#>              E     GEE_est       lowCI     upperCI   prop_low   prop_up
#> 1  (Intercept)  0.56759596  0.37129034  0.72776978 0.05034085 0.9500993
#> 2 day_in_study -0.02085567 -0.02868189 -0.01228783 0.04927559 0.9505832
#>         pvalue
#> 1 1.268071e-06
#> 2 1.651741e-04
```

#### Analysis Using Manually Created Interaction Terms

Interaction terms can be manually created and included in the $H_{t}$
and $S_{t}$ variable lists. In this example, we create an indicator for
whether the time in the study is over 14 days and include its
interactions with `logstep_30min_lag1`, `logstep_pre30min`, and
`is_at_home_or_work`.

``` r
data_mimicHeartSteps$timeover14 = as.numeric(data_mimicHeartSteps$decision_point>14)
data_mimicHeartSteps$int_lag1_timeover14 = data_mimicHeartSteps$logstep_30min_lag1 * data_mimicHeartSteps$timeover14

data_mimicHeartSteps$int_pre30_timeover14 = data_mimicHeartSteps$logstep_pre30min * data_mimicHeartSteps$timeover14

data_mimicHeartSteps$int_home_timeover14 = data_mimicHeartSteps$is_at_home_or_work * data_mimicHeartSteps$timeover14
```

We add the interaction terms in the $H_{t}$ and $S_{t}$ variable lists
and rerun the procedure.

``` r
set.seed(200)
ID = 'userid'
Ht_int = c('logstep_30min_lag1','logstep_pre30min','is_at_home_or_work', 'timeover14',
       'int_lag1_timeover14', 'int_pre30_timeover14', 'int_home_timeover14','day_in_study')
St_int = c('logstep_30min_lag1','logstep_pre30min','is_at_home_or_work', 'timeover14',
       'int_lag1_timeover14', 'int_pre30_timeover14', 'int_home_timeover14','day_in_study')
At = 'intervention'
outcome = 'logstep_30min'
prob = 'rand_prob'


# UI_return_int_python = DR_WCLS_LASSO(data = data_mimicHeartSteps, 
#                           fold = 5, ID = ID, 
#                           time = "decision_point", Ht = Ht_int, St = St_int, 
#                           At = At, prob = prob, outcome = outcome,
#                           # virtualenv_path = 'fakepath/wcls',
#                           method_pseu = "CVLASSO", 
#                           varSelect_program = "Python",
#                           standardize_x = F, standardize_y = F)
# 
# UI_return_int_python
```

``` r
UI_return_int_R = DR_WCLS_LASSO(data = data_mimicHeartSteps, 
                          fold = 5, ID = ID, 
                          time = "decision_point", Ht = Ht_int, St = St_int, 
                          At = At, prob = prob, outcome = outcome,
                          # virtualenv_path = 'fakepath/wcls',
                          method_pseu = "CVLASSO", 
                          varSelect_program = "R",
                          standardize_x = F, standardize_y = F)
#> [1] "remove 0 lines of data due to NA produced for yDR"
#> [1] "The current lambda value is: 152.861099902094"
#> [1] "select predictors: (Intercept)"       
#> [2] "select predictors: logstep_30min_lag1"
#> [3] "select predictors: logstep_pre30min"  
#> [4] "select predictors: day_in_study"      
#> [1] FALSE
#> [1] "is_at_home_or_work"   "timeover14"           "int_lag1_timeover14" 
#> [4] "int_pre30_timeover14" "int_home_timeover14"

UI_return_int_R
#>                    E      GEE_est       lowCI      upperCI   prop_low   prop_up
#> 1        (Intercept)  0.561683648 -0.66248660  0.595717013 0.05019632 0.9496353
#> 2 logstep_30min_lag1 -0.001813508 -0.03943796  0.141708464 0.04902147 0.9507998
#> 3   logstep_pre30min -0.001274280 -0.01350841  0.159601041 0.04965494 0.9502152
#> 4       day_in_study -0.020247918 -0.02870156 -0.004984289 0.05049238 0.9501175
#>       pvalue
#> 1 0.94301393
#> 2 0.55333936
#> 3 0.19788783
#> 4 0.03907502
```

### Intern Health Study

IHS is a micro-randomized trial involving 859 medical interns. The aim
of the study is to investigate how mHealth interventions affect
participants’ weekly mood, physical activity, and sleep. Each day,
participants had a 0.5 probability of receiving a tailored message. The
dataset includes a range of variables collected from wearable devices.

``` r
df_IHS = read.csv('~/Desktop/25Winter/IHS Design/Hyperparameter/dfEventsThreeMonthsTimezones.csv')
df_IHS$prob = rep(0.5, length(df_IHS$PARTICIPANTIDENTIFIER))
df_IHS$less27 = (df_IHS$Age < 27)

ID = 'PARTICIPANTIDENTIFIER'

# Ht = c('StepsPastDay', 'is_weekend', 'maxHRPast24Hours', 'HoursSinceLastMood')
# St = c('StepsPastDay', 'is_weekend', 'maxHRPast24Hours',  'HoursSinceLastMood')
Ht = c('StepsPastDay', 'is_weekend', 'maxHRPast24Hours', 'HoursSinceLastMood','less27', "Sex", "has_child","StepsPastWeek",'PHQ10above0','exercisedPast24Hours')
St = c('StepsPastDay', 'is_weekend', 'maxHRPast24Hours',  'HoursSinceLastMood', 'less27', "Sex", "has_child",'StepsPastWeek','PHQ10above0','exercisedPast24Hours')
At = 'sent'
outcome = 'LogStepsReward'
prob = 'prob'

set.seed(99)

df_IHS_cleaned = na.omit(df_IHS)

# UI_return_IHS_python = DR_WCLS_LASSO(data = df_IHS_cleaned,
#                           fold = 5, ID = ID,
#                           time = "time", Ht = Ht, St = St, At = At,
#                           prob = prob, outcome = outcome,
#                           # virtualenv_path = 'fakepath/wcls',
#                           method_pseu = "CVLASSO", lam = 55,
#                           noise_scale = NULL, splitrat = 0.8,
#                           level = 0.9, core_num=3, CI_algorithm = 'lapply',
#                           max_iterate = 10^{6}, max_tol = 10^{-3}, varSelect_program = "Python",
#                           standardize_x =  T, standardize_y = T)
# 
# UI_return_IHS_python
```

``` r
set.seed(100)
UI_return_IHS_R = DR_WCLS_LASSO(data = df_IHS_cleaned,
                          fold = 5, ID = ID,
                          time = "time", Ht = Ht, St = St, At = At,
                          prob = prob, outcome = outcome,
                          # virtualenv_path = 'fakepath/wcls',
                          method_pseu = "CVLASSO", lam = 30,
                          noise_scale = NULL, splitrat = 0.8,
                          level = 0.9, core_num=3, CI_algorithm = 'lapply',
                          max_iterate = 10^{6}, max_tol = 10^{-3}, varSelect_program = "R",
                          standardize_x =  T, standardize_y = T)

UI_return_IHS_R
```

``` r
# df_IHS$sent = as.numeric(df_IHS$sent)
# df_IHS$is_weekend = as.numeric(df_IHS$is_weekend)
# library(MRTAnalysis)
# wcls(
#   data = df_IHS,
#   id = 'PARTICIPANTIDENTIFIER',
#   outcome = 'Steps24HoursAfter',
#   treatment = 'sent',
#   rand_prob = 'prob',
#   moderator_formula= ~1,
#   control_formula = ~StepsPastDay + is_weekend + maxHRPast24Hours + HoursSinceLastMood,
#   availability = NULL,
#   numerator_prob = NULL,
#   verbose = TRUE
# )
```

``` r
# formula_str = paste(outcome, "~", paste(Ht, collapse = " + "))
# lm_model = lm(as.formula(formula_str), data = df_IHS)
# car::vif(lm_model)
```
