% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PseudoOutcome_GBM.R
\name{pseudo_outcome_generator_gbm}
\alias{pseudo_outcome_generator_gbm}
\title{pseudo_outcome_generator_gbm}
\usage{
pseudo_outcome_generator_gbm(
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
}
\arguments{
\item{fold}{number of folds to split when do corss-fitting}

\item{ID}{the name of column where participants' ID are stored}

\item{data}{dataset name}

\item{Ht}{a vector that contains column names of control variables}

\item{St}{a vector that contains column names of moderator variables; St should be a subset of Ht}

\item{At}{column names of treatment (At)}

\item{prob}{column names of \eqn{p_t(A_t = 1|H_t)}, the experiment design treatment probability}

\item{outcome}{column names of outcome variable}

\item{core_num}{number of cores will be used for calculation}
}
\value{
This function returns a dataset with pseudo outcome. It learns appropriate working models with Gradient Boosting and
generates pseudo outcome using the DR-WCLS.
}
\description{
This function uses Gradient Boosting to conduct model training on each fold then estimates \eqn{E[Y_{t+1}|H_t, A_t], E[A_t|H_t], E[A_t|S_t]} for reserved fold. Then
the function returns a dataset with column named "yDR" for DR-WCLS pseudo outcome.
}
\examples{
pseudo_outcome_generator_gbm(fold = 5,ID = "id", data, Ht, St, "action", "outcome",core_num = 5)

}
