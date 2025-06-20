% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simData.R
\name{simData}
\alias{simData}
\title{Generate synthetic longitudinal data with optional right-censoring}
\usage{
simData(
  n,
  n_visits,
  covariate_counts = rep(2, n_visits),
  amodel,
  ymodel,
  y_type = c("binary", "continuous"),
  right_censor = FALSE,
  cmodel = NULL,
  seed = NULL
)
}
\arguments{
\item{n}{Integer. Sample size.}

\item{n_visits}{Integer. Number of visits (including baseline as visit 1).}

\item{covariate_counts}{Integer vector of length `n_visits`. Number of covariates per visit (default: rep(2, n_visits)).}

\item{amodel}{List of length `n_visits`. Each element is a named numeric vector of coefficients
for the logistic model of treatment `Aj` on covariates (and `A_prev` for j > 1).}

\item{ymodel}{Named numeric vector. Coefficients for the end-of-study outcome model.
If `y_type = "binary"`, a logistic model is used; if `"continuous"`, a linear model with Gaussian noise.}

\item{y_type}{Character. One of "binary" or "continuous".}

\item{right_censor}{Logical. If TRUE, generates `Cj` using `cmodel` at each visit.}

\item{cmodel}{List of length `n_visits`. Named numeric vectors for logistic censoring models
at each visit, regressing `Cj` on covariates and current `Aj`.}

\item{seed}{Integer. Optional random seed.}
}
\value{
A `data.frame` with columns `Lk_j`, `Aj`, optional `Cj`, and `Y`.
}
\description{
This function simulates repeated measurements of normally-distributed covariates,
binary treatments, and an end-of-study outcome for longitudinal causal analyses.
When `right_censor = TRUE`, a right-censoring indicator `Cj` is generated at each visit:
if `Cj = 1`, all subsequent `L`, `A`, and `Y` values are set to `NA`.
}
