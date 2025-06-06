% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesweight_cen.R
\name{bayesweight_cen}
\alias{bayesweight_cen}
\title{Bayesian Treatment Effect Weight Estimation for Censored Data}
\usage{
bayesweight_cen(
  trtmodel.list,
  cenmodel.list,
  data,
  n.chains = 2,
  n.iter = 25000,
  n.burnin = 15000,
  n.thin = 5,
  seed = NULL,
  parallel = TRUE
)
}
\arguments{
\item{trtmodel.list}{A list of formulas corresponding to each time point with the time-specific treatment variable on the left-hand side and pre-treatment covariates to be balanced on the right-hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point. Interactions and functions of covariates are allowed.}

\item{cenmodel.list}{A list of formulas for the censored data at each time point, with censoring indicators on the left-hand side and covariates on the right-hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point.}

\item{data}{A data set in the form of a data frame containing the variables in "trtmodel.list" and "cenmodel.list". This must be a wide data set with exactly one row per unit.}

\item{n.chains}{Integer specifying the number of MCMC chains to run. Set to 1 for non-parallel computation. For parallel computation, it is required to use at least 2 chains. The default is 2.}

\item{n.iter}{Integer specifying the total number of iterations for each chain (including burn-in). The default is 25000.}

\item{n.burnin}{Integer specifying the number of burn-in iterations for each chain. The default is 15000.}

\item{n.thin}{Integer specifying the thinning rate for the MCMC sampler. The default is 5.}

\item{seed}{Starting seed for the JAGS model. The default is NULL.}

\item{parallel}{Logical scalar indicating whether to run the MCMC chains in parallel. The default is TRUE.}
}
\value{
A list of the calculated weights and the JAGS model where "weights" is a vector of posterior mean weights, computed by taking the average of the weights across all MCMC iterations and `model_string` is a character of the JAGS model based on the input of "trtmodel.list".
}
\description{
This function estimates Bayesian importance sampling weights for treatment models and censoring models across multiple time points via JAGS
}
\examples{
simdat_cen <- read.csv(system.file("extdata",
                       "sim_causal_cen.csv",
                       package = "bayesmsm"))
weights_cen <- bayesweight_cen(
                trtmodel.list = list(
                A1 ~ L11 + L21,
                A2 ~ L11 + L21 + L12 + L22 + A1,
                A3 ~ L11 + L21 + L12 + L22 + A1 + L13 + L23 + A2),
                cenmodel.list = list(C ~ L11 + L21 + A1 + L12 + L22 + A2),
                data = simdat_cen,
                n.chains = 1,
                n.iter = 20,
                n.burnin = 10,
                n.thin = 1,
                seed = 890123,
                parallel = FALSE)
summary(weights_cen)
}
