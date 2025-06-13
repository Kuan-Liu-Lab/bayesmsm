
<!-- README.md is generated from README.Rmd. Please edit this file -->

## bayesmsm

<!-- badges: start -->
<!-- badges: end -->

# Overview

*bayesmsm* is an R package that implements the Bayesian marginal
structural models to estimate average treatment effect for drawing
causal inference with time-varying treatment assignment and confounding
with extension to handle informative right-censoring. The Bayesian
marginal structural models is a semi-parametric approach and features a
two-step estimation process. The first step involves Bayesian parametric
estimation of the time-varying treatment assignment models and the
second step involves non-parametric Bayesian bootstrap to estimate the
average treatment effect between two distinct treatment sequences of
interest.

Reference paper on Bayesian marginal structural models:

- Saarela, O., Stephens, D. A., Moodie, E. E., & Klein, M. B. (2015). On
  Bayesian estimation of marginal structural models. Biometrics, 71(2),
  279-288.

- Liu, K., Saarela, O., Feldman, B. M., & Pullenayegum, E. (2020).
  Estimation of causal effects with repeatedly measured outcomes in a
  Bayesian framework. Statistical methods in medical research, 29(9),
  2507-2519.

# Installation

Install using `devtools` package:

``` r
## install.packages(devtools) ## make sure to have devtools installed 
devtools::install_github("Kuan-Liu-Lab/bayesmsm")
library(bayesmsm)
```

# Dependency

This package depends on the following R packages:

- `MCMCpack`
- `doParallel`
- `foreach`
- `parallel`
- `R2jags`
- `coda`

# Quick Start

Here are some examples demonstrating how to use the `bayesmsm` package:

``` r
# Simulating longitudinal causal data without right-censoring
# 1) Define coefficient lists for 2 visits
amodel <- list(
  # Visit 1: logit P(A1=1) = -0.3 + 0.4*L1_1 - 0.2*L2_1
  c("(Intercept)" = -0.3, "L1_1" = 0.4, "L2_1" = -0.2),
  # Visit 2: logit P(A2=1) = -0.1 + 0.3*L1_2 - 0.1*L2_2 + 0.5*A_prev
  c("(Intercept)" = -0.1, "L1_2" = 0.3, "L2_2" = -0.1, "A_prev" = 0.5)
)

# 2) Define binary outcome model: logistic on treatments and last covariates
ymodel <- c(
  "(Intercept)" = -0.8,
  "A1"          = 0.2,
  "A2"          = 0.4,
  "L1_2"        = 0.3,
  "L2_2"        = -0.3
)

# 3) Load package and simulate data without censoring
testdata <- simData(
  n                = 100,
  n_visits         = 2,
  covariate_counts = c(2, 2),
  amodel           = amodel,
  ymodel           = ymodel,
  y_type           = "binary",
  right_censor     = FALSE,
  seed             = 123
)


# Calculate Bayesian weights
weights <- bayesweight(
  trtmodel.list = list(
    A1 ~ L1_1 + L2_1,
    A2 ~ L1_2 + L2_2 + A1),
  data = simdat,
  n.chains = 1,
  n.iter = 200,
  n.burnin = 100,
  n.thin = 1,
  seed = 890123,
  parallel = FALSE)

# Perform Bayesian non-parametric bootstrap
model <- bayesmsm(ymodel = Y ~ A1 + A2,
  nvisit = 2,
  reference = c(rep(0,2)),
  comparator = c(rep(1,2)),
  family = "binomial",
  data = simdat,
  wmean = weights$weights,
  nboot = 50,
  optim_method = "BFGS",
  parallel = TRUE,
  seed = 890123,
  ncore = 2)

# View model summary
summary_bayesmsm(model)
```

# License

This package is licensed under the MIT License. See the LICENSE file for
details.

# Citation

Please cite our software using:

    @Manual{,
      title = {bayesmsm: An R package for longitudinal causal analysis using Bayesian Marginal Structural Models},
      author = {Xiao Yan and Martin Urner and Kuan Liu},
      year = {2025},
      note = {https://github.com/Kuan-Liu-Lab/bayesmsm},
      url = {https://kuan-liu-lab.github.io/bayesmsm/},
    }

# Contact

- e-mail: <kuan.liu@utoronto.ca>, <Clarence.YXA@gmail.com>
- Please report bugs by opening an
  [issue](https://github.com/Kuan-Liu-Lab/bayesmsm/issues/new).
