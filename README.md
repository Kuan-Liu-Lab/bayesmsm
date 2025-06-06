
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
# Load example data without right-censoring
testdata <- read.csv(system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm"))

# Load example data without right-censoring
testdata <- read.csv(system.file("extdata",
                                 "sim_causal.csv",
                                 package = "bayesmsm"))

# Calculate Bayesian weights
weights <- bayesweight(
  trtmodel.list = list(
    A1 ~ L11 + L21,
    A2 ~ L11 + L21 + L12 + L22 + A1,
    A3 ~ L11 + L21 + L12 + L22 + A1 + L13 + L23 + A2
  ),
  data = testdata,
  n.chains = 2,
  n.iter = 250,
  n.burnin = 150,
  n.thin = 5,
  seed = 890123,
  parallel = TRUE
)

# Perform Bayesian non-parametric bootstrap
model <- bayesmsm(
  ymodel = Y ~ A1 + A2 + A3,
  nvisit = 3,
  reference = c(rep(0,3)),
  comparator = c(rep(1,3)),
  treatment_effect_type = "sq",
  family = "binomial",
  data = testdata,
  wmean = weights$weights,
  nboot = 1000,
  optim_method = "BFGS",
  seed = 890123,
  parallel = TRUE,
  ncore = 2
)

# View model summary
summary.bayesmsm(model)
```

# License

This package is licensed under the MIT License. See the LICENSE file for
details.

# Citation

Please cite our software using:

    @Manual{,
      title = {bayesmsm: An R package for longitudinal causal analysis using Bayesian Marginal Structural Models},
      author = {Xiao Yan and Martin Urner and Kuan Liu},
      year = {2024},
      note = {https://github.com/Kuan-Liu-Lab/bayesmsm},
      url = {https://kuan-liu-lab.github.io/bayesmsm/},
    }

# Contact

- e-mail: <kuan.liu@utoronto.ca>, <Clarence.YXA@gmail.com>
- Please report bugs by opening an
  [issue](https://github.com/Kuan-Liu-Lab/bayesmsm/issues/new).
