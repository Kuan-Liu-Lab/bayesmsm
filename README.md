# R package bayesmsm

[![pkgdown](https://img.shields.io/badge/pkgdown--blue)](https://Kuan-Liu-Lab.github.io/bayesmsm/)

This package provides tools for estimating causal effects using Bayesian Marginal Structural Models. It includes functions for estimating Bayesian weights using JAGS and for Bayesian non-parametric bootstrap to calculate causal effects.


# Reference

Reference paper on Bayesian marginal structural models:

-  Saarela, O., Stephens, D. A., Moodie, E. E., & Klein, M. B. (2015). On Bayesian estimation of marginal structural models. Biometrics, 71(2), 279-288.

-  Liu, K., Saarela, O., Feldman, B. M., & Pullenayegum, E. (2020). Estimation of causal effects with repeatedly measured outcomes in a Bayesian framework. Statistical methods in medical research, 29(9), 2507-2519.


# Installation

Install using `devtools` package:

```r
## install.packages(devtools) ## make sure to have devtools installed 
devtools::install_github("Kuan-Liu-Lab/bayesmsm")
library(bayesmsm)
```

# Dependency

This package depends on the following R packages:
-  `MCMCpack`
-  `doParallel`
-  `foreach`
-  `parallel`
-  `R2jags`
-  `coda`


# Examples

Here are some examples demonstrating how to use the `bayesmsm` package:

```r
# Load example data
testdata <- read.csv(system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm"))

# Calculate Bayesian weights
weights <- bayesweight(
  trtmodel.list = list(
    a_1 ~ w1 + w2 + L1_1 + L2_1,
    a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1
  ),
  data = testdata,
  n.iter = 250,
  n.burnin = 150,
  n.thin = 5,
  n.chains = 2,
  seed = 890123,
  parallel = TRUE
)

# Perform Bayesian non-parametric bootstrap
model <- bayesmsm(
  ymodel = y ~ a_1 + a_2,
  nvisit = 2,
  reference = c(rep(0, 2)),
  comparator = c(rep(1, 2)),
  family = "gaussian",
  data = testdata,
  wmean = weights,
  nboot = 1000,
  optim_method = "BFGS",
  seed = 890123,
  parallel = TRUE,
  ncore = 2
)

# View model summary
summary.bayesmsm(model)
```


# Documentation

For comprehensive documentation and examples, visit the [pkgdown website](https://Kuan-Liu-Lab.github.io/bayesmsm/).


# License

This package is licensed under the MIT License. See the LICENSE file for details.


# Citation

If you use the `bayesmsm` package in your research, please cite the reference papers listed above.


# Developers

-  Kuan Liu
-  Xiao Yan (Maintainer)
