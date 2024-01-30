#' Bayesian Marginal Structural Model Bootstrap Estimation
#'
#' This function performs Bayesian non-parametric bootstrap to calculate causal effects
#' for Bayesian marginal structural models.
#'
#' @param wloglik_function The weighted log-likelihood function to optimize.
#' @param testdata A data frame containing the outcome variable 'y' and treatment variables 'a_1', 'a_2', etc.
#' @param wmean A vector of weights to be used in the bootstrap procedure.
#' @param nboot The number of bootstrap iterations (default is 1000).
#' @param inits1 Initial values for the optimization algorithm (default is c(rep(0.1, length(treatment_vars)+1),4).
#' @param optim_control List of control parameters to pass to the optim() function (default is list(fnscale = -1)).
#' @param optim_method The optimization method to be used in the optim() function (default is 'BFGS').
#' @param optim_hessian Logical; argument of optim() function (default is FALSE).
#'
#' @return A list of summary measures (mean, sd, quantile)
#' @export
#'
bayesm_bootstrap <- function(wloglik_function, testdata, wmean = rep(1, 1000), nboot = 1000,
                             inits1 = c(rep(0.1, length(treatment_vars)+1),4), optim_control = list(fnscale = -1),
                             optim_method = 'BFGS', optim_hessian = FALSE) {
  bootest <- numeric(nboot)
  treatment_vars <- names(testdata)[grepl("^a_", names(testdata))]

  for (j in 1:nboot) {
    alpha <- as.numeric(rdirichlet(1, rep(1.0, length(testdata$y))))

    maxim <- optim(inits1,
                   fn = wloglik_function,
                   Y = testdata$y,
                   A = as.matrix(cbind(1, testdata[, treatment_vars])),
                   weight = alpha * wmean,
                   control = optim_control,
                   method = optim_method,
                   hessian = optim_hessian)
    bootest[j] <- sum(maxim$par[2:length(treatment_vars)+1]) # Sum the effects for all treatment variables

    if (j %% 100 == 0) {
      print(j)
    }
  }

  return(list(
    mean = mean(bootest),
    sd = sqrt(var(bootest)),
    quantile = quantile(bootest, probs = c(0.025, 0.975))
  ))
}
