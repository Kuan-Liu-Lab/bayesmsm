#' Bayesian Marginal Structural Model Bootstrap Estimation
#'
#' This function performs Bayesian non-parametric bootstrap to estimate causal
#' effects in Bayesian marginal structural models. It supports both continuous
#' (gaussian) and binary (binomial) outcome variables.
#'
#' @param ymodel A formula representing the outcome model with interactions.
#' @param nvisit Number of visits or time points.
#' @param reference Reference intervention across all visits. Default is a vector of all 0's with length nvisit (i.e. never treated).
#' @param comparator Comparison intervention across all visits. Default is a vector of all 1's with length nvisit (i.e. always treated).
#' @param family Outcome distribution family; "gaussian" (default) for continuous outcomes or "binomial" for binary outcomes.
#' @param data The dataset.
#' @param wmean Vector of treatment assignment weights. Default is rep(1, 1000).
#' @param nboot Number of bootstrap iterations.
#' @param optim_method Optimization method to be used. Default is 'BFGS'.
#' @param seed A seed to ensure reproducibility.
#' @param parallel Whether parallel computation should be used. Default is TRUE.
#' @param ncore Number of cores to use for parallel computation. Default is 4.
#'
#' @return It returns an object of class `bayesmsm` that contains the information about the data, model, etc.
#'
#' An object of class `bayesmsm` is a list containing at least the following components:
#' * `mean`, the mean of the bootstrap estimates
#' * `sd`, the standard deviation of the bootstrap estimates
#' * `quantile`, the 95% quantiles of the bootstrap estimates
#' * `bootdata`, a data frame of bootstrapped estimates
#' * `reference`, the reference intervention level
#' * `comparator`, the camparison intervention level
#'
#' @importFrom foreach "%dopar%"
#' @import doParallel
#' @import parallel
#' @import MCMCpack
#' @importFrom MCMCpack rdirichlet
#'
#' @export
#'
#' @examples
#'
#' # Continuous outcome
#' testdata <- read.csv(system.file("extdata",
#'                                  "continuous_outcome_data.csv",
#'                                  package = "bayesmsm"))
#' model <- bayesmsm(ymodel = y ~ a_1+a_2,
#'                            nvisit = 2,
#'                            reference = c(rep(0,2)),
#'                            comparator = c(rep(1,2)),
#'                            family = "gaussian",
#'                            data = testdata,
#'                            wmean = rep(1, 1000),
#'                            nboot = 100,
#'                            optim_method = "BFGS",
#'                            seed = 890123,
#'                            parallel = TRUE,
#'                            ncore = 2)
#'
#'
#'
bayesmsm <- function(ymodel,
                     nvisit,
                     reference = c(rep(0,nvisit)), # An example of never treated
                     comparator = c(rep(1,nvisit)),
                     family = "gaussian", # "gaussian" or "binomial"
                     data,
                     wmean = rep(1, nrow(data)),
                     nboot = 1000,
                     optim_method = 'BFGS',
                     seed = 890123,
                     parallel = TRUE,
                     ncore = 6){

  # load all the required R packages;
  # require(foreach)
  # require(doParallel)
  # require(MCMCpack)
  # require(parallel)
  # if (!require(foreach)){
  #   install.packages("foreach",repos="http://cran.r-project.org")
  #   library(foreach)
  # }
  # if (!require(doParallel)){
  #   install.packages("doParallel",repos="http://cran.r-project.org")
  #   library(doParallel)
  # }
  # if (!require(MCMCpack)){
  #   install.packages("MCMCpack",repos="http://cran.r-project.org")
  #   library(MCMCpack)
  # }

  require(foreach)
  require(doParallel)
  require(MCMCpack)

  # return error message if the input weight vector has different length comparing to the outcome Y;
  if (length(wmean) != nrow(data)) {
    stop("The length of the weight vector does not match the length of Y.")
  }

  # load utility functions
  extract_variables <- function(formula) {
    # Get the terms of the formula
    formula_terms <- terms(formula)

    # Extract the response variable name (if there is one)
    response_variable <- attr(formula_terms, "response")
    response_name <- if (response_variable > 0) {
      all_vars <- all.vars(formula)
      all_vars[response_variable]
    } else {NA}

    # Extract predictor variable names
    predictor_names <- attr(formula_terms, "term.labels")

    # Return a list of response and predictor variables
    list(response = response_name, predictors = predictor_names)
  }

  variables <- extract_variables(ymodel) # Extract variable names from the formula
  Y_name <- variables$response

  Y <- data[[Y_name]]
  A_base <- data.frame(matrix(data = NA,
                              nrow = nrow(data),
                              ncol = length(variables$predictors)))
  for (i in 1:length(variables$predictors)){
    initial_vector <- variables$predictors[i]
    split_vector <- strsplit(initial_vector, ":")
    new_vector <- unlist(split_vector)
    if (length(new_vector)==1){
      A_base[,i] <-  data[, new_vector]
    } else if (length(new_vector)>1){
      A_base[,i] <-  apply(data[, new_vector],1,prod)
    }
  }

  A <- cbind(1, A_base)
  colnames(A)[2:ncol(A)]<- variables$predictors

  wloglik_normal<-function(param,
                           Y,
                           A,
                           weight){
    #number of observations;
    n <- length(Y)
    theta <- param[1:dim(A)[2]] #causal parameters on the mean
    #number of parameter is determined by number of treatment variables, plus intercept;
    sigma <- param[(dim(A)[2]+1)] # the remaining the parameter represent the standard deviation;
    mmat <- as.matrix(A) #design matrix of the causal outcome model, e.g., A = cbind(1, a_1, a_2);
    logl<- -0.5*log(sigma**2) - 0.5*((Y - mmat%*%theta)**2)/(sigma**2)
    wlogl<-sum(weight*logl)
    return(wlogl)
  }

  wloglik_binomial <- function(param,
                               Y,
                               A,
                               weight){
    # number of observations;
    n <- length(Y)
    beta <- param[1:dim(A)[2]] # causal parameters on the log-odds scale
    mmat <- as.matrix(A)
    eta<-mmat %*% beta # linear predictor
    logl <- Y*eta - log(1+exp(eta))
    wlogl<-sum(weight*logl)
    return(wlogl)
  }

  expit <- function(x){exp(x) / (1+exp(x))}

  if (family == "gaussian"){
    wfn = wloglik_normal
    inits1 <- c(rep(0.1, length(A)), 4)  # Default initial values, 4 is for the SD;
  } else if (family == "binomial"){
    wfn = wloglik_binomial
    inits1 <- c(rep(0.1, length(A)))
  } else if (!family %in% c("gaussian","binomial")){
    stop("Current version only handles continuous (gaussian) and binary (binomial) outcomes.")
  }


  # Parallel computing for bootstrapping
  if (parallel == TRUE){

    cl <- parallel::makeCluster(ncore)
    doParallel::registerDoParallel(cl)
    results <- foreach::foreach(i=1:nboot,
                       .combine = 'rbind',
                       .packages = 'MCMCpack') %dopar% {

      calculate_effect <- function(intervention_levels, variables, param_estimates) {
        # Start with the intercept term
        effect<-effect_intercept<-param_estimates[1]

        # Go through each predictor and add its contribution
        for (i in 1:length(variables$predictors)) {
          term <- variables$predictors[i]
          term_variables <- unlist(strsplit(term, ":"))
          term_index <- which(names(param_estimates) == term)

          # Calculate the product of intervention levels for the interaction term
          term_contribution <- param_estimates[term_index]
          for (term_variable in term_variables) {
            var_index <- which(variables$predictors == term_variable)
            term_contribution <- term_contribution * intervention_levels[var_index]
          }

          # Add the term contribution to the effect
          effect <- effect + term_contribution
        }

        return(effect)
      }

      results.it <- matrix(NA, 1, 5) # Result matrix for RD, RR, OR, effect_ref, and effect_comp

      set.seed(seed+i) #define seed;
      alpha <- as.numeric(MCMCpack::rdirichlet(1, rep(1.0, length(Y))))

      maxim <- optim(inits1,
                     fn = wfn,
                     Y = Y,
                     A = A,
                     weight = alpha * wmean,
                     control = list(fnscale = -1),
                     method = optim_method,
                     hessian = FALSE)

      names(maxim$par) <- c("(Intercept)", variables$predictors)

      # Calculate the effects
      effect_ref <- calculate_effect(reference, variables, param_estimates=maxim$par)
      effect_comp <- calculate_effect(comparator, variables, param_estimates=maxim$par)

      # Calculate the ATE
      if (family == "binomial") { # Binary outcomes
        results.it[1,1] <- expit(effect_comp) - expit(effect_ref)  # RD
        results.it[1,2] <- expit(effect_comp) / expit(effect_ref)  # RR
        results.it[1,3] <- (expit(effect_comp) / (1 - expit(effect_comp))) /
          (expit(effect_ref) / (1 - expit(effect_ref)))  # OR
      } else if (family == "gaussian"){ # Continuous outcomes
        results.it[1,1] <- effect_comp - effect_ref  # RD
        results.it[1,2] <- NA  # RR not applicable
        results.it[1,3] <- NA  # OR not applicable
      }

      # Store the reference and comparator effects
      results.it[1,4] <- effect_ref
      results.it[1,5] <- effect_comp

      # combining parallel results;
      cbind(i,results.it) #end of parallel;
    }

    parallel::stopCluster(cl)

    #saving output for the parallel setting;
    if (family == "binomial") {
      return(list(
        RD_mean = mean(results[,2]),
        RR_mean = mean(results[,3]),
        OR_mean = mean(results[,4]),
        RD_sd = sqrt(var(results[,2])),
        RR_sd = sqrt(var(results[,3])),
        OR_sd = sqrt(var(results[,4])),
        RD_quantile = quantile(results[,2], probs = c(0.025, 0.975)),
        RR_quantile = quantile(results[,3], probs = c(0.025, 0.975)),
        OR_quantile = quantile(results[,4], probs = c(0.025, 0.975)),
        bootdata = data.frame(
          effect_reference = results[,5],
          effect_comparator = results[,6],
          RD = results[,2],
          RR = results[,3],
          OR = results[,4]
        ),
        reference = reference,
        comparator = comparator
      ))
    } else if (family == "gaussian"){
      return(list(
        RD_mean = mean(results[,2]),
        RD_sd = sqrt(var(results[,2])),
        RD_quantile = quantile(results[,2], probs = c(0.025, 0.975)),
        bootdata = data.frame(
          effect_reference = results[,5],
          effect_comparator = results[,6],
          RD = results[,2]
        ),
        reference = reference,
        comparator = comparator
      ))
    }
  }

  else if (parallel == FALSE) {

    bootest_RD <- numeric(nboot)
    bootest_RR <- numeric(nboot)
    bootest_OR <- numeric(nboot)
    effect_reference <- numeric(nboot)
    effect_comparator <- numeric(nboot)

    for (j in 1:nboot) {
      alpha <- as.numeric(rdirichlet(1, rep(1.0, length(Y))))

      maxim <- optim(inits1,
                     fn = wfn,
                     Y = Y,
                     A = A,
                     weight = alpha * wmean,
                     control = list(fnscale = -1),
                     method = optim_method,
                     hessian = FALSE)

      names(maxim$par) <- c("(Intercept)", variables$predictors)

      # Calculate the effects
      effect_reference[j] <- calculate_effect(reference, variables, param_estimates=maxim$par)
      effect_comparator[j] <- calculate_effect(comparator, variables, param_estimates=maxim$par)

      # Calculate the ATE
      if (family == "binomial") { # Binary outcomes
        bootest_RD[j] <- expit(effect_comparator[j]) - expit(effect_reference[j])  # RD
        bootest_RR[j] <- expit(effect_comparator[j]) / expit(effect_reference[j])  # RR
        bootest_OR[j] <- (expit(effect_comparator[j]) / (1 - expit(effect_comparator[j]))) /
          (expit(effect_reference[j]) / (1 - expit(effect_reference[j])))  # OR
      } else if (family == "gaussian"){ # Continuous outcomes
        bootest_RD[j] <- effect_comparator[j] - effect_reference[j]  # RD
      }

    }

    #saving output for the non-parallel setting;
    if (family == "binomial") {
      return(list(
        RD_mean = mean(bootest_RD),
        RR_mean = mean(bootest_RR),
        OR_mean = mean(bootest_OR),
        RD_sd = sqrt(var(bootest_RD)),
        RR_sd = sqrt(var(bootest_RR)),
        OR_sd = sqrt(var(bootest_OR)),
        RD_quantile = quantile(bootest_RD, probs = c(0.025, 0.975)),
        RR_quantile = quantile(bootest_RR, probs = c(0.025, 0.975)),
        OR_quantile = quantile(bootest_OR, probs = c(0.025, 0.975)),
        bootdata = data.frame(
          effect_reference = effect_reference,
          effect_comparator = effect_comparator,
          RD = bootest_RD,
          RR = bootest_RR,
          OR = bootest_OR
        ),
        reference = reference,
        comparator = comparator
      ))
    } else {
      return(list(
        RD_mean = mean(bootest_RD),
        RD_sd = sqrt(var(bootest_RD)),
        RD_quantile = quantile(bootest_RD, probs = c(0.025, 0.975)),
        bootdata = data.frame(
          effect_reference = effect_reference,
          effect_comparator = effect_comparator,
          RD = bootest_RD
        ),
        reference = reference,
        comparator = comparator
      ))
    }

  }
}
