#' Bayesian Marginal Structural Model Bootstrap Estimation
#'
#' This function performs Bayesian non-parametric bootstrap to estimate causal
#' effects in Bayesian marginal structural models.It supports both continuous
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
#' @import foreach
#' @import doParallel
#' @import MCMCpack
#'
#' @export
#'
#' @examples Using the testdata in the data folder:
#'
#' # Continuous outcome
#' testdata <- readr::read_csv("inst/extdata/continuous_outcome_data.csv")
#' model <- bayesmsm(ymodel = y ~ a_1+a_2,
#'                            nvisit = 2,
#'                            reference = c(rep(0,2)),
#'                            comparator = c(rep(1,2)),
#'                            family = "gaussian",
#'                            data = testdata,
#'                            wmean = rep(1, 1000),
#'                            nboot = 1000,
#'                            optim_method = "BFGS",
#'                            parallel = FALSE,
#'                            ncore = 6)
#'
#'
#'
bayesmsm <- function(ymodel,
                     nvisit,
                     reference = c(rep(0,nvisit)), # An example of never treated
                     comparator = c(rep(1,nvisit)),
                     family = "gaussian", # "gaussian" or "binomial"
                     data,
                     wmean = rep(1, 1000),
                     nboot = 1000,
                     optim_method = 'BFGS',
                     parallel = TRUE,
                     ncore = 4){

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
    beta <- param[1:dim(A)[2]] # causal parameters on the log-odds scale (no sigma for binomial?)
    mmat <- as.matrix(A)
    eta<-mmat %*% beta # linear predictor
    p <- 1 / (1 + exp(-eta))
    logl <- Y * log(p + 0.0001) + (1 - Y) * log(1 - p + 0.0001)
    wlogl<-sum(weight*logl)
    return(wlogl)
  }

  if (family == "gaussian"){
    wfn = wloglik_normal
    inits1 <- c(rep(0.1, length(A)), 4)  # Default initial values, 4 is for the SD;
  } else if (family == "binomial"){
    wfn = wloglik_binomial
    inits1 <- c(rep(0.1, length(A)))
  } else if (!family %in% c("gaussian","binomial")){
    stop("Current version only handles continuous (gaussian) and binary (binomial) outcomes.")
  }


  #parallel computing only for this bootstrap step;
  if (parallel == TRUE){
    numCores <- ncore
    registerDoParallel(cores = numCores)

    results <- foreach(i=1:nboot, .combine = 'rbind') %dopar% {

      results.it <- matrix(NA, 1, 3) #result matrix, three columns for bootest, effect_ref, and effect_comp;

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
      results.it[1,1] <- calculate_effect(reference, variables, param_estimates=maxim$par)
      results.it[1,2] <- calculate_effect(comparator, variables, param_estimates=maxim$par)
      # Calculate the ATE
      results.it[1,3] <- results.it[1,1] - results.it[1,2]

      # combining parallel results;
      cbind(i,results.it) #end of parallel;
    }

    #saving output for the parallel setting;
    return(list(
      mean = mean(results[,4]),
      sd = sqrt(var(results[,4])),
      quantile = quantile(results[,4], probs = c(0.025, 0.975)),
      bootdata <- data.frame(results[,-1]),
      reference = reference,
      comparator = comparator
    ))

  }

  else if (parallel == FALSE) {

    bootest <- numeric(nboot)
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
      bootest[j] <- effect_comparator[j] - effect_reference[j]

    }

    #saving output for the non-parallel setting;
    return(list(
      mean = mean(bootest),
      sd = sqrt(var(bootest)),
      quantile = quantile(bootest, probs = c(0.025, 0.975)),
      bootdata = data.frame(effect_reference, effect_comparator, ATE=bootest),
      reference = reference,
      comparator = comparator
    ))

  }
}
