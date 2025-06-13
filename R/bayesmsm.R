#' Bayesian Marginal Structural Model Bootstrap Estimation
#'
#' This function performs Bayesian non-parametric bootstrap to estimate causal
#' effects in Bayesian marginal structural models. It supports both continuous
#' (Gaussian) and binary (binomial) outcome variables
#'
#' @param ymodel Model statement for the outcome variable.
#' @param nvisit Number of visits or time points to simulate.
#' @param reference Vector denoting the intervention to be used as the reference across all visits for calculating the risk ratio and risk difference. The default is a vector of all 0's with length nvisit (i.e. never treated).
#' @param comparator Vector denoting the intervention to be used as the comparator across all visits for calculating the risk ratio and risk difference. The default is a vector of all 1's with length nvisit (i.e. always treated).
#' @param treatment_effect_type Character string specifying the type of treatment effect to estimate. Options are "sq" for sequential treatment effects, which estimates effects for specific treatment sequences across visits, and "cum" for cumulative treatment effects, which assumes a single cumulative treatment variable representing the total exposure. The default is "sq".
#' @param family Character string specifying the outcome distribution family. The possible distributions are: "Gaussian" (default) for continuous outcomes, and "binomial" for binary outcomes.
#' @param data Data table containing the variable names in the outcome model.
#' @param wmean Vector of treatment assignment weights. The default is rep(1, nrow(data)).
#' @param nboot Integer specifying the number of bootstrap iterations. The default is 1000.
#' @param optim_method Character string specifying the optimization method to be used. The default is "BFGS".
#' @param seed Starting seed for simulations and bootstrapping. The default is NULL.
#' @param parallel Logical scalar indicating whether to parallelize bootstrapping to multiple cores. The default is TRUE.
#' @param ncore Integer specifying the number of CPU cores to use in parallel simulation. This argument is required when parallel is set to TRUE, and the default is 4.
#'
#' @return It returns an object of class "bayesmsm" that contains the information about the data, model, etc.
#' An object of class "bayesmsm" is a list containing at least the following components: "mean", the mean of the bootstrap estimates; "sd", the standard deviation of the bootstrap estimates; "quantile", the 95\% quantiles of the bootstrap estimates; "bootdata", a data frame of bootstrapped estimates; "reference", the reference intervention level and "comparator", the comparison intervention level
#'
#' @importFrom foreach "%dopar%"
#' @import doParallel
#' @import parallel
#' @importFrom MCMCpack rdirichlet
#' @importFrom stats as.formula density optim quantile sd terms var
#'
#' @export
#'
#' @examples
#' # 1) Specify simple treatment‐assignment models
#' amodel <- list(
#'   c("(Intercept)" =  0, "L1_1" =  0.5, "L2_1" = -0.5),
#'   c("(Intercept)" =  0, "L1_2" =  0.5, "L2_2" = -0.5, "A_prev" = 0.3)
#' )
#' # 2) Specify a continuous‐outcome model
#' ymodel <- c("(Intercept)" = 0,
#'             "A1"         = 0.2,
#'             "A2"         = 0.3,
#'             "L1_2"       = 0.1,
#'             "L2_2"       = -0.1)
#' # 3) Simulate without right‐censoring
#' testdata <- simData(
#'   n                = 200,
#'   n_visits         = 2,
#'   covariate_counts = c(2, 2),
#'   amodel           = amodel,
#'   ymodel           = ymodel,
#'   y_type           = "continuous",
#'   right_censor     = FALSE,
#'   seed             = 123)
#' model <- bayesmsm(ymodel = Y ~ A1 + A2,
#'                   nvisit = 2,
#'                   reference = c(rep(0,2)),
#'                   comparator = c(rep(1,2)),
#'                   treatment_effect_type = "sq",
#'                   family = "binomial",
#'                   data = testdata,
#'                   wmean = rep(1,200),
#'                   nboot = 10,
#'                   optim_method = "BFGS",
#'                   seed = 890123,
#'                   parallel = FALSE)
bayesmsm <- function(ymodel,
                     nvisit,
                     reference = c(rep(0,nvisit)), # An example of never treated
                     comparator = c(rep(1,nvisit)),
                     treatment_effect_type = "sq", # "sq" or "cum"
                     family = "gaussian", # "gaussian" or "binomial"
                     data,
                     wmean = rep(1, nrow(data)),
                     nboot = 1000,
                     optim_method = 'BFGS',
                     seed = NULL,
                     parallel = TRUE,
                     ncore = 4){

  # return error message if the input weight vector has different length comparing to the outcome Y;
  if (length(wmean) != nrow(data)) {
    stop("The length of the weight vector does not match the length of Y.")
  }
  # Check treatment_effect_type validity
  if (!treatment_effect_type %in% c("cum", "sq")) {
    stop("Invalid treatment_effect_type. Choose either 'cum' or 'sq'.")
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

  # Define seed and bootstrap sampling weights matrix (alpha)
  if(!is.null(seed)){
    set.seed(seed)
  }
  alpha <- matrix(NA, nrow=length(Y), ncol=nboot)
  for (i in 1:nboot){
    alpha[,i] <- as.numeric(MCMCpack::rdirichlet(1, rep(1.0, length(Y))))
  }

  # Check for cumulative treatment effect conditions
  if (treatment_effect_type == "cum") {
    # Ensure ymodel has only one predictor
    if (length(variables$predictors) != 1) {
      stop("A cumulative treatment effect is specified but the model does not have a single predictor.")
    }

    # Ensure the single predictor in the dataset contains values > 1
    predictor_var <- variables$predictors[1]
    if (!predictor_var %in% colnames(data)) {
      stop(paste("The predictor variable", predictor_var, "is not found in the dataset."))
    }
    if (!any(data[[predictor_var]] > 1)) {
      stop("A cumulative treatment effect is specified but the predictor variable does not contain any value > 1.")
    }
  }

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

      calculate_effect <- function(intervention_levels, variables, param_estimates, treatment_effect_type) {
        if (treatment_effect_type == "cum") {
          # For cumulative treatment, only consider b1
          effect <- param_estimates[1] * intervention_levels[1]
        } else {
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
        }
        return(effect)
      }

      results.it <- matrix(NA, 1, 5) # Result matrix for RD, RR, OR, effect_ref, and effect_comp

      maxim <- optim(inits1,
                     fn = wfn,
                     Y = Y,
                     A = A,
                     weight = alpha[,i] * wmean,
                     control = list(fnscale = -1),
                     method = optim_method,
                     hessian = FALSE)

      names(maxim$par) <- c("(Intercept)", variables$predictors)

      # Calculate the effects
      effect_ref <- calculate_effect(reference, variables, param_estimates=maxim$par, treatment_effect_type)
      effect_comp <- calculate_effect(comparator, variables, param_estimates=maxim$par, treatment_effect_type)

      if (treatment_effect_type == "sq") {
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
      } else if (treatment_effect_type == "cum") {
        # Calculate the ATE
        if (family == "binomial") { # Binary outcomes
          results.it[1,1] <- NA  # RD not applicable
          results.it[1,2] <- NA  # RR not applicable
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
      }

      # combining parallel results;
      cbind(i,results.it) #end of parallel;
    }

    parallel::stopCluster(cl)

    #saving output for the parallel setting;
    if (family == "binomial" && treatment_effect_type == "cum") {
      return(list(
        OR_mean = mean(results[, 4]),
        OR_sd = sqrt(var(results[, 4])),
        OR_quantile = quantile(results[, 4], probs = c(0.025, 0.975)),
        bootdata = data.frame(
          effect_reference = results[, 5],
          effect_comparator = results[, 6],
          OR = results[, 4]
        ),
        reference = reference,
        comparator = comparator
      ))
    } else if (family == "gaussian" && treatment_effect_type == "cum") {
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
    } else if (treatment_effect_type == "sq") {
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
  }

  else if (parallel == FALSE) {

    bootest_RD <- numeric(nboot)
    bootest_RR <- numeric(nboot)
    bootest_OR <- numeric(nboot)
    effect_reference <- numeric(nboot)
    effect_comparator <- numeric(nboot)

    for (j in 1:nboot) {

      maxim <- optim(inits1,
                     fn = wfn,
                     Y = Y,
                     A = A,
                     weight = alpha[,i] * wmean,
                     control = list(fnscale = -1),
                     method = optim_method,
                     hessian = FALSE)

      names(maxim$par) <- c("(Intercept)", variables$predictors)

      # Calculate the effects
      effect_reference[j] <- calculate_effect(reference, variables, param_estimates=maxim$par, treatment_effect_type)
      effect_comparator[j] <- calculate_effect(comparator, variables, param_estimates=maxim$par, treatment_effect_type)

      if (treatment_effect_type == "sq") {
        # Calculate the ATE
        if (family == "binomial") { # Binary outcomes
          bootest_RD[j] <- expit(effect_comparator[j]) - expit(effect_reference[j])  # RD
          bootest_RR[j] <- expit(effect_comparator[j]) / expit(effect_reference[j])  # RR
          bootest_OR[j] <- (expit(effect_comparator[j]) / (1 - expit(effect_comparator[j]))) /
            (expit(effect_reference[j]) / (1 - expit(effect_reference[j])))  # OR
        } else if (family == "gaussian"){ # Continuous outcomes
          bootest_RD[j] <- effect_comparator[j] - effect_reference[j]  # RD
        }
      } else if (treatment_effect_type == "cum") {
        if (family == "binomial") { # Binary outcomes
          bootest_OR[j] <- (expit(effect_comparator[j]) / (1 - expit(effect_comparator[j]))) /
            (expit(effect_reference[j]) / (1 - expit(effect_reference[j])))  # OR
        } else if (family == "gaussian") { # Continuous outcomes
          bootest_RD[j] <- effect_comparator[j] - effect_reference[j]  # RD
        }
      }
    }

    #saving output for the non-parallel setting;
    if (treatment_effect_type == "sq") {
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
    } else if (treatment_effect_type == "cum") {
      if (family == "binomial") {
        return(list(
          OR_mean = mean(bootest_OR),
          OR_sd = sqrt(var(bootest_OR)),
          OR_quantile = quantile(bootest_OR, probs = c(0.025, 0.975)),
          bootdata = data.frame(
            effect_reference = effect_reference,
            effect_comparator = effect_comparator,
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
}
