#' Bayesian Treatment Effect Weight Estimation Using JAGS
#'
#' This function estimates Bayesian importance sampling weights for time-varying treatment effects using specified models for each treatment time point via JAGS
#'
#' @param trtmodel.list A list of formulas corresponding to each time point with the time-specific treatment variable on the left hand side and pre-treatment covariates to be balanced on the right hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point. Interactions and functions of covariates are allowed.
#' @param data A data set in the form of a data frame containing the variables in "trtmodel.list". This must be a wide data set with exactly one row per unit.
#' @param n.chains Integer specifying the number of MCMC chains to run. Set to 1 for non-parallel computation. For parallel computation, it is required to use at least 2 chains. The default is 2.
#' @param n.iter Integer specifying the total number of iterations for each chain (including burn-in). The default is 25000.
#' @param n.burnin Integer specifying the number of burn-in iterations for each chain. The default is 15000.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampler. The default is 5.
#' @param seed Starting seed for the JAGS model. The default is NULL.
#' @param parallel Logical scalar indicating whether to run the MCMC chains in parallel. The default is TRUE.
#'
#' @return A list of the calculated weights and the JAGS model, where `weights` is a vector of posterior mean weights, computed by taking the average of the weights across all MCMC iterations and `model_string` is a character of the JAGS model based on the input of `trtmodel.list`.
#'
#' @importFrom R2jags jags
#' @importFrom coda mcmc as.mcmc geweke.diag
#' @import parallel
#' @import doParallel
#' @importFrom foreach "%dopar%"
#' @importFrom stats as.formula terms var
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
#' weights <- bayesweight(trtmodel.list = list(
#'                        A1 ~ L1_1 + L2_1,
#'                        A2 ~ L2_2 + L2_2 + A1),
#'                        data = testdata,
#'                        n.chains = 1,
#'                        n.iter = 20,
#'                        n.burnin = 10,
#'                        n.thin = 1,
#'                        seed = 890123,
#'                        parallel = FALSE)
#' summary(weights)
bayesweight <- function(trtmodel.list,
                        data,
                        n.chains = 2,
                        n.iter = 25000,
                        n.burnin = 15000,
                        n.thin = 5,
                        seed = NULL,
                        parallel = TRUE){

  create_marginal_treatment_models <- function(trtmodel.list) {

    # Use lapply to iterate over indices
    trtmodel.list_s <- lapply(seq_along(trtmodel.list), function(index) {
      # Extract the response variable (treatment variable) from each model
      response_var <- all.vars(trtmodel.list[[index]])[1]  # assuming the response is the first variable on the LHS

      # Create the marginal model formula
      if (index == 1) {
        # The first treatment model does not depend on any previous treatments
        formula_s <- as.formula(paste(response_var, "~ 1"))
      } else {
        # Subsequent treatment models depend on all previous treatments
        previous_treatments <- sapply(seq_len(index - 1), function(j) {
          all.vars(trtmodel.list[[j]])[1]
        })
        formula_s <- as.formula(paste(response_var, "~", paste(previous_treatments, collapse = " + ")))
      }

      return(formula_s)
    })

    return(trtmodel.list_s)
  }

  # Generate trtmodel.list_s
  trtmodel.list_s <- create_marginal_treatment_models(trtmodel.list)

  extract_variables_list <- function(input) {
    # Function to extract variables from a single formula
    extract_from_formula <- function(formula) {
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

    # Check if input is a list of formulas
    if (is.list(input)) {
      # Apply the function to each formula in the list
      lapply(input, extract_from_formula)
    } else {
      # Input is a single formula
      extract_from_formula(input)
    }
  }

  trtmodel <- extract_variables_list(trtmodel.list)
  trtmodel_s <- extract_variables_list(trtmodel.list_s)

  # Function to generate and write the JAGS model
  write_jags_model <- function(trtmodel.list) {
    # Extract variable information for each formula
    var_info <- lapply(trtmodel.list, extract_variables_list)

    # Start writing the model string
    model_string <- "model{\n#N = nobs\nfor(i in 1:N){\n"

    # Process each visit
    for (v in seq_along(var_info)) {
      visit <- var_info[[v]]
      response <- visit$response
      predictors <- visit$predictors

      # Marginal treatment assignment model
      model_string <- paste0(model_string,
                             "\n# visit ", v, ";\n",
                             "# marginal treatment assignment model, visit ", v, ";\n",
                             response, "s[i] ~ dbern(p", response, "s[i])\n")

      # Build the logistic model including all previous treatments
      if (v == 1) {
        model_string <- paste0(model_string, "p", response, "s[i] <- ilogit(bs", v, "0)\n")
      } else {
        bs_terms <- paste0("bs", v, "0")
        for (j in 1:(v-1)) {
          prev_response_s <- var_info[[j]]$response
          bs_terms <- paste(bs_terms, " + bs", v, j, "*", prev_response_s, "s[i]", sep = "")
        }
        model_string <- paste0(model_string, "p", response, "s[i] <- ilogit(", bs_terms, ")\n")
      }

      # Conditional treatment assignment model
      model_string <- paste0(model_string,
                             "\n# conditional treatment assignment model, visit ", v, ";\n",
                             response, "[i] ~ dbern(p", response, "[i])\n",
                             "p", response, "[i] <- ilogit(b", v, "0")
      for (p in seq_along(predictors)) {
        model_string <- paste0(model_string, " + b", v, p, "*", predictors[p], "[i]")
      }
      model_string <- paste0(model_string, ")\n")
    }

    # Close loop and add priors
    model_string <- paste0(model_string, "\n# export quantity in full posterior specification;\n",
                           "w[i] <- (", paste(sapply(seq_along(var_info), function(x) paste0("p", var_info[[x]]$response, "s[i]")), collapse = "*"),
                           ")/(", paste(sapply(seq_along(var_info), function(x) paste0("p", var_info[[x]]$response, "[i]")), collapse = "*"), ")\n}\n\n#prior;\n")

    # Add priors for all parameters
    for (v in seq_along(var_info)) {
      visit <- var_info[[v]]
      predictors <- visit$predictors
      num_preds <- length(predictors) + 1  # +1 for intercept

      # bs parameters
      model_string <- paste0(model_string, "bs", v, "0~dnorm(0,.01)\n")
      if (v > 1) {
        for (j in 1:(v-1)) {
          model_string <- paste0(model_string, "bs", v, j, "~dnorm(0,.01)\n")
        }
      }

      # b parameters
      for (p in 0:(num_preds - 1)) {
        model_string <- paste0(model_string, "b", v, p, "~dnorm(0,.01)\n")
      }
    }

    # Add the closing brace for the model block
    model_string <- paste(model_string, "}\n", sep="")

    return(model_string)
  }

  # extracting parameters list;
  jags_model_parameter <- function(trtmodel.list){
    all_parameters <- c()

    # Process each visit
    for (v in seq_along(trtmodel.list)) {
      var_info <- extract_variables_list(trtmodel.list[[v]])
      response <- var_info$response
      predictors <- var_info$predictors

      # Add bs parameters for marginal models
      all_parameters <- c(all_parameters, sprintf("bs%d0", v))
      if (v > 1) {
        extra_bs <- sapply(1:(v-1), function(j) sprintf("bs%d%d", v, j))
        all_parameters <- c(all_parameters, extra_bs)
      }

      # Add b parameters for conditional models
      all_parameters <- c(all_parameters, sprintf("b%d0", v))  # intercept
      for (p in seq_along(predictors)) {
        all_parameters <- c(all_parameters, sprintf("b%d%d", v, p))
      }
    }

    return(unique(all_parameters))
  }

  # Prepare data and parameters for JAGS
  prepare_jags_data <- function(data, trtmodel.list) {
    # Extract variable information from formulas
    variable_info <- lapply(trtmodel.list, extract_variables_list)

    # Initialize the list to pass to JAGS
    jags.data <- list(N = nrow(data))

    # Loop over each formula to add necessary variables
    for (info in variable_info) {
      response <- info$response
      predictors <- info$predictors

      # Add response and its 's' version
      jags.data[[response]] <- data[[response]]
      jags.data[[paste0(response, "s")]] <- data[[response]]  # Assume treatment assignment variable is same as response

      # Add predictors
      for (predictor in predictors) {
        jags.data[[predictor]] <- data[[predictor]]
      }
    }

    return(jags.data)
  }

  # Generate the model string
  model_str <- write_jags_model(trtmodel.list)

  # Prepare JAGS data
  jags.data <- prepare_jags_data(data, trtmodel.list)
  jags.params <- jags_model_parameter(trtmodel.list)

  # Define seed for each chain
  if(!is.null(seed)){
    set.seed(seed)
  }
  new.seed.by.chain <- sample.int(2^30, size=n.chains)

  # section running rjags;
  if (parallel == TRUE) {
    if (n.chains == 1) {
      stop("Parallel MCMC requires at least 2 chains. Computing is running on 1 core per chain.")
    }
    available_cores <- detectCores(logical = FALSE)
    if (n.chains >= available_cores) {
      stop(paste("Parallel MCMC requires 1 core per chain. You have", available_cores, "cores. We recommend using", available_cores - 2, "cores."))
    }
    # Run JAGS model in parallel
    cl <- parallel::makeCluster(n.chains)
    doParallel::registerDoParallel(cl)

    posterior <- foreach::foreach(chain_idx=1:n.chains, .packages=c('R2jags'),
                         .combine='rbind') %dopar%{

                           con <- textConnection(model_str)
                           jagsfit <- R2jags::jags(data = jags.data,
                                           parameters.to.save = jags.params,
                                           model.file = con,
                                           n.chains = 1,
                                           n.iter = n.iter,
                                           n.burnin = n.burnin,
                                           n.thin = n.thin,
                                           jags.seed = new.seed.by.chain[chain_idx])
                           close(con)

                           # Combine MCMC output from multiple chains
                           out.mcmc <- coda::as.mcmc(jagsfit)
                           return(do.call(rbind, lapply(out.mcmc, as.matrix)))

                         }

    parallel::stopCluster(cl)

  } else if (parallel == FALSE) {

    if (n.chains != 1) {
      stop("Non-parallel MCMC requires 1 chain.")
    }

    con <- textConnection(model_str)
    # Run JAGS model without parallel computing
    jagsfit <- R2jags::jags(data = jags.data,
                    parameters.to.save = jags.params,
                    model.file = con,
                    n.chains = 1,
                    n.iter = n.iter,
                    n.burnin = n.burnin,
                    n.thin = n.thin,
                    jags.seed = new.seed.by.chain[1])
    close(con)

    # Extract MCMC output
    out.mcmc <- coda::as.mcmc(jagsfit)
    diagnostics <- geweke.diag(out.mcmc)

    # Check diagnostics for convergence issues
    significant_indices <- which(abs(diagnostics[[1]]$z) > 1.96)
    if (length(significant_indices) > 0) {
      warning("Some parameters have not converged with Geweke index > 1.96. More iterations may be needed.")
    }

    posterior <- as.matrix(out.mcmc[[1]])

  }


  # number of parameters for this model is 5 and design matrix is 1 variables
  # looping through each treatment model 1 visit at a time;
  expit <- function(x){exp(x) / (1+exp(x))}

  n_visits <- length(trtmodel)
  n_posterior <- dim(posterior)[1]
  n_obs <- nrow(data)

  # Initialize arrays for storing probabilities
  psc <- array(dim = c(n_visits, n_posterior, n_obs))
  psm <- array(dim = c(n_visits, n_posterior, n_obs))

  # Calculate probabilities for each visit
  parameter_map <- colnames(posterior)

  for (nvisit in 1:n_visits) {
    predictors_c <- trtmodel[[nvisit]]$predictors
    predictors_s <- trtmodel_s[[nvisit]]$predictors

    design_matrix_c <- cbind(1, data[, predictors_c, drop = FALSE])
    design_matrix_s <- cbind(1, data[, predictors_s, drop = FALSE])

    beta_indices_c <- match(c(sprintf("b%d0", nvisit),
                              sapply(1:length(predictors_c), function(p) sprintf("b%d%d", nvisit, p))),
                            parameter_map) # This produces the correct indices that corresponds to the predictors in the treatment models

    beta_indices_s <- match(c(sprintf("bs%d0", nvisit),
                              if (nvisit > 1) sapply(1:(nvisit-1), function(j) sprintf("bs%d%d", nvisit, j))),
                            parameter_map)

    for (j in 1:n_posterior) {
      psc[nvisit, j, ] <- expit(posterior[j, beta_indices_c] %*% t(design_matrix_c))
      psm[nvisit, j, ] <- expit(posterior[j, beta_indices_s] %*% t(design_matrix_s))
    }
  }

  # Calculate the product of probabilities across visits for each posterior and observation
  numerator <- apply(psm, c(2, 3), prod)  # Apply 'prod' across the first dimension (visits)
  denominator <- apply(psc, c(2, 3), prod)  # Same for psc

  # Calculate weights by element-wise division
  weights <- numerator / denominator  # Resulting in a 40 (posterior samples) x 1000 (observations) matrix

  # Mean weight across all observations
  wmean <- colMeans(weights)

  # Return the weights and the model string
  return(list(weights = wmean, model_string = model_str))

}
