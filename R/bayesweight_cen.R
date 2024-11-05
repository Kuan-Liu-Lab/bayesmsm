#' Bayesian Weight Estimation for Censored Data
#'
#' This function computes posterior mean weights using Bayesian estimation for treatment models and censoring models across multiple time points. The models can be run in parallel to estimate the weights needed for causal inference with censored data.
#'
#' @param trtmodel.list A list of formulas corresponding to each time point with the time-specific treatment variable on the left-hand side and pre-treatment covariates to be balanced on the right-hand side. Interactions and functions of covariates are allowed.
#' @param cenmodel.list A list of formulas for the censoring data at each time point, with censoring indicators on the left-hand side and covariates on the right-hand side.
#' @param data A data frame containing the variables in the models (treatment, censoring, and covariates).
#' @param n.iter Number of iterations to run the MCMC algorithm in JAGS.
#' @param n.burnin Number of iterations to discard as burn-in in the MCMC algorithm.
#' @param n.thin Thinning rate for the MCMC samples.
#' @param parallel Logical, indicating whether to run the MCMC sampling in parallel (default is `FALSE`).
#' @param n.chains Number of MCMC chains to run. If parallel is `TRUE`, this specifies the number of chains run in parallel.
#' @param seed A seed for random number generation to ensure reproducibility of the MCMC.
#'
#' @return A vector of posterior mean weights, computed by taking the average of the weights across all MCMC iterations.
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
#' simdat_cen <- read.csv(system.file("extdata",
#'                                    "sim_causal.csv",
#'                                    package = "bayesmsm"))
#' weights_cen <- bayesweight_cen(trtmodel.list = list(A1 ~ L11 + L21,
#'                                                     A2 ~ L11 + L21 + L12 +
#'                                                          L22 + A1,
#'                                                     A3 ~ L11 + L21 + L12 +
#'                                                          L22 + A1 + L13 +
#'                                                          L23 + A2),
#'                                cenmodel.list = list(C1 ~ L11 + L21,
#'                                                     C2 ~ L11 + L21 + A1,
#'                                                     C3 ~ L11 + L21 + A1 +
#'                                                          L12 + L22 + A2),
#'                                data = simdat_cen,
#'                                n.iter = 1500,
#'                                n.burnin = 500,
#'                                n.thin = 1,
#'                                parallel = FALSE,
#'                                n.chains = 1,
#'                                seed = 890123)
bayesweight_cen <- function(trtmodel.list = list(A1 ~ L11 + L21,
                                                 A2 ~ L11 + L21 + L12 + L22 + A1,
                                                 A3 ~ L11 + L21 + L12 + L22 + A1 + L13 + L23 + A2),
                            cenmodel.list = list(C1 ~ L11 + L21,
                                                 C2 ~ L11 + L21 + A1,
                                                 C3 ~ L11 + L21 + A1 + L12 + L22 + A2),
                            data,
                            n.iter = 2500,
                            n.burnin = 1500,
                            n.thin = 5,
                            parallel = FALSE,
                            n.chains = 1,
                            seed = 890123) {


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

    # # Initialize the list for the marginal treatment models
    # trtmodel.list_s <- list()
    #
    # # Loop through each model in the original list
    # for (i in seq_along(trtmodel.list)) {
    #   # Extract the response variable (treatment variable) from each model
    #   response_var <- all.vars(trtmodel.list[[i]])[1]  # assuming the response is the first variable on the LHS
    #
    #   # Create the marginal model formula
    #   if (i == 1) {
    #     # The first treatment model does not depend on any previous treatments
    #     formula_s <- as.formula(paste(response_var, "~ 1"))
    #   } else {
    #     # Subsequent treatment models depend on all previous treatments
    #     previous_treatments <- sapply(seq_len(length(trtmodel.list) - 1), function(j) {
    #       all.vars(trtmodel.list[[j]])[1]
    #     })
    #     formula_s <- as.formula(paste(response_var, "~", paste(previous_treatments, collapse = " + ")))
    #   }
    #
    #   # Append the new formula to the list
    #   trtmodel.list_s[[i]] <- formula_s
    # }
    #
    return(trtmodel.list_s)
  }

  # Generate trtmodel.list_s
  trtmodel.list_s <- create_marginal_treatment_models(trtmodel.list)
  cenmodel.list_s <- create_marginal_treatment_models(cenmodel.list)

  # Extract variables from treatment and censoring models
  extract_variables_list <- function(input) {
    extract_from_formula <- function(formula) {
      formula_terms <- terms(formula)
      response_variable <- attr(formula_terms, "response")
      response_name <- if (response_variable > 0) {
        all_vars <- all.vars(formula)
        all_vars[response_variable]
      } else {NA}
      predictor_names <- attr(formula_terms, "term.labels")
      list(response = response_name, predictors = predictor_names)
    }
    if (is.list(input)) {
      lapply(input, extract_from_formula)
    } else {
      extract_from_formula(input)
    }
  }

  trtmodel <- extract_variables_list(trtmodel.list)
  trtmodel_s <- extract_variables_list(trtmodel.list_s)
  cenmodel <- extract_variables_list(cenmodel.list)
  cenmodel_s <- extract_variables_list(cenmodel.list_s)

  # Define JAGS model for treatment and censoring
  write_jags_model <- function(trtmodel.list, cenmodel.list) {
    var_info_trt <- lapply(trtmodel.list, extract_variables_list)
    var_info_cen <- lapply(cenmodel.list, extract_variables_list)

    model_string <- "model{\n"

    all_parameters <- c()

    for (v in seq_along(var_info_trt)) {
      visit_trt <- var_info_trt[[v]]
      response_trt <- visit_trt$response
      predictors_trt <- visit_trt$predictors
      visit_cen <- var_info_cen[[v]]
      response_cen <- visit_cen$response
      predictors_cen <- visit_cen$predictors

      model_string <- paste0(model_string, "\nfor (i in 1:N", v, ") {\n")

      # Conditional treatment assignment model
      model_string <- paste0(model_string,
                             "\n# conditional model;\n",
                             response_trt, "[i] ~ dbern(p", v, "[i])\n",
                             "logit(p", v, "[i]) <- b", v, "0")
      for (p in seq_along(predictors_trt)) {
        model_string <- paste0(model_string, " + b", v, p, "*", predictors_trt[p], "[i]")
        all_parameters <- c(all_parameters, sprintf("b%d%d", v, p))
      }
      model_string <- paste0(model_string, "\n")

      # Censoring model
      model_string <- paste0(model_string,
                             response_cen, "[i] ~ dbern(cp", v, "[i])\n",
                             "logit(cp", v, "[i]) <- s", v, "0")
      for (p in seq_along(predictors_cen)) {
        model_string <- paste0(model_string, " + s", v, p, "*", predictors_cen[p], "[i]")
        all_parameters <- c(all_parameters, sprintf("s%d%d", v, p))
      }
      model_string <- paste0(model_string, "\n")

      # Marginal treatment assignment model
      model_string <- paste0(model_string,
                             "\n# marginal model;\n",
                             response_trt, "s[i] ~ dbern(p", v, "s[i])\n",
                             "logit(p", v, "s[i]) <- bs", v, "0")
      all_parameters <- c(all_parameters, sprintf("bs%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          model_string <- paste0(model_string, " + bs", v, j, "*", prev_response_trt, "s[i]")
          all_parameters <- c(all_parameters, sprintf("bs%d%d", v, j))
        }
      }
      model_string <- paste0(model_string, "\n")

      # Marginal censoring model
      model_string <- paste0(model_string,
                             response_cen, "s[i] ~ dbern(cp", v, "s[i])\n",
                             "logit(cp", v, "s[i]) <- ts", v, "0")
      all_parameters <- c(all_parameters, sprintf("ts%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          model_string <- paste0(model_string, " + ts", v, j, "*", prev_response_trt, "s[i]")
          all_parameters <- c(all_parameters, sprintf("ts%d%d", v, j))
        }
      }
      model_string <- paste0(model_string, "\n}\n")
    }

    # Priors section
    model_string <- paste0(model_string, "\n# Priors\n")
    for (v in seq_along(var_info_trt)) {
      num_preds_trt <- length(var_info_trt[[v]]$predictors)
      num_preds_cen <- length(var_info_cen[[v]]$predictors)

      # Treatment priors
      for (p in 0:num_preds_trt) {
        model_string <- paste0(model_string, "b", v, p, " ~ dunif(-10, 10)\n")
        all_parameters <- c(all_parameters, sprintf("b%d%d", v, p))
      }

      # Censoring priors
      for (p in 0:num_preds_cen) {
        model_string <- paste0(model_string, "s", v, p, " ~ dunif(-10, 10)\n")
        all_parameters <- c(all_parameters, sprintf("s%d%d", v, p))
      }

      # Marginal treatment priors
      model_string <- paste0(model_string, "bs", v, "0 ~ dunif(-10, 10)\n")
      all_parameters <- c(all_parameters, sprintf("bs%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          model_string <- paste0(model_string, "bs", v, j, " ~ dunif(-10, 10)\n")
          all_parameters <- c(all_parameters, sprintf("bs%d%d", v, j))
        }
      }

      # Marginal censoring priors
      model_string <- paste0(model_string, "ts", v, "0 ~ dunif(-10, 10)\n")
      all_parameters <- c(all_parameters, sprintf("ts%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          model_string <- paste0(model_string, "ts", v, j, " ~ dunif(-10, 10)\n")
          all_parameters <- c(all_parameters, sprintf("ts%d%d", v, j))
        }
      }
    }

    # Add the closing brace for the model block
    model_string <- paste0(model_string, "}\n")

    # Write the finalized model string to a file
    cat(model_string, file = "censoring_model.txt")

    return(unique(all_parameters))
  }

  write_jags_model(trtmodel.list, cenmodel.list)

  # Prepare data for JAGS
  prepare_jags_data <- function(data, trtmodel.list, cenmodel.list) {
    variable_info_trt <- lapply(trtmodel.list, extract_variables_list)
    variable_info_cen <- lapply(cenmodel.list, extract_variables_list)

    # Collect all variables from trtmodel and cenmodel
    all_vars <- unique(unlist(lapply(c(variable_info_trt, variable_info_cen), function(info) {
      c(info$response, info$predictors)
    })))

    # Initialize the list to pass to JAGS
    jags.data <- list()

    # Add each column of data to the jags.data list, with NA values removed
    for (col in all_vars) {
      if (!is.na(col)) {
        jags.data[[col]] <- data[[col]][!is.na(data[[col]])]
      }
    }

    # Determine N for each visit based on the length of non-NA values of the response variable
    for (v in seq_along(variable_info_trt)) {
      response_trt <- variable_info_trt[[v]]$response
      jags.data[[paste0("N", v)]] <- sum(!is.na(data[[response_trt]]))
      jags.data[[paste0("A", v, "s")]] <- data[[paste0("A", v)]][!is.na(data[[paste0("A", v)]])] # Remove the NAs
      jags.data[[paste0("C", v, "s")]] <- data[[paste0("C", v)]][!is.na(data[[paste0("C", v)]])]
    }

    return(jags.data)
  }

  jags.data <- prepare_jags_data(data, trtmodel.list, cenmodel.list)
  jags.params <- write_jags_model(trtmodel.list, cenmodel.list)

  # Run JAGS model
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

    jags.model.wd <- paste(getwd(), '/censoring_model.txt',sep='')

    posterior <- foreach::foreach(chain_idx=1:n.chains, .packages=c('R2jags'),
                         .combine='rbind') %dopar%{

                           set.seed(seed+chain_idx) #define seed;
                           jagsfit <- jags(data = jags.data,
                                           parameters.to.save = jags.params,
                                           model.file = jags.model.wd,
                                           n.chains = 1,
                                           n.iter = n.iter,
                                           n.burnin = n.burnin,
                                           n.thin = n.thin)
                           # Combine MCMC output from multiple chains
                           out.mcmc <- as.mcmc(jagsfit)
                           return(do.call(rbind, lapply(out.mcmc, as.matrix)))

                         }

    parallel::stopCluster(cl)

  } else if (parallel == FALSE) {
    if (n.chains != 1) {
      stop("Non-parallel MCMC requires exactly 1 chain.")
    }
    jagsfit <- jags(data = jags.data,
                    parameters.to.save = jags.params,
                    model.file = "censoring_model.txt",
                    n.chains = 1,
                    n.iter = n.iter,
                    n.burnin = n.burnin,
                    n.thin = n.thin,
                    jags.seed = seed)

    out.mcmc <- as.mcmc(jagsfit)
    posterior <- as.matrix(out.mcmc[[1]])

    diagnostics <- geweke.diag(out.mcmc)

    # Check diagnostics for convergence issues
    significant_indices <- which(abs(diagnostics[[1]]$z) > 1.96)
    if (length(significant_indices) > 0) {
      warning("Some parameters have not converged with Geweke index > 1.96. More iterations may be needed.")
    }

  }

  expit <- function(x){exp(x) / (1+exp(x))}


  # Initialize arrays for storing probabilities
  n_visits <- length(trtmodel)
  n_posterior <- dim(posterior)[1]
  # Calculate the number of observations for the last visit (i.e. complete data)
  n_obs <- dim(data)[1]

  psc <- array(dim = c(n_visits, n_posterior, n_obs))
  psm <- array(dim = c(n_visits, n_posterior, n_obs))
  csc <- array(dim = c(n_visits, n_posterior, n_obs))
  csm <- array(dim = c(n_visits, n_posterior, n_obs))

  # Calculate probabilities for each visit
  parameter_map <- colnames(posterior)

  for (nvisit in 1:n_visits) {
    predictors_c <- trtmodel[[nvisit]]$predictors
    predictors_s <- trtmodel_s[[nvisit]]$predictors
    predictors_cen_c <- cenmodel[[nvisit]]$predictors
    predictors_cen_s <- cenmodel_s[[nvisit]]$predictors

    design_matrix_c <- cbind(1, data[predictors_c])
    design_matrix_s <- cbind(1, data[predictors_s])
    design_matrix_cen_c <- cbind(1, data[predictors_cen_c])
    design_matrix_cen_s <- cbind(1, data[predictors_cen_s])

    beta_indices_c <- match(c(sprintf("b%d0", nvisit),
                              sapply(1:length(predictors_c), function(p) sprintf("b%d%d", nvisit, p))),
                            parameter_map)

    beta_indices_s <- match(c(sprintf("bs%d0", nvisit),
                              if (nvisit > 1) sapply(1:length(predictors_s), function(p) sprintf("bs%d%d", nvisit, p)) else NULL),
                            parameter_map)

    beta_indices_cen_c <- match(c(sprintf("s%d0", nvisit),
                                  sapply(1:length(predictors_cen_c), function(p) sprintf("s%d%d", nvisit, p))),
                                parameter_map)

    beta_indices_cen_s <- match(c(sprintf("ts%d0", nvisit),
                                  if (nvisit > 1) sapply(1:length(predictors_cen_s), function(p) sprintf("ts%d%d", nvisit, p)) else NULL),
                                parameter_map)

    for (j in 1:n_posterior) {
      psc[nvisit, j, ] <- expit(posterior[j, beta_indices_c] %*% t(design_matrix_c))
      psm[nvisit, j, ] <- expit(posterior[j, beta_indices_s] %*% t(design_matrix_s))
      csc[nvisit, j, ] <- expit(posterior[j, beta_indices_cen_c] %*% t(design_matrix_cen_c))
      csm[nvisit, j, ] <- expit(posterior[j, beta_indices_cen_s] %*% t(design_matrix_cen_s))
    }
  }

  numerator_trt <- apply(psm, c(2, 3), prod)
  denominator_trt <- apply(psc, c(2, 3), prod)
  numerator_cen <- apply(csm, c(2, 3), prod)
  denominator_cen <- apply(csc, c(2, 3), prod)

  weights <- (numerator_trt / denominator_trt) * (numerator_cen / denominator_cen)
  wmean <- colMeans(weights)

  return(wmean)

}
