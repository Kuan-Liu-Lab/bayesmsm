#' Bayesian Treatment Effect Weight Estimation for Censored Data
#'
#' This function estimates Bayesian importance sampling weights for treatment models and censoring models across multiple time points.
#' The models can be run in parallel to estimate the weights for censored data.
#'
#' @param trtmodel.list A list of formulas corresponding to each time point with the time-specific treatment variable on the left-hand side and pre-treatment covariates to be balanced on the right-hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point. Interactions and functions of covariates are allowed.
#' @param cenmodel.list A list of formulas for the censored data at each time point, with censoring indicators on the left-hand side and covariates on the right-hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point.
#' @param data A data set in the form of a data frame containing the variables in `trtmodel.list` and `cenmodel.list`. This must be a wide data set with exactly one row per unit.
#' @param n.iter Integer specifying the total number of iterations for each chain (including burn-in). The default is 25000.
#' @param n.burnin Integer specifying the number of burn-in iterations for each chain. The default is 15000.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampler. The default is 5.
#' @param parallel Logical scalar indicating whether to run the MCMC chains in parallel. The default is TRUE.
#' @param n.chains Integer specifying the number of MCMC chains to run. Set to 1 for non-parallel computation. For parallel computation, it is required to use at least 2 chains. The default is 2.
#' @param seed Starting seed for the JAGS model. The default is NULL.
#' @param save_jags_model_file Logical; if TRUE, writes the model to outputfile. Default is FALSE.
#' @param output_file File name to save the JAGS model (if save_to_file = TRUE).
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
#'                                n.iter = 200,
#'                                n.burnin = 100,
#'                                n.thin = 1,
#'                                parallel = FALSE,
#'                                n.chains = 1,
#'                                seed = 890123,
#'                                save_jags_model_file = FALSE)
bayesweight_cen <- function(trtmodel.list,
                            cenmodel.list,
                            data,
                            n.iter = 25000,
                            n.burnin = 15000,
                            n.thin = 5,
                            parallel = TRUE,
                            n.chains = 2,
                            seed = NULL,
                            save_jags_model_file = FALSE,
                            output_file = "treatment_censoring_model.txt") {

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
  write_jags_model <- function(trtmodel.list,
                               cenmodel.list,
                               save_jags_model_file=FALSE,
                               output_file=output_file) {

    var_info_trt <- lapply(trtmodel.list, extract_variables_list)
    var_info_cen <- lapply(cenmodel.list, extract_variables_list)

    model_string <- "model{\n"

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
      }
      model_string <- paste0(model_string, "\n")

      # Censoring model
      model_string <- paste0(model_string,
                             response_cen, "[i] ~ dbern(cp", v, "[i])\n",
                             "logit(cp", v, "[i]) <- s", v, "0")
      for (p in seq_along(predictors_cen)) {
        model_string <- paste0(model_string, " + s", v, p, "*", predictors_cen[p], "[i]")
      }
      model_string <- paste0(model_string, "\n")

      # Marginal treatment assignment model
      model_string <- paste0(model_string,
                             "\n# marginal model;\n",
                             response_trt, "s[i] ~ dbern(p", v, "s[i])\n",
                             "logit(p", v, "s[i]) <- bs", v, "0")
      if (v > 1) {
        for (j in 1:(v - 1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          model_string <- paste0(model_string, " + bs", v, j, "*", prev_response_trt, "s[i]")
        }
      }
      model_string <- paste0(model_string, "\n")

      # Marginal censoring model
      model_string <- paste0(model_string,
                             response_cen, "s[i] ~ dbern(cp", v, "s[i])\n",
                             "logit(cp", v, "s[i]) <- ts", v, "0")
      if (v > 1) {
        for (j in 1:(v - 1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          model_string <- paste0(model_string, " + ts", v, j, "*", prev_response_trt, "s[i]")
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
      }

      # Censoring priors
      for (p in 0:num_preds_cen) {
        model_string <- paste0(model_string, "s", v, p, " ~ dunif(-10, 10)\n")
      }

      # Marginal treatment priors
      model_string <- paste0(model_string, "bs", v, "0 ~ dunif(-10, 10)\n")
      if (v > 1) {
        for (j in 1:(v - 1)) {
          model_string <- paste0(model_string, "bs", v, j, " ~ dunif(-10, 10)\n")
        }
      }

      # Marginal censoring priors
      model_string <- paste0(model_string, "ts", v, "0 ~ dunif(-10, 10)\n")
      if (v > 1) {
        for (j in 1:(v - 1)) {
          model_string <- paste0(model_string, "ts", v, j, " ~ dunif(-10, 10)\n")
        }
      }
    }

    # Add the closing brace for the model block
    model_string <- paste0(model_string, "}\n")

    # Assume model_string is being built here...

    # If saving to file, write and return file path
    if (save_jags_model_file) {
      writeLines(model_string, output_file)
      return(output_file)
    } else {
      # Otherwise, return the model as a text connection for direct use in R2jags
      return(model_string)
    }
  }

  # extracting parameters list;
  jags_model_parameter <- function(trtmodel.list,
                                   cenmodel.list){
    all_parameters <- c()
    var_info_trt <- lapply(trtmodel.list, extract_variables_list)
    var_info_cen <- lapply(cenmodel.list, extract_variables_list)

    for (v in seq_along(var_info_trt)) {
      visit_trt <- var_info_trt[[v]]
      response_trt <- visit_trt$response
      predictors_trt <- visit_trt$predictors
      visit_cen <- var_info_cen[[v]]
      response_cen <- visit_cen$response
      predictors_cen <- visit_cen$predictors

      # Treatment model
      for (p in seq_along(predictors_trt)) {
        all_parameters <- c(all_parameters, sprintf("b%d%d", v, p))
      }

      # Censoring model
      for (p in seq_along(predictors_cen)) {
        all_parameters <- c(all_parameters, sprintf("s%d%d", v, p))
      }

      all_parameters <- c(all_parameters, sprintf("bs%d0", v))

      if (v > 1) {
        for (j in 1:(v - 1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          all_parameters <- c(all_parameters, sprintf("bs%d%d", v, j))
        }
      }

      all_parameters <- c(all_parameters, sprintf("ts%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          all_parameters <- c(all_parameters, sprintf("ts%d%d", v, j))
        }
      }
    }

    # Priors section
    for (v in seq_along(var_info_trt)) {
      num_preds_trt <- length(var_info_trt[[v]]$predictors)
      num_preds_cen <- length(var_info_cen[[v]]$predictors)

      # Treatment priors
      for (p in 0:num_preds_trt) {
        all_parameters <- c(all_parameters, sprintf("b%d%d", v, p))
      }

      # Censoring priors
      for (p in 0:num_preds_cen) {
        all_parameters <- c(all_parameters, sprintf("s%d%d", v, p))
      }

      # Marginal treatment priors
      all_parameters <- c(all_parameters, sprintf("bs%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          all_parameters <- c(all_parameters, sprintf("bs%d%d", v, j))
        }
      }

      # Marginal censoring priors
      all_parameters <- c(all_parameters, sprintf("ts%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          all_parameters <- c(all_parameters, sprintf("ts%d%d", v, j))
        }
      }
    }

    return(unique(all_parameters))

  }

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
  jags.params <- jags_model_parameter(trtmodel.list, cenmodel.list)

  # Check if parallel computing is requested
  # Handle JAGS model file creation based on parallel and save_jags_model_file flags
  if (parallel == TRUE) {
    # If parallel, always use files (either user-specified or temp files)
    if (save_jags_model_file== TRUE) {
      # User wants to save the model to a file
      jags.model.file <- write_jags_model(trtmodel.list, cenmodel.list,
                                          save_jags_model_file=TRUE,
                                          output_file=output_file)
      jags.model.files <- rep(jags.model.file, n.chains)  # Use same file for all chains
    } else {
      # Create temporary model files for each chain
      jags.model.files <- sapply(1:n.chains, function(i) tempfile(fileext = ".txt"))
      sapply(jags.model.files, function(f) {
        model_str <- write_jags_model(trtmodel.list, cenmodel.list,
                                      save_jags_model_file=FALSE,
                                      output_file=f)
        writeLines(model_str, f)  # Ensure model is written
      })
    }
  } else {
    # If not parallel, handle saving logic
    if (save_jags_model_file== TRUE) {
      jags.model.file <- write_jags_model(trtmodel.list, cenmodel.list,
                                          save_jags_model_file=TRUE,
                                          output_file=output_file)
    } else {
      jags.model.file <- textConnection(write_jags_model(trtmodel.list, cenmodel.list,
                                                         save_jags_model_file=FALSE))  # textConnection
    }
  }

  # Define seed for each chain
  if(!is.null(seed)){
    set.seed(seed)
  }
  new.seed.by.chain <- sample.int(2^30, size=n.chains)

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

    # jags.model.wd <- paste(getwd(), '/censoring_model.txt',sep='')

    posterior <- foreach::foreach(chain_idx=1:n.chains, .packages=c('R2jags'),
                         .combine='rbind') %dopar%{

                           jagsfit <- jags(data = jags.data,
                                           parameters.to.save = jags.params,
                                           model.file = jags.model.files[chain_idx],
                                           n.chains = 1,
                                           n.iter = n.iter,
                                           n.burnin = n.burnin,
                                           n.thin = n.thin,
                                           jags.seed = new.seed.by.chain[chain_idx])
                           # Combine MCMC output from multiple chains
                           out.mcmc <- as.mcmc(jagsfit)
                           return(do.call(rbind, lapply(out.mcmc, as.matrix)))

                         }

    parallel::stopCluster(cl)

    # Clean up temp files if they were used
    if (!save_jags_model_file) {
      file.remove(jags.model.files)
    }

  } else if (parallel == FALSE) {
    if (n.chains != 1) {
      stop("Non-parallel MCMC requires exactly 1 chain.")
    }
    jagsfit <- jags(data = jags.data,
                    parameters.to.save = jags.params,
                    model.file = jags.model.file,
                    n.chains = 1,
                    n.iter = n.iter,
                    n.burnin = n.burnin,
                    n.thin = n.thin,
                    jags.seed = new.seed.by.chain[1])

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
