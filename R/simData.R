#' Generate synthetic longitudinal data with optional right-censoring
#'
#' This function simulates repeated measurements of normally-distributed covariates,
#' binary treatments, and an end-of-study outcome for longitudinal causal analyses.
#' When `right_censor = TRUE`, a right-censoring indicator `Cj` is generated at each visit:
#' if `Cj = 1`, all subsequent `L`, `A`, and `Y` values are set to `NA`.
#'
#' @param n Integer. Sample size.
#' @param n_visits Integer. Number of visits (including baseline as visit 1).
#' @param covariate_counts Integer vector of length `n_visits`. Number of covariates per visit (default: rep(2, n_visits)).
#' @param amodel List of length `n_visits`. Each element is a named numeric vector of coefficients
#'   for the logistic model of treatment `Aj` on covariates (and `A_prev` for j > 1).
#' @param ymodel Named numeric vector. Coefficients for the end-of-study outcome model.
#'   If `y_type = "binary"`, a logistic model is used; if `"continuous"`, a linear model with Gaussian noise.
#' @param y_type Character. One of "binary" or "continuous".
#' @param right_censor Logical. If TRUE, generates `Cj` using `cmodel` at each visit.
#' @param cmodel List of length `n_visits`. Named numeric vectors for logistic censoring models
#'   at each visit, regressing `Cj` on covariates and current `Aj`.
#' @param seed Integer. Optional random seed.
#' @importFrom stats rnorm rbinom model.matrix
#' @return A `data.frame` with columns `Lk_j`, `Aj`, optional `Cj`, and `Y`.
#' @export
simData <- function(n,
                    n_visits,
                    covariate_counts = rep(2, n_visits),
                    amodel,
                    ymodel,
                    y_type = c("binary", "continuous"),
                    right_censor = FALSE,
                    cmodel = NULL,
                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y_type <- match.arg(y_type)
  expit <- function(x) exp(x) / (1 + exp(x))

  # track censoring
  censored <- rep(FALSE, n)
  visit_data <- vector("list", n_visits)
  Aj_prev <- NULL

  for (j in seq_len(n_visits)) {
    # simulate covariates
    p_cov <- covariate_counts[j]
    L <- as.data.frame(matrix(rnorm(n * p_cov), nrow = n,
                              dimnames = list(NULL, paste0("L", seq_len(p_cov), "_", j))))
    L[censored, ] <- NA

    # treatment model
    dfA <- L
    if (j > 1) dfA$A_prev <- Aj_prev
    X_A <- model.matrix(~ ., data = dfA)
    coefA <- rep(0, ncol(X_A)); names(coefA) <- colnames(X_A)
    coefA[names(amodel[[j]])] <- amodel[[j]]
    pA <- expit(X_A %*% coefA)
    A_j <- rbinom(n, 1, pA)
    A_j[censored] <- NA

    # censoring model (if requested)
    if (right_censor) {
      if (is.null(cmodel) || length(cmodel) < j) {
        stop("Provide cmodel list of length >= n_visits when right_censor=TRUE.")
      }
      dfC <- dfA; dfC$A <- A_j
      X_C <- model.matrix(~ ., data = dfC)
      coefC <- rep(0, ncol(X_C)); names(coefC) <- colnames(X_C)
      coefC[names(cmodel[[j]])] <- cmodel[[j]]
      pC <- expit(X_C %*% coefC)
      C_j <- rbinom(n, 1, pC)
      C_j[censored] <- NA
      censored <- censored | (C_j == 1)
    }

    # store
    df_j <- L
    df_j[[paste0("A", j)]] <- A_j
    if (right_censor) df_j[[paste0("C", j)]] <- C_j
    visit_data[[j]] <- df_j
    Aj_prev <- A_j
  }

  # combine visits
  df <- do.call(cbind, visit_data)

  # outcome model
  X_Y <- model.matrix(~ ., data = df)
  coefY <- rep(0, ncol(X_Y)); names(coefY) <- colnames(X_Y)
  coefY[names(ymodel)] <- ymodel
  lpY <- X_Y %*% coefY
  if (y_type == "binary") {
    Y <- rbinom(n, 1, expit(lpY))
  } else {
    Y <- lpY + rnorm(n)
  }
  Y[censored] <- NA
  df$Y <- Y
  return(df)
}
