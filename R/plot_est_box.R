#' Error bar plots for causal treatment effects
#'
#' This function plots the point estimates and 95\% credible intervals of ATE and APO from bayesmsm output.
#'
#' @param input A data frame or model object containing bootstrap results.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot object presenting error bar plot of the mean effects and their 95\% credible intervals for comparator level, reference level, and ATE.
#' @importFrom stats density quantile
#' @importFrom grDevices rgb
#' @importFrom graphics abline arrows axis legend mtext par polygon text
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
#' plot_est_box(model)
plot_est_box <- function(input, ...) {
  # Extract bootdata from the model or use the data frame directly
  bootdata <- if (is.data.frame(input)) {
    input
  } else if ("bootdata" %in% names(input)) {
    input$bootdata
  } else {
    stop("Input must be a data frame or a model object containing 'bootdata'.")
  }

  # Identify if the family is binomial to include RR and OR
  is_binomial <- "RR" %in% names(bootdata) && "OR" %in% names(bootdata)

  # Validate bootdata
  required_columns <- c("effect_comparator", "effect_reference", "RD")
  if (is_binomial) {
    required_columns <- c(required_columns, "RR", "OR")
  }
  if (!all(required_columns %in% names(bootdata))) {
    stop("bootdata must contain 'effect_comparator', 'effect_reference', and 'RD' columns.")
  }

  # Calculate means and confidence intervals
  means <- sapply(bootdata[required_columns], mean)
  lowerbd <- sapply(bootdata[required_columns], function(x) quantile(x, probs = 0.025))
  upperbd <- sapply(bootdata[required_columns], function(x) quantile(x, probs = 0.975))

  # Create data frame for ggplot
  plot_data <- data.frame(
    Treatment = factor(names(means), levels = names(means)),
    Mean = means,
    LowerCI = lowerbd,
    UpperCI = upperbd
  )

  # ggplot2 visualization
  ggplot2::ggplot(plot_data, ggplot2::aes(x = Treatment, y = Mean)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "blue") +
    ggplot2::labs(title = "Treatment Effect Estimates", x = "Treatment Level", y = "Effect") +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(ggplot2::aes(y = UpperCI + 0.05, label = paste0("Mean: ", round(Mean, 2))), vjust = -0.5) +
    ggplot2::geom_text(ggplot2::aes(y = UpperCI + 0.50, label = paste0("95% CI: [", round(LowerCI, 2), ", ", round(UpperCI, 2), "]")), vjust = -0.5) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::expand_limits(y = max(plot_data$UpperCI) + 0.15)
}
