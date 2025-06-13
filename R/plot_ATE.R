#' Plot Average Treatment Effect Density from bayesmsm output
#'
#' This function plots the density of ATE from bayesmsm output.
#'
#' @param input A model object, data frame or vector containing the bootstrap estimates of ATE.
#' @param ATE define causal estimand of interest from RD, OR, RR.
#' @param col_density Color for the density plot (default is "blue").
#' @param fill_density Fill color for the density plot (default is "lightblue").
#' @param main Title of the plot (default is "Density of ATE Estimates").
#' @param xlab X-axis label (default is "ATE").
#' @param ylab Y-axis label (default is "Density").
#' @param xlim Limits for the x-axis (default is NULL).
#' @param ylim Limits for the y-axis (default is NULL).
#' @param ... Additional graphical parameters passed to the plot function.
#' @return A ggplot object representing the density plot for the posterior predictive distribution of the Average Treatment Effect (ATE).
#' @import ggplot2
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
#' plot_ATE(model)
plot_ATE <- function(input,
                     ATE = "RD",
                     col_density = "blue",
                     fill_density = "lightblue",
                     main = "Posterior Predictive Distribution of Average Treatment Effect",
                     xlab = "ATE", ylab = "Posterior Predictive Distribution",
                     xlim = NULL, ylim = NULL, ...) {
  # Check if input is either a data frame or part of a model object
  if (is.list(input) && "bootdata" %in% names(input)) {
    # If input is a list and has bootdata, check for ATE column within bootdata
    if (ATE %in% names(input$bootdata)) {
      ate_values <- unlist(input$bootdata[ATE])
    } else {
      stop("bootdata within the model object must have an 'ATE' column.")
    }
  }

  # Calculate density and CI
  ate_density <- stats::density(ate_values)
  ci <- stats::quantile(ate_values, probs = c(0.025, 0.975))

  # Create data frame for ggplot
  density_data <- data.frame(x = ate_density$x, y = ate_density$y)
  ci_data <- data.frame(x = c(ci[1], ci[2]), y = c(0, 0))

  # ggplot2 visualization
  ggplot2::ggplot() +
    ggplot2::geom_line(data = density_data, ggplot2::aes(x = x, y = y), color = col_density, linewidth = 1) +
    ggplot2::geom_ribbon(data = density_data, ggplot2::aes(x = x, ymin = 0, ymax = y), fill = fill_density, alpha = 0.3) +
    ggplot2::geom_vline(xintercept = mean(ate_values), color = "purple", linetype = "dashed", linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = ci, color = "darkgreen", linetype = "dotted", linewidth = 1.2) +
    ggplot2::labs(title = main, x = xlab, y = ylab) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = mean(ate_values), y = max(density_data$y) * 0.9,
                      label = paste("Mean:", round(mean(ate_values), 3)),
                      color = "purple", angle = 90, vjust = -0.5) +
    ggplot2::annotate("text", x = ci[1], y = max(density_data$y) * 0.8,
                      label = paste("95% CI Lower:", round(ci[1], 3)),
                      color = "darkgreen", angle = 90, vjust = -0.5) +
    ggplot2::annotate("text", x = ci[2], y = max(density_data$y) * 0.8,
                      label = paste("95% CI Upper:", round(ci[2], 3)),
                      color = "darkgreen", angle = 90, vjust = -0.5)
}
