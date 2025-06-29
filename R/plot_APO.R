#' Plot Average Potential Outcomes (APO)
#'
#' This function plots the density of APO for a specified effect type from bayesmsm output.
#'
#' @param input A data frame or model object containing bootstrap results.
#' @param effect_type A character string specifying which effect to plot (e.g., comparator or reference treatment sequences).
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot object representing density plot showing the distribution of the specified average potential outcome (reference or comparison).
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
#' plot_APO(model$bootdata, effect_type = "effect_comparator")
#' plot_APO(model, effect_type = "effect_reference")
plot_APO <- function(input, effect_type, ...) {
  # Validate input
  if ("bootdata" %in% names(input)) {
    bootdata <- input$bootdata
  } else if (is.data.frame(input)) {
    bootdata <- input
  } else {
    stop("Input must be a data frame or a model object containing a 'bootdata' data frame.")
  }
  if (!is.data.frame(bootdata) || !("effect_comparator" %in% names(bootdata)) || !("effect_reference" %in% names(bootdata))) {
    stop("bootdata must be a data frame containing 'effect_comparator' and 'effect_reference' columns.")
  }
  if (!is.character(effect_type) || length(effect_type) != 1) {
    stop("effect_type must be a single character string specifying the effect to plot.")
  }
  if (!effect_type %in% c("effect_comparator", "effect_reference")) {
    stop("effect_type must be either 'effect_comparator' or 'effect_reference'.")
  }

  # Extract the relevant column
  effect <- bootdata[, effect_type, drop = FALSE]

  # Calculate density
  density_effect <- stats::density(effect[[1]])

  # Define titles and colors based on effect_type
  titles <- c(effect_comparator = "Comparator Level", effect_reference = "Reference Level")
  colors <- c(effect_comparator = "blue", effect_reference = "red")

  # Calculate mean
  mean_effect <- mean(effect[[1]])

  # Calculate CI
  ci <- stats::quantile(effect[[1]], probs = c(0.025, 0.975))

  # Create data frame for ggplot
  density_data <- data.frame(x = density_effect$x, y = density_effect$y)
  ci_data <- data.frame(x = c(ci[1], ci[2]), y = c(0, 0))

  # ggplot2 visualization
  ggplot2::ggplot() +
    ggplot2::geom_line(data = density_data, ggplot2::aes(x = x, y = y), color = colors[effect_type], linewidth = 1) +
    ggplot2::geom_ribbon(data = density_data, ggplot2::aes(x = x, ymin = 0, ymax = y), fill = "lightblue", alpha = 0.3) +
    ggplot2::geom_vline(xintercept = mean_effect, color = "purple", linetype = "dashed", linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = ci, color = "darkgreen", linetype = "dotted", linewidth = 1.2) +
    ggplot2::labs(title = paste("Average Potential Outcome (APO) of", titles[effect_type]), x = "Effect", y = "Density") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = mean_effect, y = max(density_data$y) * 0.9,
                      label = paste("Mean:", round(mean_effect, 3)),
                      color = "purple", angle = 90, vjust = -0.5) +
    ggplot2::annotate("text", x = ci[1], y = max(density_data$y) * 0.8,
                      label = paste("95% CI Lower:", round(ci[1], 3)),
                      color = "darkgreen", angle = 90, vjust = -0.5) +
    ggplot2::annotate("text", x = ci[2], y = max(density_data$y) * 0.8,
                      label = paste("95% CI Upper:", round(ci[2], 3)),
                      color = "darkgreen", angle = 90, vjust = -0.5)
}
