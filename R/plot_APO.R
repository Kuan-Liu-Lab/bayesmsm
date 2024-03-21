#' Plot Average Potential Outcomes (APO)
#'
#' This function plots the density of APO for a specified effect type from bootstrap simulation results.
#'
#' @param input A data frame or model object containing 'bootdata', which include 'effect_comparator' and 'effect_reference' columns.
#' @param effect_type A character string specifying which effect to plot: 'effect_comparator' or 'effect_reference'.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A density plot showing the distribution of the specified average potential outcome (reference or comparison).
#' @export
#'
#' @examples
#' plot_APO(model$bootdata, effect_type = "effect_comparator")
#' plot_APO(model, effect_type = "effect_reference")
#'
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
  density_effect <- density(effect[[1]])

  # Define titles and colors based on effect_type
  titles <- c(effect_comparator = "Comparator Level", effect_reference = "Reference Level")
  colors <- c(effect_comparator = "blue", effect_reference = "red")

  # Calculate mean and median
  mean_effect <- mean(effect[[1]])
  median_effect <- median(effect[[1]])

  # Plotting
  plot(density_effect, main = paste("Average Potential Outcome (APO) of", titles[effect_type]), xlab = "Effect", ylab = "Density", col = colors[effect_type], lwd = 2, ...)

  # Add a vertical line for the mean and median
  abline(v = mean_effect, col = "darkgrey", lty = 3)
  abline(v = median_effect, col = "darkgrey", lty = 4)

  # Legend with mean and median
  legend_text <- c(titles[effect_type],
                   paste("Mean:", round(mean_effect, 3)),
                   paste("Median:", round(median_effect, 3)))

  legend("topright", legend = legend_text, col = c(colors[effect_type], "darkgrey", "darkgrey"),
         lwd = 2, lty = c(1, 3, 4))
}
