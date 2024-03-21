#' Error bar plots for treatment effects
#'
#' @param input A data frame or model object containing 'bootdata', which include 'effect_comparator', 'effect_reference', and 'ATE' columns.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return An error bar plot of the mean effects and their 95\% confidence intervals for comparator level, reference level, and ATE.
#' @export
#'
#' @examples
#' plot_est_box(model$bootdata) # without reference & comparator information below labels
#' plot_est_box(model) # with reference & comparator information below labels
#'
plot_est_box <- function(input, ...) {
  # Extract bootdata from the model or use the data frame directly
  bootdata <- if (is.data.frame(input)) {
    input
  } else if ("bootdata" %in% names(input)) {
    input$bootdata
  } else {
    stop("Input must be a data frame or a model object containing 'bootdata'.")
  }

  # Validate bootdata
  required_columns <- c("effect_comparator", "effect_reference", "ATE")
  if (!all(required_columns %in% names(bootdata))) {
    stop("bootdata must contain 'effect_comparator', 'effect_reference', and 'ATE' columns.")
  }

  # Adjust margins if necessary
  par(mar = c(5, 4, 4, 3) + 0.1) # Adjust the last value if text is plotted outside; bottom, left, top, and right margins

  # Calculate means and standard deviations
  means <- sapply(bootdata[required_columns], mean)
  # sds <- sapply(bootdata[required_columns], sd)
  # ses <- sds / sqrt(nrow(bootdata))
  lowerbd <- sapply(bootdata[required_columns], function(x) quantile(x, probs = 0.025))
  upperbd <- sapply(bootdata[required_columns], function(x) quantile(x, probs = 0.975))

  # Define the position for each point
  position <- 1:length(means)

  # Define some offsets for text placement
  text_offset <- (max(upperbd) - min(lowerbd)) * 0.05

  # Plotting
  plot(position, means, ylim = range(lowerbd - text_offset, upperbd + text_offset), pch = 19, xaxt = "n", # round down vs round up;
       xlab = "Treatment Level", ylab = "Effect", main = "Treatment Effect Estimates", ...)
  axis(1, at = position, labels = c("Comparator Level", "Reference Level", "ATE"))

  # Error bars
  arrows(position, lowerbd, position, upperbd, angle = 90, code = 3, length = 0.1, ...)

  # Adding text for means and CIs
  # text(position, lowerbd - 0.1, labels = paste("Mean:", round(means, 2), "\n95% CI:", round(lowerbd, 2), "-", round(upperbd, 2)), adj = c(0,1))

  # for (i in seq_along(means)) {
  #   text(position[i], upperbd[i] + text_offset, labels = paste("Mean:", round(means[i], 2), "\n95% CI:", round(lowerbd[i], 2), "-", round(upperbd[i], 2)), cex = 0.8, pos = 3)
  # }

  # Adjust the y-position offset for clarity
  offset <- 0.1
  # Adding text for the first and second items with left adjustment
  text(position[1:2], lowerbd[1:2] - offset, labels = paste("Mean:", round(means[1:2], 2), "\n95% CI: [", round(lowerbd[1:2], 2), ", ", round(upperbd[1:2], 2), "]"), adj = c(0,1), ...)
  # Adding text for the third item with right adjustment
  text(position[3], upperbd[3] + offset, labels = paste("Mean:", round(means[3], 2), "\n95% CI: [", round(lowerbd[3], 2), ", ", round(upperbd[3], 2), "]"), adj = c(1,0), ...)

  # Check if the input is a model and extract treatment sequences if they exist
  has_treatment_info <- "reference" %in% names(input) && "comparator" %in% names(input)

  # Conditional treatment sequence information below x-axis labels
  if (has_treatment_info) {
    # mtext(paste("(", paste(input$reference, collapse = ", "), ")", sep = ""), side = 1, at = position[2], line = 2.5, cex = 0.7)
    # mtext(paste("(", paste(input$comparator, collapse = ", "), ")", sep = ""), side = 1, at = position[1], line = 2.5, cex = 0.7)
    mtext(paste("(", paste(input$reference, collapse = ", "), ")", sep = ""), side = 1, at = position[2], line = 2)
    mtext(paste("(", paste(input$comparator, collapse = ", "), ")", sep = ""), side = 1, at = position[1], line = 2)
  }
}
