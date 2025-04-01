#' Boxplots for Overall Comparisons by Continuous Variables.
#'
#' @description
#' This function creates a PDF file containing box plots for the continuous
#' variables in the provided data. If the number of columns in `data` exceeds
#' `bin.size`, the function splits the plots across multiple pages.
#'
#' @param data A matrix or data frame containing the raw data to be plotted.
#' @param pdf_title A string representing the name of the PDF file to
#' be created. If set to \code{NULL}, the box plots are displayed on the current
#' graphics device. Default is \code{NULL}.
#' @param bin_size An integer specifying the maximum number of box plots to
#' display on a single page.
#' @param y_lim An optional numeric vector defining the y-axis limits
#' for the plots.
#' @param scale An optional character string. If set to "log2",
#' numeric columns are log2-transformed.
#'
#'
#' @return A PDF file containing the box plots for the continuous variables.
#'
#' @examples
#' # Loading data
#' data.df <- ExampleData1
#' # Generate box plots for log2-transformed values to check for outliers:
#' cyt_bp(data.df[,-c(1:3)], pdf_title = NULL, scale = "log2")
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
cyt_bp <- function(data, pdf_title, bin_size = 25, y_lim = NULL, scale = NULL) {
  # Ensure data is a data frame
  data <- as.data.frame(data)

  # Apply log2 transformation if requested
  if (!is.null(scale) && scale == "log2") {
    numeric_cols <- sapply(data, is.numeric)
    for (col in names(data)[numeric_cols]) {
      # Replace non-positive values with NA to avoid issues with log transformation
      data[[col]][data[[col]] <= 0] <- NA
      data[[col]] <- log2(data[[col]])
    }
  }

  n_col <- ncol(data)
  if (n_col == 0) stop("No columns to plot in 'data'.")

  n_chunks <- ceiling(n_col / bin_size)
  plot_list <- vector("list", n_chunks)

  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * bin_size + 1
    end_idx <- min(i * bin_size, n_col)
    chunk_cols <- names(data)[start_idx:end_idx]

    melted <- reshape2::melt(data[, chunk_cols, drop = FALSE],
                             measure.vars = chunk_cols,
                             variable.name = "Variable",
                             value.name = "Value")

    p <- ggplot2::ggplot(melted, ggplot2::aes(x = Variable, y = Value)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(title = "Boxplots for Variables:",
                    x = "Variable", y = "Value") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

    if (!is.null(y_lim)) {
      p <- p + ggplot2::coord_cartesian(ylim = y_lim)
    }

    plot_list[[i]] <- p
  }

  # Write all plots to the PDF if a filename is provided; each plot goes on a new page.
  if (!is.null(pdf_title)) {
    pdf(file = pdf_title, width = 7, height = 5)
    on.exit(dev.off(), add = TRUE)
  }

  for (p in plot_list) {
    print(p)
  }

  invisible(NULL)
}
