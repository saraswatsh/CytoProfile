#' Boxplot Function Enhanced for Specific Group Comparisons. `r lifecycle::badge("deprecated")`
#'
#' This function generates a PDF file containing boxplots for each combination
#' of numeric and factor variables in the provided data. It first converts
#' any character columns to factors and checks that the data contains at
#' least one numeric and one factor column. If the scale argument is set to
#' "log2", all numeric columns are log2-transformed. The function then
#' creates boxplots using ggplot2 for each numeric variable grouped by
#' each factor variable.
#'
#' @param data A matrix or data frame of raw data.
#' @param pdf_title A string representing the title
#' (and filename) of the PDF file. If \code{NULL}, the boxplots are displayed on the
#' current graphics device. Defaults to \code{NULL}.
#' @param scale Transformation option for continuous variables.
#' Options are NULL (default) and "log2". When set to "log2",
#' numeric columns are transformed using the log2 function.
#' @param y_lim An optional numeric vector defining the y-axis limits
#' for the plots.
#'
#' @return A PDF file containing the boxplots.
#'
#' @author Shubh Saraswat
#'
#' @examples
#' # Loading data
#' data_df <- ExampleData1[, -c(3, 5:28)]
#' data_df <- dplyr::filter(data_df, Group == "T2D", Treatment == "Unstimulated")
#' cyt_bp2(data_df, pdf_title = NULL, scale = "log2")
#'
#' @export
#' @import ggplot2

cyt_bp2 <- function(data, pdf_title, scale = NULL, y_lim = NULL) {
  lifecycle::deprecate_warn(
    "0.2.3", # version when deprecation begins
    "CytoProfile::cyt_bp2()",
    "CytoProfile::cyt_bp()"
  )
  # Convert any character variables to factors
  cat_vars <- sapply(data, is.character)
  if (any(cat_vars)) {
    data[cat_vars] <- lapply(data[cat_vars], as.factor)
  }

  # Identify numeric and factor columns
  num_cols <- sapply(data, is.numeric)
  fac_cols <- sapply(data, is.factor)

  # Ensure there is at least one numeric and one factor column
  if (!any(num_cols)) {
    stop("Data must contain at least one numeric column")
  }
  if (!any(fac_cols)) {
    stop("Data must contain at least one factor column")
  }

  # Apply log2 transformation if scale is "log2"
  if (!is.null(scale) && scale == "log2") {
    data[num_cols] <- lapply(data[num_cols], function(x) log2(x))
  }

  # Open PDF device if pdf_title is provided and ensure it's closed when function exits
  if (!is.null(pdf_title)) {
    pdf(file = pdf_title)
    on.exit(dev.off(), add = TRUE)
  }

  # Loop through each numeric column
  for (num_col in names(data)[num_cols]) {
    # Loop through each factor column
    for (fac_col in names(data)[fac_cols]) {
      # Create the boxplot
      p <- ggplot2::ggplot(
        data,
        aes_string(x = fac_col, y = num_col, fill = fac_col)
      ) +
        ggplot2::geom_boxplot(alpha = 0.5) +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
        ggplot2::labs(
          title = paste0(num_col, " by ", fac_col),
          x = fac_col,
          y = paste0("Values of ", num_col)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"),
          legend.background = element_rect(fill = "white", colour = "white"),
          axis.title = element_text(color = "black", size = 12, face = "bold"),
          legend.title = element_text(
            color = "black",
            size = 10,
            face = "bold"
          ),
          legend.text = element_text(color = "black")
        ) +
        ggplot2::guides(fill = guide_legend(title = NULL))

      if (!is.null(y_lim)) {
        p <- p + ggplot2::ylim(y_lim)
      }

      print(p)
    }
  }
}
