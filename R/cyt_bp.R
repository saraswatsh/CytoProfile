#' Boxplots for Continuous Variables with Optional Grouping
#'
#' @description
#' This function generates boxplots for numeric variables in a data
#' frame or matrix.  It supports optional grouping by one or more
#' categorical variables.  Numeric variables can be scaled using
#' various transformations before plotting.  When grouping is not
#' used, boxplots are arranged in pages with a specified maximum
#' number of plots per page.  Plots can be saved to a PDF file or
#' displayed on the current graphics device.
#'
#' @param data A matrix or data frame containing numeric and
#'   categorical variables.
#' @param pdf_title Optional string specifying the name of the PDF
#'   file to be created.  When `NULL` (default), plots are drawn on
#'   the current graphics device.
#' @param group_by Optional character vector specifying one or more
#'   columns to use for grouping.  If `NULL` (default) no grouping is
#'   applied.
#' @param bin_size Integer.  Maximum number of boxplots per page when
#'   grouping is not used.  Default is 25, as in the original
#'   `cyt_bp`.
#' @param y_lim Optional numeric vector giving y‑axis limits for the
#'   plots.  Applies to all plots.
#' @param scale Character specifying a transformation for numeric
#'   variables.  Accepts `"none"`, `"log2"`, `"log10"`,
#'   `"zscore"`, or `"custom"`.  When `"custom"`, supply a
#'   function via `custom_fn`.
#' @param custom_fn A user supplied function to transform numeric
#'   columns when `scale = "custom"`.
#' @return Invisibly returns a list of `ggplot` objects.  When
#'   `pdf_title` is provided, plots are written to the PDF file.
#' @examples
#' data("ExampleData1")
#' # Boxplots without grouping
#' cyt_bp(ExampleData1[, -c(1:3)], pdf_title = NULL, scale = "log2")
#' # Boxplots grouped by Group
#' cyt_bp(ExampleData1[, -c(3,5:28)], group_by = "Group", scale = "zscore")
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
cyt_bp <- function(
  data,
  pdf_title = NULL,
  group_by = NULL,
  bin_size = 25,
  y_lim = NULL,
  scale = c("none", "log2", "log10", "zscore", "custom"),
  custom_fn = NULL
) {
  scale <- match.arg(scale)
  df <- as.data.frame(data)
  # Identify numeric variables (exclude grouping columns)
  if (!is.null(group_by)) {
    if (!all(group_by %in% names(df))) {
      stop("All group_by columns must exist in data.")
    }
  }
  numeric_cols <- names(df)[sapply(df, is.numeric) & !(names(df) %in% group_by)]
  if (length(numeric_cols) == 0) {
    stop("No numeric columns to plot.")
  }
  # Apply scaling to numeric variables
  if (scale != "none") {
    df <- apply_scale(
      df,
      columns = numeric_cols,
      scale = scale,
      custom_fn = custom_fn
    )
  }
  # Convert character grouping columns to factors
  if (!is.null(group_by)) {
    for (g in group_by) {
      if (is.character(df[[g]])) {
        df[[g]] <- factor(df[[g]])
      }
    }
    # If multiple grouping columns, create interaction term
    if (length(group_by) > 1) {
      df$._group <- interaction(df[, group_by], sep = ":")
      grp_col <- "._group"
    } else {
      grp_col <- group_by
    }
  } else {
    grp_col <- NULL
  }
  # List to store plots
  plot_list <- list()
  # Grouping case
  if (!is.null(grp_col)) {
    # Loop through each numeric variable and create grouped boxplot
    for (var in numeric_cols) {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = .data[[grp_col]],
          y = .data[[var]],
          fill = .data[[grp_col]]
        )
      ) +
        ggplot2::geom_boxplot(alpha = 0.6) +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) +
        ggplot2::labs(
          title = paste0(var, " by ", paste(group_by, collapse = ":")),
          x = paste(group_by, collapse = ":"),
          y = var
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
      if (!is.null(y_lim)) {
        p <- p + ggplot2::coord_cartesian(ylim = y_lim)
      }
      plot_list[[var]] <- p
    }
    # Write to PDF if requested
    if (!is.null(pdf_title)) {
      pdf(file = pdf_title, width = 7, height = 5)
      on.exit(dev.off(), add = TRUE)
      for (p in plot_list) {
        print(p)
      }
    } else {
      # Print each plot to the current graphics device
      for (p in plot_list) {
        print(p) # displays plots on screen
      }
    }
    return(invisible(plot_list))
  }
  # No grouping: original cyt_bp logic
  n_col <- length(numeric_cols)
  if (n_col == 0) {
    stop("No numeric columns to plot in 'data'.")
  }
  n_chunks <- ceiling(n_col / bin_size)
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * bin_size + 1
    end_idx <- min(i * bin_size, n_col)
    chunk_cols <- numeric_cols[start_idx:end_idx]
    melted <- reshape2::melt(
      df[, chunk_cols, drop = FALSE],
      measure.vars = chunk_cols,
      variable.name = "Variable",
      value.name = "Value"
    )
    p <- ggplot2::ggplot(melted, ggplot2::aes(x = Variable, y = Value)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(
        title = "Boxplots for Variables:",
        x = "Variable",
        y = "Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
    if (!is.null(y_lim)) {
      p <- p + ggplot2::coord_cartesian(ylim = y_lim)
    }
    plot_list[[paste0("page", i)]] <- p
  }
  if (!is.null(pdf_title)) {
    pdf(file = pdf_title, width = 7, height = 5)
    on.exit(dev.off(), add = TRUE)
    for (p in plot_list) {
      print(p)
    }
  } else {
    for (p in plot_list) {
      print(p)
    }
  }
  invisible(plot_list)
}
