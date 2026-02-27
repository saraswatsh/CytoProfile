#' Violin Plots for Continuous Variables with Optional Grouping
#'
#' @description
#' `cyt_violin` produces violin plots for each numeric variable in
#' `data`, optionally grouped by one or more categorical variables.
#' When grouping is not specified, the function behaves similarly to
#' `cyt_bp` but uses violins instead of boxplots and supports
#' pagination via the `bin_size` argument.  When grouping is
#' provided, a separate violin is drawn for each level (or
#' interaction of levels) of the grouping variables.  Users may
#' optionally overlay boxplots within each violin to visualize the
#' median and interquartile range.
#'
#' @param data A matrix or data frame containing numeric and
#'   categorical variables.
#' @param output_file Optional string specifying the name of the file
#'   to be created.  When `NULL` (default), plots are drawn on
#'   the current graphics device. Ensure that the file
#'   extension matches the desired format (e.g., ".pdf" for PDF output
#'   or ".png" for PNG output or .tiff for TIFF output).
#' @param group_by Optional character vector specifying one or more
#'   columns to use for grouping.  If `NULL` (default) no grouping is
#'   applied.
#' @param bin_size Integer.  Maximum number of violins per page when
#'   grouping is not used.  Default is 25, mirroring `cyt_bp`.
#' @param y_lim Optional numeric vector giving y‑axis limits for the
#'   plots.  Applies to all plots.
#' @param scale Character specifying a transformation for numeric
#'   variables.  Accepts "none", "log2", "log10",
#'   "zscore", or "custom".  When "custom", supply a
#'   function via `custom_fn`.
#' @param custom_fn A user supplied function to transform numeric
#'   columns when `scale = "custom"`.
#' @param boxplot_overlay Logical.  When `TRUE`, a narrow boxplot is
#'   drawn inside each violin to summarize the median and quartiles.
#'   Default is `FALSE`.
#' @return Invisibly returns a list of `ggplot` objects.  When
#'   `output_file` is provided, plots are written to the PDF file.
#' @examples
#' # Violin plots without grouping
#' cyt_violin(ExampleData1[, -c(1:3)], output_file = NULL, scale = "zscore")
#' # Violin plots grouped by Group with boxplot overlay
#' cyt_violin(ExampleData1[, -c(3,5:28)], group_by = "Group",
#'                         boxplot_overlay = TRUE, scale = "log2")
#' @author Shubh Saraswat
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
cyt_violin <- function(
  data,
  output_file = NULL,
  group_by = NULL,
  bin_size = 25,
  y_lim = NULL,
  scale = c("none", "log2", "log10", "zscore", "custom"),
  custom_fn = NULL,
  boxplot_overlay = FALSE
) {
  names(data) <- make.names(names(data), unique = TRUE)
  scale <- match.arg(scale)
  df <- as.data.frame(data)
  # Validate grouping columns
  if (!is.null(group_by)) {
    if (!all(group_by %in% names(df))) {
      stop("All group_by columns must exist in data.")
    }
  }
  # Identify numeric variables excluding grouping columns
  numeric_cols <- names(df)[sapply(df, is.numeric) & !(names(df) %in% group_by)]
  if (length(numeric_cols) == 0) {
    stop("No numeric columns to plot.")
  }
  # Apply scaling
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
    # Interaction term for multiple grouping variables
    if (length(group_by) > 1) {
      df$._group <- interaction(df[, group_by], sep = ":")
      grp_col <- "._group"
    } else {
      grp_col <- group_by
    }
  } else {
    grp_col <- NULL
  }
  plot_list <- list()
  # Grouping case
  if (!is.null(grp_col)) {
    for (var in numeric_cols) {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = .data[[grp_col]],
          y = .data[[var]],
          fill = .data[[grp_col]]
        )
      ) +
        ggplot2::geom_violin(trim = FALSE, alpha = 0.6)
      if (boxplot_overlay) {
        p <- p +
          ggplot2::geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA)
      }
      p <- p +
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
    if (!is.null(output_file)) {
      ext <- tolower(tools::file_ext(output_file))
      if (ext == "pdf") {
        pdf(file = output_file, width = 7, height = 5)
        on.exit(dev.off(), add = TRUE)
        for (p in plot_list) {
          print(p)
        }
      } else {
        # For PNG, SVG, etc., save each page as a separate file
        for (i in seq_along(plot_list)) {
          fname <- if (length(plot_list) > 1) {
            sub(paste0("\\.", ext, "$"), paste0("_", i, ".", ext), output_file)
          } else {
            output_file
          }
          ggplot2::ggsave(fname, plot_list[[i]], width = 7, height = 5)
        }
      }
    } else {
      for (p in plot_list) {
        print(p)
      }
    }
    return(invisible(plot_list))
  }
  # No grouping: paginate violins across numeric variables
  n_col <- length(numeric_cols)
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
      ggplot2::geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.7)
    if (boxplot_overlay) {
      p <- p +
        ggplot2::geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA)
    }
    y_label <- if (scale == "none") {
      "Value"
    } else {
      paste0("Value (", scale, "-transformed)")
    }
    p <- p +
      ggplot2::labs(
        title = "Violin Plots for Variables:",
        x = "Variable",
        y = y_label
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust = 1)
      ) +
      ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 1))
    if (!is.null(y_lim)) {
      p <- p + ggplot2::coord_cartesian(ylim = y_lim)
    }
    plot_list[[paste0("page", i)]] <- p
  }
  if (!is.null(output_file)) {
    ext <- tolower(tools::file_ext(output_file))
    if (ext == "pdf") {
      pdf(file = output_file, width = 7, height = 5)
      on.exit(dev.off(), add = TRUE)
      for (p in plot_list) {
        print(p)
      }
    } else {
      # For PNG, SVG, etc., save each page as a separate file
      for (i in seq_along(plot_list)) {
        fname <- if (length(plot_list) > 1) {
          sub(paste0("\\.", ext, "$"), paste0("_", i, ".", ext), output_file)
        } else {
          output_file
        }
        ggplot2::ggsave(fname, plot_list[[i]], width = 7, height = 5)
      }
    }
  } else {
    for (p in plot_list) {
      print(p)
    }
  }
  invisible(plot_list)
}
