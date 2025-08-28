#' Distribution of the Data Set Shown by Skewness and Kurtosis.
#'
#' @description This function computes summary statistics --- including sample
#'   size, mean, standard error, skewness, and kurtosis --- for each numeric
#'   measurement column in a data set. If grouping columns are provided via
#'   `group_cols`, the function computes the metrics separately for each group
#'   defined by the combination of these columns (using the first element as
#'   the treatment variable and the second as the grouping variable, or
#'   the same column for both if only one is given). When no grouping columns
#'   are provided, the entire data set is treated as a single group ("Overall").
#'   A log2 transformation (using a cutoff equal to one-tenth of the smallest
#'   positive value in the data) is applied to generate alternative metrics.
#'   Histograms showing the distribution of skewness and kurtosis for both raw
#'   and log2-transformed data are then generated and saved to a PDF if a file
#'   name is provided.
#'
#' @param data A matrix or data frame containing the raw data. If
#'   `group_cols` is specified, the columns with names in `group_cols` are
#'   treated as grouping variables and all other columns are assumed to be
#'   numeric measurement variables.
#' @param group_cols A character vector specifying the names of the grouping
#'   columns. When provided, the first element is treated as the treatment
#'   variable and the second as the group variable. If not provided, the
#'   entire data set is treated as one group.
#' @param pdf_title A character string specifying the file name for the PDF file in
#'   which the histograms will be saved. If \code{NULL}, the histograms are
#'   displayed on the current graphics device. Default is \code{NULL}.
#' @param print_res_raw Logical. If \code{TRUE}, the function returns and prints
#'   the computed summary statistics for the raw data. Default is
#'   \code{FALSE}.
#' @param print_res_log Logical. If \code{TRUE}, the function returns and prints
#'   the computed summary statistics for the log2-transformed data. Default is
#'   \code{FALSE}.
#'
#' @return The function generates histograms of skewness and kurtosis for both
#'   raw and log2-transformed data. Additionally, if either
#'   \code{printResRaw} and/or \code{printResLog} is \code{TRUE}, the function
#'   returns the corresponding summary statistics as a data frame or a list of
#'   data frames.
#' @author Xiaohua Douglas Zhang and Shubh Saraswat
#'
#' @details A cutoff is computed as one-tenth of the minimum positive value
#'   among all numeric measurement columns to avoid taking logarithms of zero.
#'   When grouping columns are provided, the function loops over unique
#'   grouping columns and computes the metrics for each
#'   measurement column within each subgroup. Without grouping columns, the
#'   entire data set is analyzed as one overall group.
#'
#' @examples
#' # Example with grouping columns (e.g., "Group" and "Treatment")
#' data(ExampleData1)
#' cyt_skku(ExampleData1[, -c(2:3)], pdf_title = NULL,
#'   group_cols = c("Group")
#' )
#'
#' # Example without grouping columns (analyzes the entire data set)
#' cyt_skku(ExampleData1[, -c(1:3)], pdf_title = NULL)
#'
#' @export
#'
#' @importFrom e1071 skewness kurtosis
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#'

cyt_skku <- function(
  data,
  group_cols = NULL,
  pdf_title = NULL,
  print_res_raw = FALSE,
  print_res_log = FALSE
) {
  # Identify measurement columns.
  # If grouping columns are provided, exclude them from measurement columns.
  if (!is.null(group_cols)) {
    measure_cols <- setdiff(names(data), group_cols)
  } else {
    measure_cols <- names(data)
  }

  # Calculate a cutoff for log2 transformation:
  measure_mat <- data[, measure_cols, drop = FALSE]
  condt <- !is.na(measure_mat) & (measure_mat > 0)
  cutoff <- min(measure_mat[condt], na.rm = TRUE) / 10

  # Helper function: compute metrics for a given numeric vector Y
  # and grouping variable.
  compute_metrics <- function(Y, groups) {
    n <- tapply(Y, groups, function(x) sum(!is.na(x)))
    center <- tapply(Y, groups, mean, na.rm = TRUE)
    spread <- tapply(Y, groups, function(x) {
      sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
    })
    skew <- tapply(Y, groups, e1071::skewness, na.rm = TRUE)
    kurt <- tapply(Y, groups, e1071::kurtosis, na.rm = TRUE)

    data.frame(
      group = names(n),
      n = as.numeric(n),
      center = as.numeric(center),
      spread = as.numeric(spread),
      skewness = as.numeric(skew),
      kurtosis = as.numeric(kurt),
      stringsAsFactors = FALSE
    )
  }

  # Initialize lists to store results.
  raw_list <- list()
  log_list <- list()

  if (!is.null(group_cols)) {
    # If grouping columns are provided.
    # Determine treatment and group columns.
    if (length(group_cols) == 1) {
      treatment_col <- group_cols
      group_col <- group_cols
    } else {
      treatment_col <- group_cols[1]
      group_col <- group_cols[2]
    }

    # Get unique treatment levels.
    treatment_levels <- unique(data[[treatment_col]])

    # Loop over each treatment level.
    for (trt in treatment_levels) {
      subset_df <- data[data[[treatment_col]] == trt, ]
      groups <- subset_df[[group_col]]

      for (col in measure_cols) {
        # Raw data metrics.
        Y_raw <- subset_df[[col]]
        df_raw <- compute_metrics(Y_raw, groups)
        df_raw$measurement <- col
        df_raw$treatment <- trt
        raw_list[[paste(trt, col, sep = "_")]] <- df_raw

        # Log2-transformed data metrics.
        Y_log <- log2(Y_raw + cutoff)
        df_log <- compute_metrics(Y_log, groups)
        df_log$measurement <- col
        df_log$treatment <- trt
        log_list[[paste(trt, col, sep = "_")]] <- df_log
      }
    }

    # Combine results.
    raw_results <- do.call(rbind, raw_list)
    log_results <- do.call(rbind, log_list)
  } else {
    # No grouping columns provided: treat the entire dataset as one group.
    groups <- rep("overall", nrow(data))

    for (col in measure_cols) {
      # Raw data metrics.
      Y_raw <- data[[col]]
      df_raw <- compute_metrics(Y_raw, groups)
      df_raw$measurement <- col
      raw_list[[col]] <- df_raw

      # Log2-transformed data metrics.
      Y_log <- log2(Y_raw + cutoff)
      df_log <- compute_metrics(Y_log, groups)
      df_log$measurement <- col
      log_list[[col]] <- df_log
    }

    raw_results <- do.call(rbind, raw_list)
    log_results <- do.call(rbind, log_list)
  }

  # If a pdf title is provided, generate histograms using ggplot2.
  if (!is.null(pdf_title)) {
    pdf(file = pdf_title)
    on.exit(dev.off(), add = TRUE)
  }
  df_skew <- rbind(
    data.frame(value = raw_results$skewness, transformation = "Raw"),
    data.frame(value = log_results$skewness, transformation = "Log2")
  )
  df_kurt <- rbind(
    data.frame(value = raw_results$kurtosis, transformation = "Raw"),
    data.frame(value = log_results$kurtosis, transformation = "Log2")
  )

  p_skew <- ggplot2::ggplot(df_skew, aes(x = value, fill = transformation)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    ggplot2::labs(x = "Skewness", title = "Distribution of Skewness") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      legend.background = element_rect(fill = "white", colour = "white"),
      axis.title = element_text(color = "black", size = 12, face = "bold"),
      legend.title = element_text(color = "black", size = 10, face = "bold"),
      legend.text = element_text(color = "black")
    )

  p_kurt <- ggplot2::ggplot(df_kurt, aes(x = value, fill = transformation)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    ggplot2::labs(x = "Kurtosis", title = "Distribution of Kurtosis") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      legend.background = element_rect(fill = "white", colour = "white"),
      axis.title = element_text(color = "black", size = 12, face = "bold"),
      legend.title = element_text(color = "black", size = 10, face = "bold"),
      legend.text = element_text(color = "black")
    )

  gridExtra::grid.arrange(p_skew, p_kurt, ncol = 1)

  # Return results based on flags.
  if (print_res_raw && print_res_log) {
    return(list(raw = raw_results, log2 = log_results))
  } else if (print_res_raw) {
    return(raw_results)
  } else if (print_res_log) {
    return(log_results)
  } else {
    return(invisible(NULL))
  }
}
