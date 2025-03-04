#' Distribution of the Data Set Shown by Skewness and Kurtosis
#'
#' @description This function computes summary statistics—including sample size, mean, standard error, skewness, and kurtosis—for each numeric measurement column in a data set. If grouping columns are provided via `group.cols`, the function computes the metrics separately for each group defined by the combination of these columns (using the first element as the treatment variable and the second as the grouping variable, or the same column for both if only one is given). When no grouping columns are provided, the entire data set is treated as a single group ("Overall"). A log2 transformation (using a cutoff equal to one-tenth of the smallest positive value in the data) is applied to generate alternative metrics. Histograms showing the distribution of skewness and kurtosis for both raw and log2-transformed data are then generated and saved to a PDF if a file name is provided.
#'
#' @param x.df A matrix or data frame containing the raw data. If `group.cols` is specified, the columns with names in `group.cols` are treated as grouping variables and all other columns are assumed to be numeric measurement variables.
#' @param group.cols A character vector specifying the names of the grouping columns. When provided, the first element is treated as the treatment variable and the second as the group variable. If not provided, the entire data set is treated as one group.
#' @param Title A character string specifying the file name for the PDF file in which the histograms will be saved. If omitted, the plots are generated on the current graphics device.
#' @param printResRaw Logical. If \code{TRUE}, the function returns and prints the computed summary statistics for the raw data. Default is \code{FALSE}.
#' @param printResLog Logical. If \code{TRUE}, the function returns and prints the computed summary statistics for the log2-transformed data. Default is \code{FALSE}.
#'
#' @return The function generates histograms of skewness and kurtosis for both raw and log2-transformed data. Additionally, if either \code{printResRaw} and/or \code{printResLog} is \code{TRUE}, the function returns the corresponding summary statistics as a data frame or a list of data frames.
#'
#' @details A cutoff is computed as one-tenth of the minimum positive value among all numeric measurement columns to avoid taking logarithms of zero. When grouping columns are provided, the function loops over unique treatment levels (using the first element of `group.cols` as treatment and the second as group, if available) and computes the metrics for each measurement column within each subgroup. Without grouping columns, the entire data set is analyzed as one overall group.
#'
#' @examples
#' \dontrun{
#' # Example with grouping columns (e.g., "Group" and "Treatment")
#' data(cytodata)
#' cyt.skku(cytodata[,-c(1,3,4)], Title = "Skew_and_Kurtosis.pdf", group.cols = c("Group"))
#'
#' # Example without grouping columns (analyzes the entire data set)
#' cyt.skku(cytodata[,-c(1,4)], Title = "Skew_and_Kurtosis_Overall.pdf")
#' }
#'
#' @export
#'
#' @import e1071
#' @import ggplot2
#' @import gridExtra

cyt.skku <- function(data.df, group.cols = NULL, Title = NULL,
                     printResRaw = FALSE, printResLog = FALSE) {

  # Identify measurement columns.
  # If grouping columns are provided, exclude them from measurement columns.
  if (!is.null(group.cols)) {
    measure.cols <- setdiff(names(data.df), group.cols)
  } else {
    measure.cols <- names(data.df)
  }

  # Calculate a cutoff for log2 transformation:
  measure.mat <- data.df[, measure.cols, drop = FALSE]
  condt <- !is.na(measure.mat) & (measure.mat > 0)
  cutoff <- min(measure.mat[condt], na.rm = TRUE) / 10

  # Helper function: compute metrics for a given numeric vector Y and grouping variable.
  compute_metrics <- function(Y, groups) {
    n      <- tapply(Y, groups, function(x) sum(!is.na(x)))
    center <- tapply(Y, groups, mean, na.rm = TRUE)
    spread <- tapply(Y, groups, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    skew   <- tapply(Y, groups, skewness, na.rm = TRUE)
    kurt   <- tapply(Y, groups, kurtosis, na.rm = TRUE)

    data.frame(
      Group = names(n),
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

  if (!is.null(group.cols)) {
    # If grouping columns are provided.
    # Determine treatment and group columns.
    if (length(group.cols) == 1) {
      treatment.col <- group.cols
      group.col <- group.cols
    } else {
      treatment.col <- group.cols[1]
      group.col <- group.cols[2]
    }

    # Get unique treatment levels.
    treatment_levels <- unique(data.df[[treatment.col]])

    # Loop over each treatment level.
    for (trt in treatment_levels) {
      subset_df <- data.df[data.df[[treatment.col]] == trt, ]
      groups <- subset_df[[group.col]]

      for (col in measure.cols) {
        # Raw data metrics.
        Y_raw <- subset_df[[col]]
        df_raw <- compute_metrics(Y_raw, groups)
        df_raw$Measurement <- col
        df_raw$Treatment <- trt
        raw_list[[paste(trt, col, sep = "_")]] <- df_raw

        # Log2-transformed data metrics.
        Y_log <- log2(Y_raw + cutoff)
        df_log <- compute_metrics(Y_log, groups)
        df_log$Measurement <- col
        df_log$Treatment <- trt
        log_list[[paste(trt, col, sep = "_")]] <- df_log
      }
    }

    # Combine results.
    raw_results <- do.call(rbind, raw_list)
    log_results <- do.call(rbind, log_list)

  } else {
    # No grouping columns provided: treat the entire dataset as one group.
    groups <- rep("Overall", nrow(data.df))

    for (col in measure.cols) {
      # Raw data metrics.
      Y_raw <- data.df[[col]]
      df_raw <- compute_metrics(Y_raw, groups)
      df_raw$Measurement <- col
      raw_list[[col]] <- df_raw

      # Log2-transformed data metrics.
      Y_log <- log2(Y_raw + cutoff)
      df_log <- compute_metrics(Y_log, groups)
      df_log$Measurement <- col
      log_list[[col]] <- df_log
    }

    raw_results <- do.call(rbind, raw_list)
    log_results <- do.call(rbind, log_list)
  }

  # If a Title is provided, generate histograms using ggplot2.
  if (!is.null(Title)) {
    pdf(file = Title)
    library(ggplot2)
    library(gridExtra)

    df_skew <- rbind(
      data.frame(Value = raw_results$skewness, Transformation = "Raw"),
      data.frame(Value = log_results$skewness, Transformation = "Log2")
    )
    df_kurt <- rbind(
      data.frame(Value = raw_results$kurtosis, Transformation = "Raw"),
      data.frame(Value = log_results$kurtosis, Transformation = "Log2")
    )

    p_skew <- ggplot(df_skew, aes(x = Value, fill = Transformation)) +
      geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
      labs(x = "Skewness", title = "Distribution of Skewness") +
      theme_minimal()

    p_kurt <- ggplot(df_kurt, aes(x = Value, fill = Transformation)) +
      geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
      labs(x = "Kurtosis", title = "Distribution of Kurtosis") +
      theme_minimal()

    grid.arrange(p_skew, p_kurt, ncol = 1)
    dev.off()
  }

  # Return results based on flags.
  if (printResRaw & printResLog) {
    return(list(raw = raw_results, log2 = log_results))
  } else if (printResRaw) {
    return(raw_results)
  } else if (printResLog) {
    return(log_results)
  } else {
    return(invisible(NULL))
  }
}
