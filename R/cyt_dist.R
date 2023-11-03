###############################################################################
############### Skewness and Kurtosis #########################################
###############################################################################

#' Distribution of the data set
#'
#' @param x.df A matrix or data frame.
#' @description
#' The function takes in a data frame and subsets the numeric columns from the data which is then
#' used to calculate the skewed and kurtosis values. The values are then plotted using histograms
#' to visualize the distribution.
#'
#' @return Prints histograms of Skewness and Kurtosis of the continuous variables using raw data and log2 transformation.
#' @export
#'
cyt_dist <- function(x.df) {
  # Create a list to store column names with numeric data
  numeric_columns <- sapply(x.df, is.numeric)

  # Create vectors to store skewness and kurtosis values
  skewness_values <- numeric(length(numeric_columns))
  kurtosis_values <- numeric(length(numeric_columns))
  skewness_log2 = numeric(length(numeric_columns))
  kurtosis_log2 = numeric(length(numeric_columns))
  # Iterate through each numeric column and calculate skewness and kurtosis
  for (i in 1:sum(numeric_columns)) {
    col <- names(x.df)[numeric_columns][i]
    skewness_values[i] <- skewness(x.df[[col]])
    kurtosis_values[i] <- kurtosis(x.df[[col]])
    skewness_log2[i] = skewness(log2(x.df[[col]]))
    kurtosis_log2[i] = kurtosis(log2(x.df[[col]]))
  }

  # Plot histograms of skewness and kurtosis
  par(mfrow = c(1, 2))
  hist(skewness_values, main = "Skewness Histogram", xlab = "Skewness", col = "gray")
  hist(skewness_log2, main = "Skewness log2 Histogram", xlab = "Skewness", col = "gray")
  hist(kurtosis_values, main = "Kurtosis Histogram", xlab = "Kurtosis", col = "gray")
  hist(kurtosis_log2, main = "Kurtosis log2 Histogram", xlab = "Kurtosis", col = "gray")
}
