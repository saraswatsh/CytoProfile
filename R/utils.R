# Utility functions used throughout the CytoProfile package.

#' Apply a scale transformation to numeric columns
#'
#' This helper function applies a chosen scaling or transformation to
#' specified numeric columns in a data frame.  Supported built‑in
#' transformations include no transformation ("none"), log2, log10,
#' and z‑score scaling.  A custom function can also be supplied to
#' perform arbitrary transformations.
#'
#' @param data A data.frame or matrix containing the data to be
#'   transformed.
#' @param columns A character vector of column names to transform.  If
#'   NULL (default) all numeric columns will be transformed.
#' @param scale A character string specifying the transformation to
#'   apply.  Possible values are "none", "log2", "log10",
#'   "zscore", or "custom".  When set to "custom" the function
#'   specified in `custom_fn` will be applied to the columns.
#' @param custom_fn A function that takes a numeric vector and returns a
#'   transformed numeric vector.  Only used when `scale = "custom"`.
#' @return A data.frame with the same dimensions as `data` with
#'   transformed numeric columns.
#' @importFrom stats p.adjust
#' @export
apply_scale <- function(
  data,
  columns = NULL,
  scale = c("none", "log2", "log10", "zscore", "custom"),
  custom_fn = NULL
) {
  data <- as.data.frame(data)
  scale <- match.arg(scale)

  # Identify numeric columns if none supplied
  if (is.null(columns)) {
    num_cols <- names(data)[sapply(data, is.numeric)]
  } else {
    num_cols <- columns
  }

  # Bail early if no numeric columns
  if (length(num_cols) == 0) {
    return(data)
  }

  # Transformation functions
  trans_fun <- switch(
    scale,
    none = function(x) x,
    log2 = function(x) log2(x),
    log10 = function(x) log10(x),
    zscore = function(x) {
      if (sd(x, na.rm = TRUE) == 0) {
        return(x - mean(x, na.rm = TRUE))
      }
      (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    },
    custom = {
      if (is.null(custom_fn) || !is.function(custom_fn)) {
        stop(
          "When scale = 'custom', a valid function must be provided to 'custom_fn'."
        )
      }
      custom_fn
    }
  )

  for (col in num_cols) {
    data[[col]] <- trans_fun(data[[col]])
  }
  data
}

#' Adjust p-values using a specified method
#'
#' A thin wrapper around `stats::p.adjust` that defaults to the
#' Benjamini–Hochberg procedure.  Useful for unifying multiple
#' testing adjustments across the package.
#'
#' @param p_values A numeric vector of raw p-values.
#' @param method A character string specifying the p‑value adjustment
#'   method.  Passed directly to `p.adjust`.  Defaults to "BH"
#'   (Benjamini–Hochberg).  See `p.adjust.methods` for other
#'   options.
#' @return A numeric vector of adjusted p-values of the same length as
#'   `p_values`.
#' @export
adjust_p <- function(p_values, method = "BH") {
  p.adjust(p_values, method = method)
}
