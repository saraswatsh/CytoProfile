#' Pairwise Univariate Tests Between Two Groups
#'
#' @description
#' `cyt_univariate` supports additional scaling options and explicit
#' choice of statistical test.
#' For each categorical predictor with exactly two levels and each
#' numeric outcome, a two‑sample t‑test or Wilcoxon rank–sum test is
#' performed.  Results are returned either as a list of test objects
#' or, if `format_output = TRUE`, as a tidy data frame with one
#' row per comparison.
#'
#' @param data A data frame or matrix containing both categorical
#'   and numeric variables.
#' @param scale A character specifying a transformation to apply to
#'   numeric variables prior to testing.  Choices are `NULL` (no
#'   transformation), "log2", "log10", "zscore", or
#'   "custom".  When set to "custom", supply a function via
#'   `custom_fn`.
#' @param method Character specifying the test to perform.  Use
#'   "auto" (default) to select between t‑test and Wilcoxon based
#'   on Shapiro–Wilk normality tests for each outcome; "ttest" to
#'   always use Student’s t‑test; or "wilcox" to always use the
#'   Wilcoxon rank–sum test.
#' @param verbose Logical indicating whether to return the results.
#'   Provided for backward compatibility but has no effect on printing.
#' @param format_output Logical.  If `TRUE`, returns the results as
#'   a tidy data frame; if `FALSE` (default), returns a list of
#'   test objects similar to the original function.
#' @param custom_fn A function to apply when `scale = "custom"`.
#' @return If `format_output = FALSE`, a named list of test objects
#'   keyed by "Outcome_Categorical".  If `format_output = TRUE`, a
#'   data frame with columns `Outcome`, `Categorical`, `Comparison`,
#'   `Test`, `Estimate`, `Statistic`, and `P_value`.
#' @examples
#' data_df <- ExampleData1[, -c(3)]
#' data_df <- dplyr::filter(data_df, Group != "ND", Treatment != "Unstimulated")
#' cyt_univariate(data_df[, c(1:2, 5:6)], scale = "log2",
#'                method = "auto", format_output = TRUE)
#' @author Shubh Saraswat
#' @import stats
#' @export
cyt_univariate <- function(
  data,
  scale = NULL,
  method = c("auto", "ttest", "wilcox"),
  verbose = TRUE,
  format_output = FALSE,
  custom_fn = NULL
) {
  method <- match.arg(method)
  # Convert to data frame
  x1_df <- as.data.frame(data)
  # Convert character variables to factors
  cat_vars <- sapply(x1_df, is.character)
  if (any(cat_vars)) {
    x1_df[cat_vars] <- lapply(x1_df[cat_vars], as.factor)
  }
  # Identify categorical predictors and continuous variables
  cat_preds <- sapply(x1_df, is.factor)
  cont_vars <- sapply(x1_df, is.numeric)
  # Apply scaling to numeric variables
  if (!is.null(scale)) {
    if (scale == "log2") {
      x1_df[cont_vars] <- lapply(x1_df[cont_vars], function(x) log2(x))
    } else if (scale == "log10") {
      x1_df[cont_vars] <- lapply(x1_df[cont_vars], function(x) log10(x))
    } else if (scale == "zscore") {
      x1_df[cont_vars] <- lapply(x1_df[cont_vars], function(x) {
        mu <- mean(x, na.rm = TRUE)
        sdv <- stats::sd(x, na.rm = TRUE)
        if (sdv == 0) {
          return(rep(0, length(x)))
        }
        (x - mu) / sdv
      })
    } else if (scale == "custom") {
      if (is.null(custom_fn) || !is.function(custom_fn)) {
        stop("When scale = 'custom', a valid custom_fn must be provided.")
      }
      x1_df[cont_vars] <- lapply(x1_df[cont_vars], custom_fn)
    }
  }
  # Empty list to store test results
  test_results <- list()
  # Loop over categorical predictors and continuous outcomes
  for (cat_var in names(x1_df)[cat_preds]) {
    # Only consider predictors with exactly two levels
    if (length(levels(x1_df[[cat_var]])) != 2) {
      next
    }
    for (outcome in names(x1_df)[cont_vars]) {
      group_levels <- levels(x1_df[[cat_var]])
      group1 <- x1_df[[outcome]][x1_df[[cat_var]] == group_levels[1]]
      group2 <- x1_df[[outcome]][x1_df[[cat_var]] == group_levels[2]]
      # Check for sufficient data and variance in both groups
      if (length(group1) < 2 || length(group2) < 2) {
        warning(
          "Skipping test for ",
          outcome,
          " in ",
          cat_var,
          " due to insufficient data in one of the groups."
        )
        next
      } else if (stats::sd(group1) == 0 || stats::sd(group2) == 0) {
        warning(
          "Skipping test for ",
          outcome,
          " in ",
          cat_var,
          " due to low variance."
        )
        next
      }
      # Determine test to use
      test_used <- switch(
        method,
        auto = {
          p1 <- tryCatch(
            stats::shapiro.test(group1)$p.value,
            error = function(e) 0
          )
          p2 <- tryCatch(
            stats::shapiro.test(group2)$p.value,
            error = function(e) 0
          )
          if (p1 > 0.05 && p2 > 0.05) "ttest" else "wilcox"
        },
        ttest = "ttest",
        wilcox = "wilcox"
      )
      # Perform the test
      if (test_used == "ttest") {
        tt <- stats::t.test(
          stats::as.formula(paste(outcome, "~", cat_var)),
          data = x1_df
        )
      } else {
        tt <- stats::wilcox.test(
          stats::as.formula(paste(outcome, "~", cat_var)),
          data = x1_df,
          conf.int = TRUE,
          exact = FALSE
        )
      }
      key <- paste(outcome, cat_var, sep = "_")
      test_results[[key]] <- tt
    }
  }
  if (length(test_results) == 0) {
    return("No valid tests were performed.")
  }
  # Return results in tidy format if requested
  if (!format_output) {
    return(test_results)
  }
  # Format into data frame
  out_df <- do.call(
    rbind,
    lapply(names(test_results), function(key) {
      tt <- test_results[[key]]
      parts <- strsplit(key, "_")[[1]]
      outcome <- parts[1]
      cat_var <- parts[2]
      lvls <- levels(x1_df[[cat_var]])
      comp <- paste(lvls[1], "vs", lvls[2])
      est <- if (!is.null(tt$estimate)) unname(tt$estimate)[1] else NA_real_
      stat <- if (!is.null(tt$statistic)) unname(tt$statistic)[1] else NA_real_
      data.frame(
        Outcome = outcome,
        Categorical = cat_var,
        Comparison = comp,
        Test = tt$method,
        Estimate = round(est, 3),
        Statistic = round(stat, 3),
        P_value = round(tt$p.value, 3),
        stringsAsFactors = FALSE
      )
    })
  )
  return(out_df)
}
