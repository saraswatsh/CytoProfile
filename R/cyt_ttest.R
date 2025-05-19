#' Two Sample T-test Comparisons.
#'
#' This function performs pairwise comparisons between two groups for each combination
#' of a categorical predictor (with exactly two levels) and a continuous outcome variable.
#' It first converts any character variables in \code{data} to factors and, if specified,
#' applies a log2 transformation to the continuous variables. Depending on the value of
#' \code{scale}, the function conducts either a two-sample t-test (if \code{scale = "log2"})
#' or a Mann-Whitney U test (if \code{scale} is \code{NULL}). The resulting p-values are printed
#' and returned.
#'
#' @param data A matrix or data frame containing continuous and categorical variables.
#' @param scale A character specifying a transformation for continuous variables.
#'   Options are \code{NULL} (default) and \code{"log2"}. When \code{scale = "log2"},
#'   a log2 transformation is applied and a two-sample t-test is used; when \code{scale} is \code{NULL},
#'   a Mann-Whitney U test is performed.
#' @param verbose A logical indicating whether to print the p-values of the statistical tests.
#'   Default is \code{TRUE}.
#' @param format_output Logical. If TRUE, returns the results as a tidy data frame.
#'   Default is \code{FALSE}.
#'
#' @return If \code{format_output} is FALSE, returns a list of p-values (named by Outcome and Categorical variable).
#'   If TRUE, returns a data frame in a tidy format.
#'
#' @export
#'
#' @examples
#' data_df <- ExampleData1[, -c(3)]
#' data_df <- dplyr::filter(data_df, Group != "ND", Treatment != "Unstimulated")
#' # Test example
#' cyt_ttest(
#'  data_df[, c(1:2, 5:6)],
#'  scale = "log2",
#'  verbose = TRUE,
#'  format_output = TRUE
#' )
cyt_ttest <- function(
  data,
  scale = NULL,
  verbose = TRUE,
  format_output = FALSE
) {
  # Take input and store it as its own data frame
  x1_df <- as.data.frame(data)

  # Convert any character variables to factors
  cat_vars <- sapply(x1_df, is.character)
  if (any(cat_vars)) {
    x1_df[cat_vars] <- lapply(x1_df[cat_vars], as.factor)
  }

  # Identify categorical predictors and continuous variables
  cat_preds <- sapply(x1_df, is.factor)
  cont_vars <- sapply(x1_df, is.numeric)

  # Empty list to store test results
  test_results <- list()

  # Apply log2 transformation if scale is "log2"
  if (!is.null(scale) && scale == "log2") {
    x1_df[cont_vars] <- lapply(x1_df[cont_vars], function(x) log2(x))
  }

  # Loop over categorical predictors and continuous outcomes
  for (cat_var in names(x1_df)[cat_preds]) {
    for (outcome in names(x1_df)[cont_vars]) {
      # Check that the categorical variable has exactly two levels
      if (length(levels(x1_df[[cat_var]])) == 1) {
        stop("Categorical variable ", cat_var, " must have exactly two levels.")
      } else if (length(levels(x1_df[[cat_var]])) > 2) {
        stop("Categorical variable ", cat_var, " must have exactly two levels.")
      } else {
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
        } else if (sd(group1) == 0 || sd(group2) == 0) {
          warning(
            "Skipping test for ",
            outcome,
            " in ",
            cat_var,
            " due to low variance."
          )
          next
        }

        comparison_name <- paste(group_levels[1], "vs", group_levels[2])
        test_formula <- as.formula(paste(outcome, "~", cat_var))

        # normality check
        p1 <- tryCatch(
          stats::shapiro.test(group1)$p.value,
          error = function(e) 0
        )
        p2 <- tryCatch(
          stats::shapiro.test(group2)$p.value,
          error = function(e) 0
        )

        # pick the test
        if (p1 > 0.05 && p2 > 0.05) {
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
        # Store the p-value in the results list
        key <- paste(outcome, cat_var, sep = "_")
        test_results[[key]] <- tt
      }
    }
  }
  if (length(test_results) == 0) {
    return("No valid tests were performed.")
  }
  # Return results in tidy format if requested, otherwise as a list
  if (!format_output && verbose) {
    return(test_results)
  }

  if (format_output && verbose) {
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
        stat <- if (!is.null(tt$statistic)) unname(tt$statistic)[1] else
          NA_real_

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
    res <- list(results = out_df)
    return(res)
  }
}
