#' ANOVA Analysis on Continuous Variables. `r lifecycle::badge("deprecated")`
#'
#' This function performs an analysis of variance (ANOVA) for each continuous
#' variable against every categorical predictor in the input data. Character
#' columns are automatically converted to factors; all factor columns are used
#' as predictors while numeric columns are used as continuous outcomes.
#' For each valid predictor (i.e., with more than one level and no more than 10 levels),
#' Tukey's Honest Significant Difference (HSD) test is conducted and the adjusted
#' p-values for pairwise comparisons are extracted.
#'
#' @param data A data frame or matrix containing both categorical and continuous variables.
#'   Character columns will be converted to factors and used as predictors, while numeric columns
#'   will be used as continuous outcomes.
#' @param format_output Logical. If TRUE, returns the results as a tidy data frame instead of a list.
#'   Default is FALSE.
#'
#' @return If \code{format_output} is FALSE (default), a list of adjusted p-values from Tukey's HSD tests
#'   for each combination of continuous outcome and categorical predictor. List elements are named
#'   in the format "Outcome_Categorical".
#'   If \code{format_output} is TRUE, a data frame in a tidy format.
#' @author Shubh Saraswat
#' @export
#'
#' @examples
#' data("ExampleData1")
#' cyt_anova(ExampleData1[, c(1:2, 5:6)], format_output = TRUE)
#'

cyt_anova <- function(data, format_output = FALSE) {
  lifecycle::deprecate_warn(
    "0.2.4", # version when deprecation begins
    "CytoProfile::cyt_anova()",
    "CytoProfile::cyt_univariate_multi()"
  )
  names(data) <- make.names(names(data), unique = TRUE)
  # Convert input data to a data frame
  x1_df <- as.data.frame(data)

  # Convert character variables to factors
  cat_vars <- sapply(x1_df, is.character)
  if (any(cat_vars)) {
    x1_df[cat_vars] <- lapply(x1_df[cat_vars], as.factor)
  }

  # Identify categorical predictors and continuous outcomes
  cat_preds <- sapply(x1_df, is.factor)
  cont_vars <- sapply(x1_df, is.numeric)

  # List to store Tukey results
  tukey_results <- list()

  # Loop over each categorical predictor and continuous outcome
  for (cat_var in names(x1_df)[cat_preds]) {
    # Skip factors with only one level or more than 10 levels
    num_levels <- length(levels(x1_df[[cat_var]]))
    if (num_levels <= 1 || num_levels > 10) {
      next
    }

    for (outcome in names(x1_df)[cont_vars]) {
      # Build the ANOVA model and perform Tukey's HSD test
      model <- aov(as.formula(paste(outcome, "~", cat_var)), data = x1_df)
      tukey_result <- TukeyHSD(model)
      # Extract the adjusted p-values for the current predictor
      p_vals <- tukey_result[[cat_var]][, "p adj"]

      # Store results using a key in the format Outcome_Categorical
      result_key <- paste(outcome, cat_var, sep = "_")
      tukey_results[[result_key]] <- round(p_vals, 4)
    }
  }

  # Check if any comparisons were performed
  if (length(tukey_results) == 0) {
    return(
      "No valid comparisons were performed. Check that your data has numeric columns and factors with 2-10 levels."
    )
  }

  # Return tidy data frame if requested
  if (!format_output) {
    tukey_results
  } else {
    out_df <- data.frame(
      Outcome = character(),
      Categorical = character(),
      Comparison = character(),
      P_adj = numeric(),
      stringsAsFactors = FALSE
    )

    # Loop through each result in the list and add rows to the data frame
    for (key in names(tukey_results)) {
      parts <- unlist(strsplit(key, "_"))
      outcome <- parts[1]
      cat_var <- parts[2]
      p_vals <- tukey_results[[key]]
      for (comp in names(p_vals)) {
        out_df <- rbind(
          out_df,
          data.frame(
            Outcome = outcome,
            Categorical = cat_var,
            Comparison = comp,
            P_adj = p_vals[comp],
            stringsAsFactors = FALSE
          )
        )
      }
    }
    out_df
  }
}
