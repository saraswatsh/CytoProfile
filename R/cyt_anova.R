#' ANOVA Analysis on Continuous Variables.
#'
#' This function performs an analysis of variance (ANOVA) for each continuous
#' variable against every categorical predictor in the input data. Character
#' columns are automatically converted to factors; all factor columns are used
#' as predictors while numeric columns are used as continuous outcomes.
#' For each valid predictor (i.e., with more than one level and no more
#' than 10 levels), Tukey's Honest Significant Difference (HSD) test is
#' conducted and the adjusted p-values for pairwise comparisons are extracted.
#'
#' @param data A data frame or matrix containing both categorical and
#' continuous variables. Character columns will be converted to factors and
#' used as predictors, while numeric columns will be used as continuous outcomes.
#'
#' @return A list of adjusted p-values from Tukey's HSD tests for each
#' combination of continuous outcome and categorical predictor. List
#' elements are named in the format "outcome_predictor".
#'
#' @examples
#' # Loading data
#' data("cytodata")
#' # Perform ANOVA on selected columns of the cytodata dataset
#' anova_results <- cyt_anova(cytodata[, c(2:4, 5:6)])
#' print(anova_results)
#'
#' @export
cyt_anova <- function(data) {
  # Take input and store it as its own data frame
  x1_df <- data

  # Convert any char variables to factors
  cat_vars <- sapply(x1_df, is.character)
  if (any(cat_vars)) {
    x1_df[cat_vars] <- lapply(x1_df[cat_vars], as.factor)
  }

  # Categorical Predictors
  cat_preds <- sapply(x1_df, is.factor)

  # Create a list to store column names with numeric data
  cont_vars <- sapply(x1_df, is.numeric)

  # Empty list to store ANOVA and Tukey results
  tukey_results <- list()

  # ANOVA and Tukey Comparisons
  for (cat_var in names(x1_df)[cat_preds]) {
    for (outcome in names(x1_df)[cont_vars]) {
      if (length(levels(x1_df[[cat_var]])) == 1) {
        next
      } else if (length(levels(x1_df[[cat_var]])) > 10) {
        next
      } else {
        # Perform ANOVA tests on each continuous variable with the
        # current categorical variable
        model <- aov(as.formula(paste(outcome, "~", cat_var)), data = x1_df)

        # Tukey summary
        tukey_result <- TukeyHSD(model)

        # Extract p-values and store them in the list
        p_values_cat_var <- tukey_result[[cat_var]][, "p adj"]

        # Store the p-values in the results list
        result_key <- paste(outcome, cat_var, sep = "_")
        tukey_results[[result_key]] <- p_values_cat_var
      }
    }
  }

  # Implicit return (if you prefer not to remove, keep return(tukey_results))
  tukey_results
}
