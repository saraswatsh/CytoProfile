#' ANOVA analysis on all continuous variables within the data.
#' @param data A data frame or matrix containing both categorical and continuous variables.
#'   Character columns are converted to factors; all factor columns are used as predictors,
#'   while numeric columns are used as continuous outcomes.
#' @description
#' This function performs an ANOVA for each continuous variable against every categorical predictor
#' in the input data. For each valid predictor (i.e., with more than one level and no more than 10 levels),
#' it conducts Tukey's HSD test and extracts the adjusted p-values for pairwise comparisons.
#' @return A list of adjusted p-values from Tukey's HSD tests for each combination of continuous outcome
#'   and categorical predictor. The list elements are named in the format "outcome_predictor".
#' @examples
#' \dontrun{
#' data("cytodata")
#' cyt.anova(cytodata[, c(2:4, 5:6)])
#' }
#' @export
cyt.anova <- function(data) {
  # Take input and store it as its own data frame
  x1.df <- data
  # Convert any char variables to factors
  cat_vars <- sapply(x1.df, is.character)
  if (any(cat_vars)) {
    x1.df[cat_vars] <- lapply(x1.df[cat_vars], as.factor)
  }
  # Categorical Predictors
  cat_preds <- sapply(x1.df, is.factor)
  # Create a list to store column names with numeric data
  cont_vars <- sapply(x1.df, is.numeric)
  # Empty list to store ANOVA and Tukey results
  tukey_results <- list()

  # ANOVA and Tukey Comparisons
  for (cat_var in names(x1.df)[cat_preds]) {
    for (outcome in names(x1.df)[cont_vars]) {
      if (length(levels(x1.df[[cat_var]])) == 1) {
        next
      } else if (length(levels(x1.df[[cat_var]])) > 10) {
        next
      } else {
        # Perform ANOVA tests on each continuous variable with the current categorical variable
        model <- aov(as.formula(paste(outcome, "~", cat_var)), data = x1.df)
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
  return(tukey_results)
}
