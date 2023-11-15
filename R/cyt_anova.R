#########################################################################################################
# Function to generate ANOVA analysis and results
# Author: Shubh Saraswat
# Arguments:
#   x.df: a matrix or data frame with groups, stimulation, and continuous variables
#########################################################################################################

#' ANOVA analysis on all continuous variables within the data.
#'
#' @param x.df A matrix or data frame consisting of continuous and categorical variables.
#' @description
#' This function produces and prints list of p-values obtained from Tukey comparisons. It assumes that
#' the first two columns are categorical variables whether it is group name or treatment name and uses
#' those as predictors. The rest of the columns in the data set are assumed to be continuous variables to
#' be used as the outcomes.
#' @return Prints the p-values of comparisons conducted using the Tukey test of the ANOVA model.
#' @examples
#' data(iris)
#' cyt_anova(iris)
#' @export
cyt.anova = function(x.df) {
  # Take input and store it as it's own data frame
  x1.df = x.df
  # Convert any char variables to factors
  cat_vars = sapply(x1.df, is.character)
  if(any(cat_vars)){
    x1.df[cat_vars] = lapply(x1.df[cat_vars], as.factor)
  }
  # Categorical Predictors
  cat_preds = sapply(x1.df, is.factor)
  # Create a list to store column names with numeric data
  cont_vars <- sapply(x1.df, is.numeric)
  # Empty list to store ANOVA and Tukey results
  tukey_results = list()

  # ANOVA and TUkey Comparisons
  for(cat_var in names(x1.df)[cat_preds]){
    for(outcome in names(x1.df)[cont_vars]){
      if(length(levels(cat_var)) == 1){
        next
      }
      else if(length(levels(cat_var)) > 10){
        next
      }
      else{
        # Perform ANOVA tests on each continuous variable with the current categorical variable
        model <- aov(as.formula(paste(outcome, "~", cat_var)), data = x1.df)
        # Tukey summary
        tukey_result = TukeyHSD(model)

        # Extract p-values and store them in the list
        p_values_cat_var <- tukey_result[[cat_var]][, "p adj"]

        # Store the p-values in the results list
        result_key <- paste(outcome, cat_var, sep = "_")
        tukey_results[[result_key]] = p_values_cat_var
      }
    }
  }
  return(tukey_results)
}
