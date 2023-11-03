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
#' cyt_anova(data.df)
#'
#' @export

cyt_anova = function(x.df) {

  # Setting the rest of the variables as outcome
  cont_vars = names(x.df)[3:ncol(x.df)]

  # Empty list to store ANOVA and Tukey results
  tukey_results = list()

  # Getting variables from the first two columns
  if(length(levels(as.factor(x.df[,1]))) > 1 & length(levels(as.factor(x.df[,2]))) > 1 ){
    pred1 = names(x.df)[1]
    pred2 = names(x.df)[2]
    for(outcome in cont_vars){
      # Perform ANOVA tests on each continuous variable
      model = aov(as.formula(paste(outcome, "~", pred1, "+", pred2)), data = x.df)
      # Tukey summary
      tukey_result = TukeyHSD(model)
      # Extract p-values and store them in the list
      p_values_pred1 = tukey_result[[pred1]][, "p adj"]
      p_values_pred2 = tukey_result[[pred2]][, "p adj"]
      # Store the p-values in the results list
      tukey_results[[outcome]] = list(p_values_pred1, p_values_pred2)
    }
  }
  if(length(levels(as.factor(x.df[,1]))) == 1 & length(levels(as.factor(x.df[,2]))) > 1){
    pred1 = names(x.df)[2]
    for(outcome in cont_vars){
      # Perform ANOVA tests on each continuous variable
      model = aov(as.formula(paste(outcome, "~", pred1)), data = x.df)
      # Tukey summary
      tukey_result = TukeyHSD(model)
      # Extract p-values and store them in the list
      p_values_pred1 = tukey_result[[pred1]][, "p adj"]
      # Store the p-values in the results list
      tukey_results[[outcome]] = list(p_values_pred1)
    }
  }
  if(length(levels(as.factor(x.df[,1]))) > 1 & length(levels(as.factor(x.df[,2]))) == 1){
    pred1 = names(x.df)[1]
    for(outcome in cont_vars){
      # Perform ANOVA tests on each continuous variable
      model = aov(as.formula(paste(outcome, "~", pred1)), data = x.df)
      # Tukey summary
      tukey_result = TukeyHSD(model)
      # Extract p-values and store them in the list
      p_values_pred1 = tukey_result[[pred1]][, "p adj"]
      # Store the p-values in the results list
      tukey_results[[outcome]] = list(p_values_pred1)
    }
  }
  return(tukey_results)
}
