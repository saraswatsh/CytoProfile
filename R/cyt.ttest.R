#' Two Sample T-test Comparisons
#'
#' @param x.df A matrix or data frame consisting of continuous and categorical variables.
#' @param scale Transformation option for continuous variables. Options are NULL (default) and "log2".
#' When scale = "log2". the function will use two sample t-test, when scale = NULL, the function will use Mann-Whitney U test.
#'
#' @return Prints the p-values of comparisons conducted.
#' @export
#'
#' @examples
#' data.df = cytodata[,-c(1,4)]
#' data.df = filter(data.df, Group != "ND", Treatment != "Unstimulated")
#' cyt.ttests(data.df, scale = "log2")
#' cyt.ttests(data.df)
#'
cyt.ttests = function(x.df, scale = NULL) {
  # Take input and store it as its own data frame
  x1.df = x.df
  # Convert any character variables to factors
  cat_vars = sapply(x1.df, is.character)
  if(any(cat_vars)){
    x1.df[cat_vars] = lapply(x1.df[cat_vars], as.factor)
  }
  # Categorical Predictors
  cat_preds = sapply(x1.df, is.factor)
  # Create a list to store column names with numeric data
  cont_vars <- sapply(x1.df, is.numeric)
  # Empty list to store test results
  test_results = list()
  # Apply log2 transformation if scale is "log2"
  if(!is.null(scale) && scale == "log2"){
    x1.df[cont_vars] = lapply(x1.df[cont_vars], function(x) log2(x))
  }
  # Perform tests based on user input
  for(cat_var in names(x1.df)[cat_preds]){
    for(outcome in names(x1.df)[cont_vars]){
      if(length(levels(x1.df[[cat_var]])) == 1){
        stop("Must have two levels in the categorical variable")
      }
      else if(length(levels(x1.df[[cat_var]])) > 2){
        stop("Must have two levels in the categorical variable")
      }
      else{

        if(length(levels(x1.df[[cat_var]])) == 2){
          # Two-sample t-test or Mann-Whitney U test
          group1 <- x1.df[[outcome]][x1.df[[cat_var]] == levels(x1.df[[cat_var]])[1]]
          group2 <- x1.df[[outcome]][x1.df[[cat_var]] == levels(x1.df[[cat_var]])[2]]

          if(length(group1) < 2 || length(group2) < 2){
            cat("Skipping test due to insufficient data in one of the groups\n")
            next
          }
          # Check for low variance
          else if(sd(group1) == 0 || sd(group2) == 0){
            cat("Skipping test due to low variance in", outcome, "\n")
            next
          }
          comparison_name <- paste(levels(x1.df[[cat_var]])[1], "vs", levels(x1.df[[cat_var]])[2])
          formula <- as.formula(paste(outcome, "~", cat_var))

          if(!is.null(scale) && scale == "log2"){
            # Perform two-sample t-test
            test_result <- t.test(formula, data = x1.df)
            cat(paste0("T-test p-value for ", comparison_name, " on ", outcome, ": ", test_result$p.value, "\n"))
          }
          else if(is.null(scale)) {
            # Perform Mann-Whitney U test
            test_result <- wilcox.test(formula, data = x1.df)
            cat(paste0("Mann-Whitney U test p-value for ", comparison_name, " on ", outcome, ": ", test_result$p.value, "\n"))
          }

          # Extract p-value and store in results list
          result_key <- paste(outcome, cat_var, sep = "_")
          test_results[[result_key]] <- test_result$p.value
        }
      }
    }
  }
}
