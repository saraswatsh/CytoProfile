#' Distribution of the groups within the continuous variable.
#'
#' @param x.df A matrix or data frame consisting of continuous and categorical variables.
#' @param cont_var Continuous variable.
#' @param group Column name to be color coded by.
#'
#' @return Prints the histograms.
#' @export
plot_group_distribution = function(x.df, continuous_var, group_var) {
  p = ggplot(x.df, aes(x = !!sym(continuous_var), fill = !!sym(group_var))) +
    geom_histogram(binwidth = 1, position = "dodge") +
    labs(title = paste("Distribution of", continuous_var, "by", group_var),
         x = continuous_var, y = "Frequency") +
    theme_minimal()

  print(p)
}
