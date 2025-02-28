#' Boxplot Function Enhanced for Specific Group Comparisons
#'
#' @param x.df A matrix or data frame of raw data.
#' @param Title A string representing the title of the PDF file.
#' @param mfRow A numeric vector of length two specifying the layout (rows and columns) for the plots on each page. Defaults to c(1,1).
#' @param scale Transformation option for continuous variables. Options are NULL (default) and "log2". When set to "log2", numeric columns are transformed using the log2 function.
#' @param yLim An optional numeric vector defining the y-axis limits for the plots.
#'
#' @description
#' This function generates a PDF file containing boxplots for each combination of numeric and factor variables in the provided data.
#' It first converts any character columns to factors and checks that the data contains at least one numeric and one factor column.
#' If the scale argument is set to "log2", all numeric columns are log2-transformed.
#' The function then creates boxplots using ggplot2 for each numeric variable grouped by each factor variable.
#'
#' @return A PDF file containing the boxplots.
#'
#' @examples
#' data.df <- cytodata[,-c(1,4)]
#' cyt.bp2(data.df, Title = "boxplot2.test2.pdf", scale = "log2")
#'
#' @export
#' @import ggplot2

cyt.bp2 <- function(x.df, Title, mfRow=c(1,1), scale = NULL, yLim=NULL) {

  # Convert any character variables to factors
  cat_vars <- sapply(x.df, is.character)
  if(any(cat_vars)){
    x.df[cat_vars] <- lapply(x.df[cat_vars], as.factor)
  }

  # Identify numeric and factor columns
  num_cols <- sapply(x.df, is.numeric)
  fac_cols <- sapply(x.df, is.factor)

  # Ensure there is at least one numeric and one factor column
  if(!any(num_cols)) stop("Data must contain at least one numeric column")
  if(!any(fac_cols)) stop("Data must contain at least one factor column")

  # Apply log2 transformation if log_scale is TRUE
  if(!is.null(scale) && scale == "log2"){
    x.df[num_cols] <- lapply(x.df[num_cols], function(x) log2(x))
  }

  # Generate boxplots
  pdf(file=Title)
  par(mfrow=mfRow, cex.axis=0.75)

  # Loop through each numeric column
  for(num_col in names(x.df)[num_cols]) {
    # Loop through each factor column
    for(fac_col in names(x.df)[fac_cols]) {
      # Create the boxplot
      p <- ggplot(x.df, aes_string(x = fac_col, y = num_col, fill = fac_col)) +
        geom_boxplot(alpha = 0.5) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        labs(title = paste0(num_col, " by ", fac_col),
             x = fac_col, y = paste0("Values of ", num_col)) +
        theme_minimal() +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
        guides(fill = guide_legend(title = NULL))

      if(!is.null(yLim)) {
        p <- p + ylim(yLim)
      }

      print(p)
    }
  }

  dev.off()
}
