#############################################################################
## Function to generate a PDF file to show the measured value by            #
## column of data                                                           #
## Author: Xiaohua Douglas Zhang, January 2022                              #
## Arguments:                                                               #
##   data: a matrix of measured continuous value with row and column names  #
##   Title: names for the PDF file                                          #
##   bin.size: the number of boxplots shown in a single pane                #
##   yLim: a range for the y-axis to be plotted                             #
## Output:                                                                  #
##   None                                                                   #
#############################################################################

#' Boxtplots for Overall Comparisons by Continous Variables.
#'
#' @description
#' This function creates a PDF file containing box plots for the continuous
#' variables in the provided data. If the number of columns in `data` exceeds
#' `bin.size`, the function splits the plots across multiple pages.
#'
#' @param data A matrix or data frame containing the raw data to be plotted.
#' @param pdf_title A string representing the name of the PDF file to
#' be created.
#' @param bin_size An integer specifying the maximum number of box plots to
#' display on a single page.
#' @param mf_row A numeric vector of length two specifying the layout
#'  (rows and columns) for the plots on each page.
#' @param y_lim An optional numeric vector defining the y-axis limits
#' for the plots.
#'
#' @return A PDF file containing the box plots for the continuous variables.
#'
#' @examples
#' # Loading data
#' data.df <- cytodata
#' # Generate box plots for log2-transformed values to check for outliers:
#' cyt_bp(log2(data.df[, -c(1:4)]), pdf_title = "boxplot_by_cytokine_log2.pdf")
#'
#' @export
cyt_bp <- function(data, pdf_title, bin_size = 25,
                   mf_row = c(1, 1), y_lim = NULL) {
  pdf(file = pdf_title)
  par(mfrow = mf_row, cex.axis = 0.75)
  n_col <- ncol(data)
  if (n_col > bin_size) {
    for (i in 1:ceiling(n_col / bin_size)) {
      the_cols <- ((i - 1) * bin_size + 1):min(i * bin_size, n_col)
      if (is.null(y_lim)) {
        boxplot(data[, the_cols], las = 2, cex = 0.75)
      } else {
        boxplot(data[, the_cols], las = 2, cex = 0.75, ylim = y_lim)
      }
    }
  } else {
    if (is.null(y_lim)) {
      boxplot(data, las = 2, cex = 0.75)
    } else {
      boxplot(data, las = 2, cex = 0.75, ylim = y_lim)
    }
  }
  # Remove the ".pdf" from the pdf_title for display purposes
  display_title <- sub("\\.pdf$", "", pdf_title, ignore.case = TRUE)
  title(main = display_title, cex.main = 1.5, line = 2.5)

  dev.off()
}
