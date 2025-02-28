#################################################################################
## Function to generate a PDF file to show the measured value by column of x.df #
## Author: Xiaohua Douglas Zhang, January 2022                                  #
## Arguments:                                                                   #
##   x.df: a matrix of measured continuous value with row and column names      #
##   Title: names for the PDF file                                              #
##   bin.size: the number of boxplots shown in a single pane                    #
##   yLim: a range for the y-axis to be plotted                                 #
## Output:                                                                      #
##   None                                                                       #
#################################################################################

#' Generating a PDF file to show the measured value by column of the data frame.
#' @param x.df A matrix or data frame containing the raw data to be plotted.
#' @param Title A string representing the name of the PDF file.
#' @param bin.size The maximum number of box plots to display on a single page.
#' @param mfRow A numeric vector of length two specifying the layout (rows and columns) for the plots on each page.
#' @param yLim An optional numeric vector defining the y-axis limits for the plots.
#' @description
#' This function creates a PDF file containing box plots for the continuous variables in the provided data.
#' If the number of columns in `x.df` exceeds `bin.size`, the function splits the plots across multiple pages.
#' @return A PDF file containing the box plots of the continuous variables.
#' @export
cyt.bp <- function(x.df, Title, bin.size=25, mfRow=c(1,1), yLim=NULL) {
  pdf(file=Title)
  par(mfrow=mfRow, cex.axis=0.75)
  nCol <- ncol(x.df)
  if(nCol > bin.size) {
    for(i in 1:ceiling(nCol/bin.size)) {
      theCols <- ((i-1)*bin.size+1):min(i*bin.size, nCol)
      if(is.null(yLim)) {
        boxplot(x.df[, theCols], las=2, cex=0.75)
      } else {
        boxplot(x.df[, theCols], las=2, cex=0.75, ylim=yLim)
      }
    }
  } else {
    if(is.null(yLim)) {
      boxplot(x.df, las=2, cex=0.75)
    } else {
      boxplot(x.df, las=2, cex=0.75, ylim=yLim)
    }
  }
  # Remove the ".pdf" from the Title for display purposes
  displayTitle <- sub("\\.pdf$", "", Title, ignore.case = TRUE)
  title(main=displayTitle, cex.main=1.5, line=2.5)

  dev.off()
}

