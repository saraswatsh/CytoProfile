##*********************************************************************************
## function to draw an error-bar plot enriched with p-value and/or effect size    #
##   for the comparison of each of selected groups to the baseline group          #
## Author: Xiaohua Douglas Zhang, January 2022                                    #
## Arguments:                                                                     #
##   center.df: a data frame containing the following columns for each group:     #
##            "name" for group names                                              #
##            "center" for mean or median,                                        #
##            "spread" for standard deviation, MAD or standard error,             #
##            "p.value" for p-value,                                              #
##            "effect.size" for effect size based on SSMD                         #
##            note: the first column of center.df must be for the baseline group  #
##   pLab: whether to label the p-value                                           #
##   esLab: whether to label the effect size usually in SSMD                      #
## Output:                                                                        #
##   None                                                                         #
##*********************************************************************************

#' Error-bar plot for comparison.
#'
#' @param center.df A data frame containing the following columns for each group:
#'                  "name" for group names
#'                  "center" for mean or median,
#'                  "spread" for standard deviation, MAD or standard error,
#'                  "p.value" for p-value,
#'                  "effect.size" for effect size based on SSMD.
#' @param pLab Whether to label the p-value , default set to TRUE.
#' @param esLab whether to label the effect size usually in SSMD, default set to TRUE.
#' @param xlab Labeling x axis.
#' @param ylab Labeling y axis.
#' @param main Title of the graph.
#' @description This function draws an error-bar plot enriched with p-value and/or effect size
#'              for the comparison of each selected groups to the baseline group.
#'              Note: The first column for center.df must be for the baseline group.
#' @return Prints error-bar plots.
#' @export
#'
cyt.errbp <- function(center.df, pLab=TRUE, esLab=TRUE, classSymbol=TRUE,
                      xlab="", ylab="", main="") {
  # Significance mark function
  significanceMark.fn <- function(p.value) {
    ##**************************************************************************
    ## function to mark the significant level based on p-value
    ## Input
    ##   p.value:  a vector for p-value
    ## Author: Xiaohua Douglas Zhang, 2022
    ## Reference:
    ##   Zhang XHD, 2011. Optimal High-Throughput Screening: Practical
    ##   Experimental Design and Data Analysis for Genome-scale RNAi Research.
    ##   Cambridge University Press, Cambridge, UK
    ##**************************************************************************
    significanceClass.vec <- rep(NA, length(p.value))
    significanceClass.vec[p.value <= 0.05 & p.value > 0.01] <- "*"
    significanceClass.vec[p.value <= 0.01 & p.value > 0.001] <- "**"
    significanceClass.vec[p.value <= 0.001 & p.value > 0.0001] <- "***"
    significanceClass.vec[p.value <= 0.0001 & p.value > 0.00001] <- "****"
    significanceClass.vec[p.value <= 0.00001] <- "*****"
    if(sum(is.na(p.value)) > 0) significanceClass.vec[is.na(p.value)] <- NA
    significanceClass.vec
  }

  # Effect size mark function
  effectSizeMark.fn <- function(SSMD) {
    ##**************************************************************************
    ## function to mark SSMD-based effect size being large or beyond
    ## Input
    ##   SSMD:  a vector for SSMD value
    ## Author: Xiaohua Douglas Zhang, 2022
    ## Reference:
    ##   Zhang XHD, 2011. Optimal High-Throughput Screening: Practical
    ##   Experimental Design and Data Analysis for Genome-scale RNAi Research.
    ##   Cambridge University Press, Cambridge, UK
    ##**************************************************************************
    effectClass.vec <- rep(NA, length(SSMD))
    effectClass.vec[SSMD >= 5] <- ">>>>>"
    effectClass.vec[SSMD >= 3 & SSMD < 5] <- ">>>>"
    effectClass.vec[SSMD >= 1.645 & SSMD < 3] <- ">>>"
    effectClass.vec[SSMD >= 1 & SSMD < 1.645] <- ">>"
    effectClass.vec[SSMD > 0.25 & SSMD < 1] <- ">"
    effectClass.vec[SSMD < 0.25 & SSMD > -0.25] <- " "
    effectClass.vec[SSMD <= -0.25 & SSMD > -1] <- "<"
    effectClass.vec[SSMD <= -1 & SSMD > -1.645] <- "<<"
    effectClass.vec[SSMD <= -1.645 & SSMD > -3] <- "<<<"
    effectClass.vec[SSMD <= -3 & SSMD > -5] <- "<<<<"
    effectClass.vec[SSMD <= -5] <- "<<<<<"
    if(sum(is.na(SSMD)) > 0) effectClass.vec[is.na(SSMD)] <- NA
    effectClass.vec
  }

  Center <- center.df[, "center"]
  names(Center) <- center.df[, "name"]
  Spread <- center.df[, "spread"]
  mid <- barplot(Center, plot=FALSE)
  pvalue <- center.df[-1, "p.value"]
  effect.size <- center.df[-1, "effect.size"]

  # Plot barplot
  yLim0 <- c(min(0, min(Center-Spread)), max(0, max(Center+Spread))) * 1.4
  yRange <- diff(yLim0)
  yLim <- range(c(yLim0, Center[-1]+sign(Center[-1])*(Spread[-1]+yRange/8)))
  barplot(Center, ylim = yLim, las=2, xlab=xlab, ylab=ylab, main=main)

  # Add error bars
  arrows(x0=mid, y0=Center, x1=mid, y1=Center+sign(Center)*Spread, code=3, angle=90, length=0.1)

  # Add text for p-value
  if(pLab) {
    if(classSymbol) {
      text(x=mid[-1], y=Center[-1]+sign(Center[-1])*(Spread[-1]+yRange/20),
           labels=significanceMark.fn(pvalue), cex=1)
    } else {
      text(x=mid[-1], y=Center[-1]+sign(Center[-1])*(Spread[-1]+yRange/20),
           labels=paste("p=", ifelse(pvalue>0.001, round(pvalue, 3),
                                     formatC(pvalue, format="e", digits=1)), sep=""),
           cex=0.75)
    }
  }

  # Add text for effect size
  if(esLab) {
    if(classSymbol) {
      text(x=mid[-1], y=Center[-1]+sign(Center[-1])*(Spread[-1]+yRange/8),
           labels=effectSizeMark.fn(effect.size), cex=1)
    } else {
      text(x=mid[-1], y=Center[-1]+sign(Center[-1])*(Spread[-1]+yRange/8),
           labels=paste("s=", round(effect.size, 3), sep=""), cex=0.75)
    }
  }
}
