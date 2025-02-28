###############################################################################
############### Skewness and Kurtosis #########################################
###############################################################################

#' Distribution of the data set shown by skewness and kurtosis
#'
#' @param x.df A matrix or data frame of raw data. The first two columns should contain grouping variables
#' (e.g., "Treatment" and "Group" or "Stimulation" and "Group"), while the remaining columns are assumed to be numeric cytokine measurements.
#' @param Title A character string specifying the name for the PDF file. If provided, the histograms will be saved to this PDF file;
#' if not, the plots are produced on the current graphics device.
#' @param printResRaw Logical. If \code{TRUE}, the function prints and returns the computed summary array (count, mean, standard error, skewness, and kurtosis)
#' for the raw data. Default is \code{FALSE}.
#' @param printResLog Logical. If \code{TRUE}, the function prints and returns the computed summary array for the log2-transformed data. Default is \code{FALSE}.
#'
#' @description
#' This function subsets the numeric columns from the input data (excluding the first two grouping columns) and computes summary statistics
#' (including count, central tendency, spread as standard error, skewness, and kurtosis) for each group defined by a combination of the first two columns.
#' A small cutoff (one-tenth of the minimum positive value in the data) is added to each numeric value before applying the log2 transformation
#' to handle non-positive values. Histograms are then generated to visualize the distribution of skewness and kurtosis for both the raw and log2-transformed data.
#'
#' @return The function prints histograms of skewness and kurtosis for both raw and log2-transformed data.
#' Optionally, if \code{printResRaw} and/or \code{printResLog} are \code{TRUE}, the function returns the corresponding summary arrays.
#'
#' @examples
#' \dontrun{
#' data(cytodata)
#' cyt.skku(cytodata[,-c(1,4)], Title = "Skew and Kurtosis.pdf")
#' }
#'
#' @export
#' @import e1071

cyt.skku <- function(x.df, Title = NULL, printResRaw = FALSE, printResLog = FALSE) {
  if(!is.null(Title)){
    pdf(file = Title)
    cytokine.mat <- x.df[, -c(1:2)]
    cytokineNames <- colnames(cytokine.mat)
    nCytokine <- length(cytokineNames)

    condt <- !is.na(cytokine.mat) & cytokine.mat > 0
    min(cytokine.mat[condt], na.rm=TRUE)  # [1] 0.01
    # quantile( cytokine.mat[condt], probs=0:100/100, na.rm=TRUE)
    Cutoff <- min(cytokine.mat[condt], na.rm=TRUE)/10

    outcomes <- c("n", "center", "spread", "skewness", "kurtosis")
    nOutcome <- length(outcomes)

    ## raw value
    if ("Stimulation" %in% names(x.df[,c(1:2)])){
      Treatment.Group.vec <- paste(x.df[,"Stimulation"], x.df[,"Group"], sep=".")
    } else {
      Treatment.Group.vec <- paste(x.df[,"Treatment"], x.df[,"Group"], sep=".")
    }
    Treatment.Groups <- names(tapply(x.df[,3], INDEX=Treatment.Group.vec, mean))
    nTrtGroup <- length(Treatment.Groups)
    result.arr <- array(NA, dim=c(nTrtGroup, nOutcome, nCytokine))
    dimnames(result.arr) <- list(Treatment.Groups, outcomes, cytokineNames)

    for(k in 1:nCytokine) {
      Y <- cytokine.mat[,k]
      idx <- Treatment.Group.vec
      result.arr[,1,k] <- tapply(Y, INDEX=idx, function(x){sum(!is.na(x))})
      result.arr[,2,k] <- tapply(Y, INDEX=idx, mean, na.rm=TRUE)
      #  result.arr[,3,k] = tapply( Y, INDEX=idx, sd, na.rm=TRUE)       # standard deviation
      result.arr[,3,k] <- tapply(Y, INDEX=idx,
                                 function(x){sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))} ) # standard error
      result.arr[,4,k] <- tapply(Y, INDEX=idx, skewness, na.rm=TRUE)  # skewness
      result.arr[,5,k] <- tapply(Y, INDEX=idx, kurtosis, na.rm=TRUE)  # kurtosis
    }
    result.noComb.raw.arr <- result.arr

    ## Log2 value
    if ("Stimulation" %in% names(x.df[,c(1:2)])){
      Treatment.Group.vec <- paste(x.df[,"Stimulation"], x.df[,"Group"], sep=".")
    } else {
      Treatment.Group.vec <- paste(x.df[,"Treatment"], x.df[,"Group"], sep=".")
    }
    Treatment.Groups <- names(tapply(x.df[,3], INDEX=Treatment.Group.vec, mean))
    nTrtGroup <- length(Treatment.Groups)
    result.arr <- array(NA, dim=c(nTrtGroup, nOutcome, nCytokine))
    dimnames(result.arr) <- list(Treatment.Groups, outcomes, cytokineNames)

    for(k in 1:nCytokine) {
      Y <- log2(cytokine.mat[,k] + Cutoff)
      idx <- Treatment.Group.vec
      result.arr[,1,k] <- tapply(Y, INDEX=idx, function(x){sum(!is.na(x))})
      result.arr[,2,k] <- tapply(Y, INDEX=idx, mean, na.rm=TRUE)
      #  result.arr[,3,k] = tapply( Y, INDEX=idx, sd, na.rm=TRUE)       # standard deviation
      result.arr[,3,k] <- tapply(Y, INDEX=idx,
                                 function(x){sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))} ) # standard error
      result.arr[,4,k] <- tapply(Y, INDEX=idx, skewness, na.rm=TRUE)  # skewness
      result.arr[,5,k] <- tapply(Y, INDEX=idx, kurtosis, na.rm=TRUE)  # kurtosis
    }
    result.noComb.arr <- result.arr

    par(mfrow=c(1,2))
    hist(result.noComb.raw.arr[,"skewness",], xlab="skewness of raw data", main="raw data: skewness")
    hist(result.noComb.arr[,"skewness",], xlab="skewness of log2 data", main="log2 transformed: skewness")

    hist(result.noComb.raw.arr[,"kurtosis",], xlab="kurtosis of raw data", main="raw data: kurtosis")
    hist(result.noComb.arr[,"kurtosis",], xlab="kurtosis of log2 data", main="log2 transformed:kurtosis")

    dev.off()

    if(printResRaw == TRUE & printResLog == TRUE){
      print("Results for Raw Values:/n")
      return(result.noComb.raw.arr)
      print("Results for Log2 Transformed Values:/n")
      return(result.noComb.arr)
    } else if (printResRaw == TRUE){
      print("Results for Raw Values:/n")
      return(result.noComb.raw.arr)
    } else if (printResLog == TRUE){
      print("Results for Log2 Transformed Values:/n")
      return(result.noComb.arr)
    }
  }
  else {
    cytokine.mat <- x.df[, -c(1:2)]
    cytokineNames <- colnames(cytokine.mat)
    nCytokine <- length(cytokineNames)

    condt <- !is.na(cytokine.mat) & cytokine.mat > 0
    min(cytokine.mat[condt], na.rm=TRUE)  # [1] 0.01
    # quantile( cytokine.mat[condt], probs=0:100/100, na.rm=TRUE)
    Cutoff <- min(cytokine.mat[condt], na.rm=TRUE)/10

    outcomes <- c("n", "center", "spread", "skewness", "kurtosis")
    nOutcome <- length(outcomes)

    ## raw value
    if ("Stimulation" %in% names(x.df[,c(1:2)])){
      Treatment.Group.vec <- paste(x.df[,"Stimulation"], x.df[,"Group"], sep=".")
    } else {
      Treatment.Group.vec <- paste(x.df[,"Treatment"], x.df[,"Group"], sep=".")
    }
    Treatment.Groups <- names(tapply(x.df[,3], INDEX=Treatment.Group.vec, mean))
    nTrtGroup <- length(Treatment.Groups)
    result.arr <- array(NA, dim=c(nTrtGroup, nOutcome, nCytokine))
    dimnames(result.arr) <- list(Treatment.Groups, outcomes, cytokineNames)

    for(k in 1:nCytokine) {
      Y <- cytokine.mat[,k]
      idx <- Treatment.Group.vec
      result.arr[,1,k] <- tapply(Y, INDEX=idx, function(x){sum(!is.na(x))})
      result.arr[,2,k] <- tapply(Y, INDEX=idx, mean, na.rm=TRUE)
      #  result.arr[,3,k] = tapply( Y, INDEX=idx, sd, na.rm=TRUE)       # standard deviation
      result.arr[,3,k] <- tapply(Y, INDEX=idx,
                                 function(x){sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))} ) # standard error
      result.arr[,4,k] <- tapply(Y, INDEX=idx, skewness, na.rm=TRUE)  # skewness
      result.arr[,5,k] <- tapply(Y, INDEX=idx, kurtosis, na.rm=TRUE)  # kurtosis
    }
    result.noComb.raw.arr <- result.arr

    ## Log2 value
    if ("Stimulation" %in% names(x.df[,c(1:2)])){
      Treatment.Group.vec <- paste(x.df[,"Stimulation"], x.df[,"Group"], sep=".")
    } else {
      Treatment.Group.vec <- paste(x.df[,"Treatment"], x.df[,"Group"], sep=".")
    }
    Treatment.Groups <- names(tapply(x.df[,3], INDEX=Treatment.Group.vec, mean))
    nTrtGroup <- length(Treatment.Groups)
    result.arr <- array(NA, dim=c(nTrtGroup, nOutcome, nCytokine))
    dimnames(result.arr) <- list(Treatment.Groups, outcomes, cytokineNames)

    for(k in 1:nCytokine) {
      Y <- log2(cytokine.mat[,k] + Cutoff)
      idx <- Treatment.Group.vec
      result.arr[,1,k] <- tapply(Y, INDEX=idx, function(x){sum(!is.na(x))})
      result.arr[,2,k] <- tapply(Y, INDEX=idx, mean, na.rm=TRUE)
      #  result.arr[,3,k] = tapply( Y, INDEX=idx, sd, na.rm=TRUE)       # standard deviation
      result.arr[,3,k] <- tapply(Y, INDEX=idx,
                                 function(x){sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))} ) # standard error
      result.arr[,4,k] <- tapply(Y, INDEX=idx, skewness, na.rm=TRUE)  # skewness
      result.arr[,5,k] <- tapply(Y, INDEX=idx, kurtosis, na.rm=TRUE)  # kurtosis
    }
    result.noComb.arr <- result.arr

    if(printResRaw == TRUE & printResLog == TRUE){
      print("Results for Raw Values:/n")
      return(result.noComb.raw.arr)
      print("Results for Log2 Transformed Values:/n")
      return(result.noComb.arr)
    } else if (printResRaw == TRUE){
      print("Results for Raw Values:/n")
      return(result.noComb.raw.arr)
    } else if (printResLog == TRUE){
      print("Results for Log2 Transformed Values:/n")
      return(result.noComb.arr)
    }
  }
}

