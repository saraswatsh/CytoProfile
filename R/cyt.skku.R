###############################################################################
############### Skewness and Kurtosis #########################################
###############################################################################

#' Distribution of the data set
#'
#' @param x.df A matrix or data frame of raw data.
#' @param Title Name for the PDF file.
#' @description
#' The function takes in a data frame and subsets the numeric columns from the data which is then
#' used to calculate the skewed and kurtosis values. The values are then plotted using histograms
#' to visualize the distribution of raw skewed values and log2 transformed values.
#'
#' @return Prints histograms of Skewness and Kurtosis of the continuous variables using raw data and log2 transformation.
#' @examples
#' data(cytdata.df)
#' cyt.skku(cytdata.df[,-c(1,3)], Title = "Skew and Kurtosis.pdf")
#' cyt.skku(cytdata.df[,-c(1,4)], Title = "Skew and Kurtosis.pdf")
#' @export
#'
cyt.skku = function(x.df, Title) {
  pdf(file = Title)
  cytokine.mat = x.df[, -c(1:2)]
  cytokineNames = colnames(cytokine.mat)
  nCytokine = length(cytokineNames)

  condt = !is.na(cytokine.mat) & cytokine.mat >0
  min(cytokine.mat[condt], na.rm=TRUE)  # [1] 0.01
  # quantile( cytokine.mat[condt], probs=0:100/100, na.rm=TRUE)
  Cutoff = min(cytokine.mat[condt], na.rm=TRUE)/10

  outcomes = c("n", "center", "spread", "skewness", "kurtosis")
  nOutcome = length( outcomes )

  ## raw value
  if ("Stimulation" %in% names(x.df[,c(1:2)])){
    Treatment.Group.vec = paste(x.df[,"Stimulation"], x.df[,"Group"], sep=".")
  }else{
    Treatment.Group.vec = paste(x.df[,"Treatment"], x.df[,"Group"], sep=".")
  }
  Treatment.Groups = names( tapply( x.df[,3], INDEX=Treatment.Group.vec, mean ))
  nTrtGroup = length( Treatment.Groups )
  result.arr = array(NA, dim=c(nTrtGroup, nOutcome, nCytokine) )
  dimnames(result.arr) = list( Treatment.Groups, outcomes, cytokineNames )

  for( k in 1:nCytokine ) {
    Y = cytokine.mat[,k]
    idx = Treatment.Group.vec
    result.arr[,1,k] = tapply( Y, INDEX=idx, function(x){sum(!is.na(x))} )
    result.arr[,2,k] = tapply( Y, INDEX=idx, mean, na.rm=TRUE)
    #  result.arr[,3,k] = tapply( Y, INDEX=idx, sd, na.rm=TRUE)       # standard deviation
    result.arr[,3,k] = tapply( Y, INDEX=idx,
                               function(x){sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))} ) # standard error
    result.arr[,4,k] = tapply( Y, INDEX=idx, skewness, na.rm=TRUE)  # skewness
    result.arr[,5,k] = tapply( Y, INDEX=idx, kurtosis, na.rm=TRUE)  #kurtosis
  }
  result.noComb.raw.arr = result.arr

  ## Log2 value
  if ("Stimulation" %in% names(x.df[,c(1:2)])){
    Treatment.Group.vec = paste(x.df[,"Stimulation"], x.df[,"Group"], sep=".")
  }else{
    Treatment.Group.vec = paste(x.df[,"Treatment"], x.df[,"Group"], sep=".")
  }
  Treatment.Groups = names( tapply( x.df[,3], INDEX=Treatment.Group.vec, mean ))
  nTrtGroup = length( Treatment.Groups )
  result.arr = array(NA, dim=c(nTrtGroup, nOutcome, nCytokine) )
  dimnames(result.arr) = list( Treatment.Groups, outcomes, cytokineNames )

  for( k in 1:nCytokine ) {
    Y = log2(cytokine.mat[,k]+Cutoff)
    idx = Treatment.Group.vec
    result.arr[,1,k] = tapply( Y, INDEX=idx, function(x){sum(!is.na(x))} )
    result.arr[,2,k] = tapply( Y, INDEX=idx, mean, na.rm=TRUE)
    #  result.arr[,3,k] = tapply( Y, INDEX=idx, sd, na.rm=TRUE)       # standard deviation
    result.arr[,3,k] = tapply( Y, INDEX=idx,
                               function(x){sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))} ) # standard error
    result.arr[,4,k] = tapply( Y, INDEX=idx, skewness, na.rm=TRUE)  # skewness
    result.arr[,5,k] = tapply( Y, INDEX=idx, kurtosis, na.rm=TRUE)  #kurtosis
  }
  result.noComb.arr = result.arr

  par(mfrow=c(1,2))
  hist(result.noComb.raw.arr[,"skewness",], xlab="skewness of raw data (every combination)", main="raw data: skewness")
  hist(result.noComb.arr[,"skewness",], xlab="skewness of log2 data (every combination)", main="log2 transformed: skewness")

  hist(result.noComb.raw.arr[,"kurtosis",], xlab="kurtosis of raw data", main="raw data: kurtosis")
  hist(result.noComb.arr[,"kurtosis",], xlab="kurtosis of log2 data", main="log2 transformed:kurtosis")

  dev.off()
}
