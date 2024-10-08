% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cyt.skku.R
\name{cyt.skku}
\alias{cyt.skku}
\title{Distribution of the data set shown by skewness and kurtosis}
\usage{
cyt.skku(x.df, Title = NULL, printResRaw = FALSE, printResLog = FALSE)
}
\arguments{
\item{x.df}{A matrix or data frame of raw data.}

\item{Title}{Name for the PDF file.}

\item{printResRaw}{Print tibbles of skewness and kurtosis for raw values. Default set to FALSE.}

\item{printResLog}{Print tibbles of skewness and kurtosis for log2 transformed values. Default set to FALSE.}
}
\value{
Prints histograms of Skewness and Kurtosis of the continuous variables using raw data and log2 transformation.
}
\description{
The function takes in a data frame and subsets the numeric columns from the data which is then
used to calculate the skewed and kurtosis values. The values are then plotted using histograms
to visualize the distribution of raw skewed values and log2 transformed values.
}
\examples{
\dontrun{
data(cytodata)
cyt.skku(cytodata[,-c(1,4)], Title = "Skew and Kurtosis.pdf")
}
}
