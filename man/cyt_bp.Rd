% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cyt_bp.R
\name{cyt_bp}
\alias{cyt_bp}
\title{Boxtplots for Overall Comparisons by Continous Variables.}
\usage{
cyt_bp(data, pdf_title, bin_size = 25, mf_row = c(1, 1), y_lim = NULL)
}
\arguments{
\item{data}{A matrix or data frame containing the raw data to be plotted.}

\item{pdf_title}{A string representing the name of the PDF file to
be created.}

\item{bin_size}{An integer specifying the maximum number of box plots to
display on a single page.}

\item{mf_row}{A numeric vector of length two specifying the layout
(rows and columns) for the plots on each page.}

\item{y_lim}{An optional numeric vector defining the y-axis limits
for the plots.}
}
\value{
A PDF file containing the box plots for the continuous variables.
}
\description{
This function creates a PDF file containing box plots for the continuous
variables in the provided data. If the number of columns in \code{data} exceeds
\code{bin.size}, the function splits the plots across multiple pages.
}
\examples{
# Loading data
data.df <- cytodata
# Generate box plots for log2-transformed values to check for outliers:
cyt_bp(log2(data.df[, -c(1:4)]), pdf_title = "boxplot_by_cytokine_log2.pdf")

}
