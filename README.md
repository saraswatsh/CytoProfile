
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CytProfile

<!-- badges: start -->
<!-- badges: end -->

The goal of CytProfile is to conduct quality control using biological
meaningful cutoff on raw measured values of cytokines. Test on
distributional symmetry to suggest the adopt of transformation.
Exploratory analysis including summary statistics, enriched boxplot, and
barplots. Univariate analysis and Multivariate analysis options are
further available to dive deeper

## Installation

Before installation of the CytProfile package, make sure to install
BiocManager and mix0mics packages using:

``` r
## install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
## install mixOmics 
BiocManager::install('mixOmics')
```

You can install the development version of CytProfile from
[GitHub](https://github.com/saraswatsh/CytProfile) with:

``` r
# install.packages("devtools")
devtools::install_github("saraswatsh/CytProfile")
```

## Example

This is a basic example which shows you how to analyze a data set. The
PDF files of the plots and PLS-DA analysis will be under the output
folder.

``` r
library(CytProfile)
## basic example code
# Loading in data
data("cytdata.df")
data.df = cytdata.df
## Setting working directory to output folder to save the PDF files. 
setwd("C:/Users/shubh/Desktop/RA/R Package/CytProfile/output")
## Exploratory Data Analysis
# Generating boxplots to check for outliers for raw values
cyt.bp(data.df[,-c(1:4)], Title = "Boxplot.byCytokine.Raw.pdf") # We are removing the first 4 columns as we only want the continuous variables. 
#> png 
#>   2
# Generating boxplots to check for outliers for log2 values
cyt.bp(log2(data.df[,-c(1:4)]), Title = "Boxplot.byCytokine.log2.pdf") # Make sure to use log2 to transform the cytokines and same reason as above for removing initial columns.
#> png 
#>   2
# Generating histograms for skewness and kurtosis based on raw values and log2 transformation
cyt.skku(data.df[,-c(1,4)], Title = "Skew and Kurtosis.pdf")
# Generating Error Bar Plot
cytokine.mat = cytdata.df[, -c(1:4)] # Extracting all cytokines to be stored in one object
cytokineNames = colnames(cytokine.mat) # Extracting the cytokine names
nCytokine = length(cytokineNames) # Obtaining the total number of cytokines
results = cyt.skku(cytdata.df[,-c(1,4)], printResLog = TRUE) # Extracting values
#> [1] "Results for Log2 Transformed Values:/n"
pdf( "barErrorPlot.pdf" )
par(mfrow=c(2,2), mar=c(8.1,  4.1, 4.1, 2.1) )
for( k in 1:nCytokine ) {
  center.df = data.frame( "name"=rownames(results[,,k]), results[,,k] )
  cyt.errbp(center.df, pLab=FALSE, esLab=FALSE, classSymbol=TRUE,
               ylab="Concentration in log2 scale",  main=cytokineNames[k] )
}
dev.off()
#> png 
#>   2
# Performing ANOVA comparisons test for univariate analysis
cyt.anova(data.df[,c(2:3,5:6)]) # This only considers 2 cytokines for this example only
#> $IL.17F_Group
#>  PreT2D-ND     T2D-ND T2D-PreT2D 
#>  0.6445189  0.1638310  0.6223573 
#> 
#> $GM.CSF_Group
#>  PreT2D-ND     T2D-ND T2D-PreT2D 
#>  0.7730980  0.5373287  0.1893654 
#> 
#> $IL.17F_Treatment
#>          LPS-CD3/CD28 Unstimulated-CD3/CD28      Unstimulated-LPS 
#>          7.229772e-13          7.214229e-13          9.990797e-01 
#> 
#> $GM.CSF_Treatment
#>          LPS-CD3/CD28 Unstimulated-CD3/CD28      Unstimulated-LPS 
#>          7.183143e-13          6.974421e-13          3.481621e-01
## Partial Least Squares Discriminant Analysis (PLS-DA) 
# In this code, we will have background predict to be turned on to see the classification areas and 
# we will also print out the confusion matrix based on classification. 
# Note this takes into account all groups and treatment. 
x.df = data.df[,-c(1,4)]
cyt.plsda(x.df, title = "Example PLS-DA Analysis.pdf", bg = TRUE, conf.mat = TRUE)
#> Confusion Matrix for PLS-DA Comparison 
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  17               16
#> PreT2D               0                  31                2
#> T2D                  0                  12               21
#> Confusion Matrix for PLS-DA Comparison with VIP Score > 1 
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  14               19
#> PreT2D               0                  28                5
#> T2D                  0                  10               23

# Filtering data for specific groups and treatment
filt.data = filter(data.df, Group != "ND", Treatment != "Unstimulated")
cyt.plsda(filt.data[,-c(1,4)], colors = c("black", "purple"), title = "Example PLS-DA Analysis 2.pdf", bg = TRUE, conf.mat = TRUE)
#> Confusion Matrix for PLS-DA Comparison 
#>        predicted.as.PreT2D predicted.as.T2D
#> PreT2D                  17               16
#> T2D                      6               27
#> Confusion Matrix for PLS-DA Comparison with VIP Score > 1 
#>        predicted.as.PreT2D predicted.as.T2D
#> PreT2D                  17               16
#> T2D                      3               30
```
