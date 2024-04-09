
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
# Generating Error Bar Plot enriched with p-value and effect size 
data.df = cytdata.df[,-1]
cyt.mat = log2(data.df[,-c(1:3)])
data.df1 = data.frame(data.df[,c(1:3)], cyt.mat)
cytokineNames = colnames(cyt.mat)
nCytokine = length(cytokineNames)
condt = !is.na(cyt.mat) & cyt.mat >0
Cutoff = min(cyt.mat[condt], na.rm=TRUE)/10
# Creating a matrix for p-values from anova tests
p.aov.mat = matrix(NA, nrow=nCytokine, ncol=3)
# Changing column names
dimnames(p.aov.mat) = list( cytokineNames, c("Group", "Treatment", "Interaction") )
# Matrix to extract p-values from Tukey group comparison
p.groupComp.mat = matrix(NA, nrow=nCytokine, ncol=3)
# Changing column names
dimnames(p.groupComp.mat) = list( cytokineNames, c("2-1", "3-1", "3-2") )
# Matrix for SSMD same size as other matrices
ssmd.groupComp.stm.mat = mD.groupComp.stm.mat = p.groupComp.stm.mat = p.groupComp.mat

for( i in 1:nCytokine ) {
   #i = 1 # i=2
  Cytokine = (cyt.mat[,i]+Cutoff)
  cytokine.aov = aov( Cytokine ~ Group * Treatment, data=data.df)
  aov.table = summary(cytokine.aov)[[1]]
  p.aov.mat[i,] = aov.table[1:3,5]
  p.groupComp.mat[i,] = TukeyHSD(cytokine.aov)$Group[1:3,4]
  p.groupComp.stm.mat[i,] = TukeyHSD(cytokine.aov)$`Group:Treatment`[c(1:3),4]
  mD.groupComp.stm.mat[i,] = TukeyHSD(cytokine.aov)$`Group:Treatment`[c(1:3),1]
  ssmd.groupComp.stm.mat[i,]=mD.groupComp.stm.mat[i,]/sqrt(2*aov.table["Residuals","Mean Sq"])
}

# p.aov.mat
# p.groupComp.mat
# p.groupComp.stm.mat
results = cyt.skku(cytdata.df[,-c(1,4)], printResLog = TRUE)
#> [1] "Results for Log2 Transformed Values:/n"
pdf( "barErrorPlot.enriched.pdf" )
par(mfrow=c(2,3), mar=c(8.1,  4.1, 4.1, 2.1) )
for( k in 1:nCytokine ) {
  #k = 1
  result.mat = results[1:9,,k]
  center.df =
    data.frame( "name"=rownames(result.mat), result.mat[, c("center", "spread")],
                "p.value"= c(1,p.groupComp.stm.mat[k,1:2]),
                "effect.size"=c(0,ssmd.groupComp.stm.mat[k,1:2])
    )
  cyt.errbp(center.df, pLab=TRUE, esLab=TRUE, classSymbol=TRUE,
               ylab="Concentration in log2 scale", main=cytokineNames[k])
}
dev.off()
#> png 
#>   2

# Performing ANOVA comparisons test for univariate analysis
cyt.anova(data.df[,c(2:3,5:6)]) # This only considers 2 cytokines for this example only
#> $Time_Treatment
#>          LPS-CD3/CD28 Unstimulated-CD3/CD28      Unstimulated-LPS 
#>                     1                     1                     1 
#> 
#> $GM.CSF_Treatment
#>          LPS-CD3/CD28 Unstimulated-CD3/CD28      Unstimulated-LPS 
#>          7.183143e-13          6.974421e-13          3.481621e-01 
#> 
#> $IFN.G_Treatment
#>          LPS-CD3/CD28 Unstimulated-CD3/CD28      Unstimulated-LPS 
#>          7.445156e-13          7.394085e-13          9.987624e-01
## Partial Least Squares Discriminant Analysis (PLS-DA) 
# In this code, we will have background predict to be turned on to see the classification areas and 
# we will also print out the confusion matrix based on classification. 
# Note this takes into account all groups and treatment and all values are log transformed through 
# cyt.plsda function. 
data.df = cytdata.df
x.df = data.df[,-c(1,4)]
cyt.plsda(x.df, title = "Example PLS-DA Analysis.pdf", bg = TRUE, conf.mat = TRUE, var.num = 25, cv.opt = "loocv")
#> [1] "CD3/CD28 T2DvsPreT2D LOOCV Accuracy: 0.424242424242424"
#> [1] "CD3/CD28 T2DvsPreT2D LOOCV Accuracy (VIP) Cytokines: 0.484848484848485"
#> [1] "CD3/CD28 T2DvsPreT2D Confusion Matrix for PLS-DA Comparison"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  15               18
#> PreT2D               0                  19               14
#> T2D                  0                   4               29
#> [1] "CD3/CD28 T2DvsPreT2D Confusion Matrix for PLS-DA Comparison with VIP Score > 1"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  15               18
#> PreT2D               0                  19               14
#> T2D                  0                   4               29
#> [1] "LPS T2DvsPreT2D LOOCV Accuracy: 0.383838383838384"
#> [1] "LPS T2DvsPreT2D LOOCV Accuracy (VIP) Cytokines: 0.454545454545455"
#> [1] "LPS T2DvsPreT2D Confusion Matrix for PLS-DA Comparison"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  10               23
#> PreT2D               0                  15               18
#> T2D                  0                   5               28
#> [1] "LPS T2DvsPreT2D Confusion Matrix for PLS-DA Comparison with VIP Score > 1"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  10               23
#> PreT2D               0                  16               17
#> T2D                  0                   3               30
#> [1] "Unstimulated T2DvsPreT2D LOOCV Accuracy: 0.373737373737374"
#> [1] "Unstimulated T2DvsPreT2D LOOCV Accuracy (VIP) Cytokines: 0.474747474747475"
#> [1] "Unstimulated T2DvsPreT2D Confusion Matrix for PLS-DA Comparison"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  17               16
#> PreT2D               0                  31                2
#> T2D                  0                  12               21
#> [1] "Unstimulated T2DvsPreT2D Confusion Matrix for PLS-DA Comparison with VIP Score > 1"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  14               19
#> PreT2D               0                  28                5
#> T2D                  0                  10               23
#> png 
#>   2

# Filtering data for specific groups and treatment
filt.data = filter(data.df, Group != "ND", Treatment != "Unstimulated")
cyt.plsda(filt.data[,-c(1,4)], colors = c("black", "purple"), title = "Example PLS-DA Analysis 2.pdf", bg = TRUE, conf.mat = TRUE, var.num = 25, cv.opt = "Mfold", fold.num = 5)
#> [1] "CD3/CD28 T2DvsPreT2D Mfold Accuracy: 0.683560606060606"
#> [1] "CD3/CD28 T2DvsPreT2D Mfold Accuracy (VIP) Cytokines: 0.726030303030303"
#> [1] "CD3/CD28 T2DvsPreT2D Confusion Matrix for PLS-DA Comparison"
#>        predicted.as.PreT2D predicted.as.T2D
#> PreT2D                  19               14
#> T2D                      4               29
#> [1] "CD3/CD28 T2DvsPreT2D Confusion Matrix for PLS-DA Comparison with VIP Score > 1"
#>        predicted.as.PreT2D predicted.as.T2D
#> PreT2D                  19               14
#> T2D                      4               29
#> [1] "LPS T2DvsPreT2D Mfold Accuracy: 0.643378787878788"
#> [1] "LPS T2DvsPreT2D Mfold Accuracy (VIP) Cytokines: 0.694424242424242"
#> [1] "LPS T2DvsPreT2D Confusion Matrix for PLS-DA Comparison"
#>        predicted.as.PreT2D predicted.as.T2D
#> PreT2D                  17               16
#> T2D                      6               27
#> [1] "LPS T2DvsPreT2D Confusion Matrix for PLS-DA Comparison with VIP Score > 1"
#>        predicted.as.PreT2D predicted.as.T2D
#> PreT2D                  17               16
#> T2D                      3               30
#> png 
#>   2
```
