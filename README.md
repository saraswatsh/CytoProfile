
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
[GitHub](https://github.com/steal123/CytProfile) with:

``` r
# install.packages("devtools")
devtools::install_github("steal123/CytProfile")
```

## Example

This is a basic example which shows you how to analyze a data set:

``` r
library(CytProfile)
#> Loading required package: BiocManager
#> Bioconductor version '3.17' is out-of-date; the current release version '3.18'
#>   is available with R version '4.3'; see https://bioconductor.org/install
#> Loading required package: mixOmics
#> Loading required package: MASS
#> Loading required package: lattice
#> Loading required package: ggplot2
#> 
#> Loaded mixOmics 6.24.0
#> Thank you for using mixOmics!
#> Tutorials: http://mixomics.org
#> Bookdown vignette: https://mixomicsteam.github.io/Bookdown
#> Questions, issues: Follow the prompts at http://mixomics.org/contact-us
#> Cite us:  citation('mixOmics')
#> Loading required package: moments
#> Loading required package: tidyverse
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.4
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ lubridate 1.9.3     ✔ tibble    3.2.1
#> ✔ purrr     1.0.2     ✔ tidyr     1.3.0
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ✖ purrr::map()    masks mixOmics::map()
#> ✖ dplyr::select() masks MASS::select()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
## basic example code
# Loading in data
data("cytdata.df")
data.df = cytdata.df

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
cyt.skku(data.df[,-c(1,3)], Title = "Skew and Kurtosis.pdf")
#> png 
#>   2
# Performing ANOVA comparisons test for univariate analysis
cyt.anova(data.df[,c(2:4,5:6)]) # This only considers 2 cytokines for this example only
#> $IL.17F_Group
#> ObeseNGT-LeanNGT   preT2D-LeanNGT      T2D-LeanNGT  preT2D-ObeseNGT 
#>        0.7032340        0.8858652        0.4831059        0.2876549 
#>     T2D-ObeseNGT       T2D-preT2D 
#>        0.9868528        0.1478753 
#> 
#> $GM.CSF_Group
#> ObeseNGT-LeanNGT   preT2D-LeanNGT      T2D-LeanNGT  preT2D-ObeseNGT 
#>      0.812501426      0.747090933      0.001642992      0.999884541 
#>     T2D-ObeseNGT       T2D-preT2D 
#>      0.037723063      0.031770058 
#> 
#> $IL.17F_Treatment
#>   Etomoxir-BPTES     UK5099-BPTES    Vehicle-BPTES  UK5099-Etomoxir 
#>     6.112997e-06     4.769375e-08     1.078340e-01     7.632088e-01 
#> Vehicle-Etomoxir   Vehicle-UK5099 
#>     2.709334e-02     1.070054e-03 
#> 
#> $GM.CSF_Treatment
#>   Etomoxir-BPTES     UK5099-BPTES    Vehicle-BPTES  UK5099-Etomoxir 
#>     5.479244e-09     9.588992e-09     4.995809e-02     9.988347e-01 
#> Vehicle-Etomoxir   Vehicle-UK5099 
#>     6.937026e-04     1.095577e-03

## Partial Least Squares Discriminant Analysis (PLS-DA) 
# In this code, we will have background predict to be turned on to see the classification areas and 
# we will also print out the confusion matrix based on classification. 
# Note this takes into account all groups and treatment. 
x.df = data.df[,-c(1,4)]
cyt.plsda(x.df, title = "Example PLS-DA Analysis.pdf", bg = TRUE, conf.mat = TRUE)
#> Warning in cor(A[[k]], variates.A[[k]]): the standard deviation is zero
#> Confusion Matrix for PLS-DA Comparison 
#>          predicted.as.LeanNGT predicted.as.ObeseNGT predicted.as.preT2D
#> LeanNGT                     3                     0                  25
#> ObeseNGT                    4                     0                  18
#> preT2D                      3                     0                  25
#> T2D                         2                     0                   9
#>          predicted.as.T2D
#> LeanNGT                15
#> ObeseNGT               14
#> preT2D                 16
#> T2D                    25
#> Confusion Matrix for PLS-DA Comparison with VIP Score > 1 
#>          predicted.as.LeanNGT predicted.as.ObeseNGT predicted.as.preT2D
#> LeanNGT                     1                     0                  25
#> ObeseNGT                    1                     0                  19
#> preT2D                      0                     0                  30
#> T2D                         1                     0                   8
#>          predicted.as.T2D
#> LeanNGT                17
#> ObeseNGT               16
#> preT2D                 14
#> T2D                    27

# Filtering data for specific groups and treatment
filt.data = filter(data.df, Group != "LeanNGT" & Group != "preT2D")
cyt.plsda(filt.data[,-c(1,3)], colors = c("black", "purple"), title = "Example PLS-DA Analysis 2.pdf", bg = TRUE, conf.mat = TRUE)
#> Warning in cor(A[[k]], variates.A[[k]]): the standard deviation is zero

#> Warning in cor(A[[k]], variates.A[[k]]): the standard deviation is zero

#> Warning in cor(A[[k]], variates.A[[k]]): the standard deviation is zero

#> Warning in cor(A[[k]], variates.A[[k]]): the standard deviation is zero
#> Confusion Matrix for PLS-DA Comparison 
#>          predicted.as.ObeseNGT predicted.as.T2D
#> ObeseNGT                     5                4
#> T2D                          3                6
#> Confusion Matrix for PLS-DA Comparison with VIP Score > 1 
#>          predicted.as.ObeseNGT predicted.as.T2D
#> ObeseNGT                     5                4
#> T2D                          2                7
```
