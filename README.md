
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CytoProfile

<!-- badges: start -->
<!-- badges: end -->

The goal of CytoProfile is to conduct quality control using biological
meaningful cutoff on raw measured values of cytokines. Specifically,
test on distributional symmetry to suggest the adopt of transformation.
Conduct exploratory analysis including summary statistics, generate
enriched barplots, and boxplots. Further, conduct univariate analysis
and multivariate analysis for advance analysis.

## Installation

Before installation of the CytoProfile package, make sure to install
BiocManager and mix0mics packages using:

``` r
## install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
## install mixOmics 
BiocManager::install('mixOmics')
```

You can install the development version of CytoProfile from
[GitHub](https://github.com/saraswatsh/CytoProfile) with:

``` r
# install.packages("devtools")
devtools::install_github("saraswatsh/CytoProfile")
```

## Example

This is a basic example which shows you how to analyze a data set. The
PDF files of the plots and PLS-DA analysis will be under the output
folder.

``` r
library(CytoProfile)
## basic example code
# Loading in data
data("cytdata.df")
data.df = cytdata.df
## Setting working directory to output folder to save the PDF files. 
setwd("C:/Users/shubh/Desktop/RA/R Package/CytoProfile/output")
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
cyt.plsda(x.df, title = "Example PLS-DA Analysis.pdf", bg = TRUE, conf.mat = TRUE, scale = "log2",
var.num = 25, cv.opt = "loocv", comp.num = 2, colors = c("black", "purple", "red2"), 
pch.values = c(16,4,3), style = NULL)
#> [1] "Results based on log2 transformation:"
#> [1] "CD3/CD28 T2D vs PreT2D LOOCV Accuracy: 42"
#> [1] "CD3/CD28 T2DvsPreT2D LOOCV Accuracy (VIP) Cytokines: 48"
#> [1] "CD3/CD28 T2D vs PreT2D Confusion Matrix for PLS-DA Comparison"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  15               18
#> PreT2D               0                  19               14
#> T2D                  0                   4               29
#> [1] "CD3/CD28 T2D vs PreT2D Confusion Matrix for PLS-DA Comparison with VIP Score > 1"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  15               18
#> PreT2D               0                  19               14
#> T2D                  0                   4               29
#> [1] "LPS T2D vs PreT2D LOOCV Accuracy: 38"
#> [1] "LPS T2DvsPreT2D LOOCV Accuracy (VIP) Cytokines: 45"
#> [1] "LPS T2D vs PreT2D Confusion Matrix for PLS-DA Comparison"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  10               23
#> PreT2D               0                  15               18
#> T2D                  0                   5               28
#> [1] "LPS T2D vs PreT2D Confusion Matrix for PLS-DA Comparison with VIP Score > 1"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  10               23
#> PreT2D               0                  16               17
#> T2D                  0                   3               30
#> [1] "Unstimulated T2D vs PreT2D LOOCV Accuracy: 37"
#> [1] "Unstimulated T2DvsPreT2D LOOCV Accuracy (VIP) Cytokines: 47"
#> [1] "Unstimulated T2D vs PreT2D Confusion Matrix for PLS-DA Comparison"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  17               16
#> PreT2D               0                  31                2
#> T2D                  0                  12               21
#> [1] "Unstimulated T2D vs PreT2D Confusion Matrix for PLS-DA Comparison with VIP Score > 1"
#>        predicted.as.ND predicted.as.PreT2D predicted.as.T2D
#> ND                   0                  14               19
#> PreT2D               0                  28                5
#> T2D                  0                  10               23
#> png 
#>   2

# Filtering data for specific groups and treatment
filt.data = filter(data.df, Group != "ND", Treatment != "Unstimulated")
cyt.plsda(filt.data[,-c(1,4)], colors = c("black", "purple"), 
          title = "Example PLS-DA Analysis 2.pdf", bg = TRUE, scale = "log2", 
          conf.mat = TRUE, var.num = 25, 
          cv.opt = "Mfold", fold.num = 5, 
          comp.num = 3, pch.values = c(3,4), style = "3d")
#> [1] "Results based on log2 transformation:"
#> [1] "CD3/CD28 T2D vs PreT2D Mfold Accuracy: 68"
#> [1] "CD3/CD28 T2DvsPreT2D Mfold Accuracy (VIP) Cytokines: 73"
#> [1] "CD3/CD28 T2D vs PreT2D Confusion Matrix for PLS-DA Comparison"
#>        predicted.as.PreT2D predicted.as.T2D
#> PreT2D                  19               14
#> T2D                      4               29
#> [1] "CD3/CD28 T2D vs PreT2D Confusion Matrix for PLS-DA Comparison with VIP Score > 1"
#>        predicted.as.PreT2D predicted.as.T2D
#> PreT2D                  19               14
#> T2D                      4               29
#> [1] "LPS T2D vs PreT2D Mfold Accuracy: 64"
#> [1] "LPS T2DvsPreT2D Mfold Accuracy (VIP) Cytokines: 69"
#> [1] "LPS T2D vs PreT2D Confusion Matrix for PLS-DA Comparison"
#>        predicted.as.PreT2D predicted.as.T2D
#> PreT2D                  17               16
#> T2D                      6               27
#> [1] "LPS T2D vs PreT2D Confusion Matrix for PLS-DA Comparison with VIP Score > 1"
#>        predicted.as.PreT2D predicted.as.T2D
#> PreT2D                  17               16
#> T2D                      3               30
#> png 
#>   2

# Generating Volcano Plot
data.df = cytdata.df[,-4]
volc_plot = cyt.volc(data.df, "Group", cond1 = "T2D", cond2 = "ND", fold_change_thresh = 2.0, top_labels= 15)
#>                    Cytokine      FC_Log      P_Log Significant
#> IL.12.P70         IL.12.P70 -2.60117683 2.18641971        TRUE
#> IL.6                   IL.6 -0.95013174 3.94758527       FALSE
#> IL.27                 IL.27 -0.67878724 2.33099419       FALSE
#> IL.23                 IL.23 -0.87320747 1.95290632       FALSE
#> CCL.20.MIP.3A CCL.20.MIP.3A -0.48569948 1.40917287       FALSE
#> IL.2                   IL.2 -0.80577278 1.22848122       FALSE
#> IL.17F               IL.17F -0.93024059 1.16938373       FALSE
#> IL.10                 IL.10 -0.48121242 1.01734902       FALSE
#> IL.28A               IL.28A -0.31081278 0.98351262       FALSE
#> IL.17A               IL.17A -0.80415853 0.90173665       FALSE
#> IL.1B                 IL.1B -0.61564856 0.83381951       FALSE
#> GM.CSF               GM.CSF  0.45980342 0.62042612       FALSE
#> IL.21                 IL.21 -0.62254771 0.51843946       FALSE
#> IL.17E.IL.25   IL.17E.IL.25  0.01449957 0.49515782       FALSE
#> IL.22                 IL.22 -0.30363695 0.47550506       FALSE
#> IL.9                   IL.9 -0.32752255 0.43117675       FALSE
#> TNF.A                 TNF.A -0.15647551 0.21142412       FALSE
#> IL.31                 IL.31  0.21056699 0.20929529       FALSE
#> IL.4                   IL.4  0.21161574 0.20542291       FALSE
#> IL.5                   IL.5 -0.20808037 0.17512546       FALSE
#> IL.15                 IL.15 -0.05298748 0.12055764       FALSE
#> IL.13                 IL.13 -0.07717527 0.06654004       FALSE
#> IFN.G                 IFN.G -0.09088794 0.06221451       FALSE
#> TNF.B                 TNF.B  0.07037796 0.05224667       FALSE
#> IL.33                 IL.33  0.01213249 0.01719622       FALSE
ggsave("VolcanoPlot.png", plot = volc_plot$`T2D vs ND`, dpi = 300)

# Printing table (This is usually printed by default when the function is called and not saved as an object)
print(volc_plot$`T2D vs ND`$data)
#>                    Cytokine      FC_Log      P_Log Significant         Label
#> IL.12.P70         IL.12.P70 -2.60117683 2.18641971        TRUE     IL.12.P70
#> IL.6                   IL.6 -0.95013174 3.94758527       FALSE          IL.6
#> IL.27                 IL.27 -0.67878724 2.33099419       FALSE         IL.27
#> IL.23                 IL.23 -0.87320747 1.95290632       FALSE         IL.23
#> CCL.20.MIP.3A CCL.20.MIP.3A -0.48569948 1.40917287       FALSE CCL.20.MIP.3A
#> IL.2                   IL.2 -0.80577278 1.22848122       FALSE          IL.2
#> IL.17F               IL.17F -0.93024059 1.16938373       FALSE        IL.17F
#> IL.10                 IL.10 -0.48121242 1.01734902       FALSE         IL.10
#> IL.28A               IL.28A -0.31081278 0.98351262       FALSE        IL.28A
#> IL.17A               IL.17A -0.80415853 0.90173665       FALSE        IL.17A
#> IL.1B                 IL.1B -0.61564856 0.83381951       FALSE         IL.1B
#> GM.CSF               GM.CSF  0.45980342 0.62042612       FALSE        GM.CSF
#> IL.21                 IL.21 -0.62254771 0.51843946       FALSE         IL.21
#> IL.17E.IL.25   IL.17E.IL.25  0.01449957 0.49515782       FALSE  IL.17E.IL.25
#> IL.22                 IL.22 -0.30363695 0.47550506       FALSE         IL.22
#> IL.9                   IL.9 -0.32752255 0.43117675       FALSE              
#> TNF.A                 TNF.A -0.15647551 0.21142412       FALSE              
#> IL.31                 IL.31  0.21056699 0.20929529       FALSE              
#> IL.4                   IL.4  0.21161574 0.20542291       FALSE              
#> IL.5                   IL.5 -0.20808037 0.17512546       FALSE              
#> IL.15                 IL.15 -0.05298748 0.12055764       FALSE              
#> IL.13                 IL.13 -0.07717527 0.06654004       FALSE              
#> IFN.G                 IFN.G -0.09088794 0.06221451       FALSE              
#> TNF.B                 TNF.B  0.07037796 0.05224667       FALSE              
#> IL.33                 IL.33  0.01213249 0.01719622       FALSE

# Generating Heat map
cyt.heatmap(data = data.df,
                    scale = "log2",        # Optional scaling
                    annotation_col_name = "Group",
                    title = "Heatmap.png")
#> png 
#>   2

# Generating dual flashlights plot
data.df = cytdata.df[,-c(1,3:4)]

dfp = cyt.dualflashplot(data.df, group_var = "Group", group1 = "T2D", group2 = "ND", 
                  ssmd_thresh = -0.2, log2fc_thresh = 1, top_labels = 10)
#> # A tibble: 25 × 11
#>    cytokine         mean_ND mean_PreT2D   mean_T2D variance_ND variance_PreT2D
#>    <chr>              <dbl>       <dbl>      <dbl>       <dbl>           <dbl>
#>  1 CCL.20.MIP.3A   634.        404.       887.        6.72e+ 5         2.74e+5
#>  2 GM.CSF            2.65        3.11       1.92      2.63e+ 1         3.14e+1
#>  3 IFN.G         57730.      18303.     61484.        2.86e+10         2.30e+9
#>  4 IL.10           979.        836.      1366.        1.99e+ 6         1.19e+6
#>  5 IL.12.P70        13.0        39.1       78.9       4.15e+ 2         2.56e+4
#>  6 IL.13          1064.       1543.      1122.        5.60e+ 6         1.11e+7
#>  7 IL.15             7.92        4.29       8.22      3.54e+ 1         2.58e+1
#>  8 IL.17A          352.        653.       615.        9.40e+ 5         2.88e+6
#>  9 IL.17E.IL.25      0.0101      0.0163     0.01      1.01e- 6         3.88e-3
#> 10 IL.17F            1.63        2.35       3.11      1.56e+ 1         3.37e+1
#> 11 IL.1B          2806.       2977.      4299.        6.63e+ 7         3.76e+7
#> 12 IL.2           9227.      10718.     16129.        2.60e+ 8         4.10e+8
#> 13 IL.21           205.        210.       316.        3.15e+ 5         2.49e+5
#> 14 IL.22             0.0513      0.0684     0.0633    4.58e- 3         4.51e-3
#> 15 IL.23             0.147       0.243      0.269     3.13e- 2         9.37e-2
#> 16 IL.27             0.0662      0.0834     0.106     6.18e- 3         5.66e-3
#> 17 IL.28A            0.0537      0.0710     0.0666    2.45e- 3         5.10e-3
#> 18 IL.31             0.0409      0.0905     0.0354    6.62e- 3         4.88e-2
#> 19 IL.33             1.17        1.43       1.16      2.09e+ 0         2.71e+0
#> 20 IL.4              0.344       0.707      0.297     4.24e- 1         2.96e+0
#> 21 IL.5            134.        340.       155.        1.09e+ 5         9.88e+5
#> 22 IL.6           4620.       5197.      8925.        2.86e+ 7         5.72e+7
#> 23 IL.9            203.        256.       254.        1.34e+ 5         2.11e+5
#> 24 TNF.A          5046.       3069.      5624.        7.02e+ 7         1.63e+7
#> 25 TNF.B             0.641       0.709      0.610     2.37e+ 0         2.76e+0
#> # ℹ 5 more variables: variance_T2D <dbl>, ssmd <dbl>, log2FC <dbl>,
#> #   SSMD_Category <chr>, Significant <lgl>
ggsave("DualFlashlightPlot.png", plot = dfp$plot_env$p, dpi = 300, width = 3000, height = 2000, units = "px")

# Printing table (This is usually printed by default when the function is called and not saved as an object)
print(dfp$data)
#> # A tibble: 25 × 11
#>    cytokine         mean_ND mean_PreT2D mean_T2D variance_ND variance_PreT2D
#>    <chr>              <dbl>       <dbl>    <dbl>       <dbl>           <dbl>
#>  1 CCL.20.MIP.3A   634.        404.       887.      6.72e+ 5         2.74e+5
#>  2 GM.CSF            2.65        3.11       1.92    2.63e+ 1         3.14e+1
#>  3 IFN.G         57730.      18303.     61484.      2.86e+10         2.30e+9
#>  4 IL.10           979.        836.      1366.      1.99e+ 6         1.19e+6
#>  5 IL.12.P70        13.0        39.1       78.9     4.15e+ 2         2.56e+4
#>  6 IL.13          1064.       1543.      1122.      5.60e+ 6         1.11e+7
#>  7 IL.15             7.92        4.29       8.22    3.54e+ 1         2.58e+1
#>  8 IL.17A          352.        653.       615.      9.40e+ 5         2.88e+6
#>  9 IL.17E.IL.25      0.0101      0.0163     0.01    1.01e- 6         3.88e-3
#> 10 IL.17F            1.63        2.35       3.11    1.56e+ 1         3.37e+1
#> # ℹ 15 more rows
#> # ℹ 5 more variables: variance_T2D <dbl>, ssmd <dbl>, log2FC <dbl>,
#> #   SSMD_Category <chr>, Significant <lgl>
```
