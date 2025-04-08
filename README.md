
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CytoProfile <a href="https://saraswatsh.github.io/CytoProfile/"><img src="man/figures/logo.png" align="right" height="150" alt="CytoProfile website" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/saraswatsh/CytoProfile/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/saraswatsh/CytoProfile/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/CytoProfile)](https://CRAN.R-project.org/package=CytoProfile)
[![CRAN
checks](https://badges.cranchecks.info/summary/CytoProfile.svg)](https://cran.r-project.org/web/checks/check_results_CytoProfile.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/CytoProfile)](https://cran.r-project.org/package=CytoProfile)
[![License: GPL (\>=
2)](https://img.shields.io/badge/license-GPL%20(%3E=%202)-blue.svg)](https://cran.r-project.org/web/licenses/GPL%20(%3E=%202))
[![Last
Commit](https://img.shields.io/github/last-commit/saraswatsh/CytoProfile.svg)](https://github.com/saraswatsh/CytoProfile/commits/main)
<!-- badges: end -->

The goal of CytoProfile is to conduct quality control using biological
meaningful cutoff on raw measured values of cytokines. Specifically,
test on distributional symmetry to suggest the adopt of transformation.
Conduct exploratory analysis including summary statistics, generate
enriched barplots, and boxplots. Further, conduct univariate analysis
and multivariate analysis for advance analysis.

## Installation

Before installation of the CytoProfile package, make sure to install
BiocManager and mixOmics packages using:

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

Install CytoProfile from CRAN with:

``` r
install.packages("CytoProfile")
```

See change log for the latest updates and changes at [News](NEWS.md)

## Example

Below are examples of using the functions provided in CytoProfile. Any
saved or generated files that are PDF or PNG format will be found at in
the [Output
Folder](https://github.com/saraswatsh/CytoProfile/tree/main/output).

## 1. Data Loading and set up

``` r
# Loading all packages required
# Data manipulation and reshaping
library(dplyr)       # For data filtering, grouping, and summarising.
library(tidyr)       # For reshaping data (e.g., pivot_longer, pivot_wider).

# Plotting and visualization
library(ggplot2)     # For creating all the ggplot-based visualizations.
library(gridExtra)   # For arranging multiple plots on a single page.
library(ggrepel)     # For improved label placement in plots (e.g., volcano plots).
library(gplots)      # For heatmap.2, which is used to generate heatmaps.
library(plot3D)      # For creating 3D scatter plots in PCA and sPLS-DA analyses.
library(reshape2)    # For data transformation (e.g., melt) in cross-validation plots.

# Statistical analysis
library(mixOmics)    # For multivariate analyses (PCA, sPLS-DA, etc.).
library(e1071)       # For computing skewness and kurtosis.
library(pROC)        # For ROC curve generation in machine learning model evaluation.

# Machine learning
library(xgboost)     # For building XGBoost classification models.
library(randomForest) # For building Random Forest classification models.
library(caret)       # For cross-validation and other machine learning utilities.

# Package development and document rendering
library(knitr)       # For knitting RMarkdown files and setting chunk options.
library(devtools)    # For installing the development version of the package from GitHub.

# Load in the CytoProfile package
library(CytoProfile)

# Loading in data
data("ExampleData1")
data_df <- ExampleData1

## Setting working directory to output folder to save the PDF files. 
opts_knit$set(root.dir = "E:/Desktop/RA/R Package/CytoProfile/output")
```

## 2. Exploratory Data Analysis

### Boxplots

``` r
# Generating boxplots to check for outliers for raw values
cyt_bp(data_df[, -c(1:3)], 
       pdf_title = "boxplot_by_cytokine_raw.pdf")  
# Removing the first 3 columns to retain only continuous variables.

# Generating boxplots to check for outliers for log2 values
cyt_bp(data_df[, -c(1:3)], 
       pdf_title = "boxplot_by_cytokine_log2.pdf",
       scale = "log2")  
# Using log2 transformation for cytokine values.
```

### Group-Specific Boxplots

``` r
# Raw values for group-specific boxplots
cyt_bp2(data_df[, -c(3)], 
        pdf_title = "boxplot_by_group_and_treatment_raw.pdf", 
        scale = NULL)

# Log2-transformed group-specific boxplots
cyt_bp2(data_df[, -c(3)], 
        pdf_title = "boxplot_by_group_and_treatment_log2.pdf", 
        scale = "log2")
```

## 3. Skewness and Kurtosis

``` r
# Histogram of skewness and kurtosis for raw data
cyt_skku(data_df[, -c(1:3)], 
         pdf_title = "skew_and_kurtosis.pdf", 
         group_cols = NULL)

# Histogram of skewness and kurtosis with grouping (e.g., "Group")
cyt_skku(ExampleData1[, -c(2:3)], 
         pdf_title = "skew_and_kurtosis_2.pdf", 
         group_cols = c("Group"))
```

## 4. Error Bar Plots

### Basic Error Bar Plot

``` r
# Generating basic error bar plots
cytokine_mat <- ExampleData1[, -c(1:3)]  # Extract all cytokines
cytokineNames <- colnames(cytokine_mat)  # Extract cytokine names
nCytokine <- length(cytokineNames)  # Total number of cytokines
results <- cyt_skku(ExampleData1[, -c(3)], print_res_log = TRUE, 
                    group_cols = c("Group", "Treatment"))
```

<img src="output/EDA 4-1.png" width="100%" />

``` r
pdf("bar_error_plot.pdf")

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2,2), mar = c(4, 4, 2, 1))

for (k in 1:nCytokine) {
  center_df <- data.frame(name = rownames(results[,,k]), results[,,k])
  cyt_errbp(center_df,
  p_lab = FALSE, es_lab = FALSE, class_symbol = TRUE,
  y_lab = "Concentration in log2 scale", main = cytokineNames[k])
}

par(oldpar)
if (dev.cur() > 1) dev.off()
#> png 
#>   2
```

### Enriched Error Bar Plot with p-values and Effect Sizes

``` r
# Generating Error Bar Plot enriched with p-value and effect size 
data_df <- ExampleData1[, -3]
cyt_mat <- log2(data_df[, -c(1:2)])
data_df1 <- data.frame(data_df[, 1:2], cyt_mat)
cytokineNames <- colnames(cyt_mat)
nCytokine <- length(cytokineNames)
condt <- !is.na(cyt_mat) & (cyt_mat > 0)
Cutoff <- min(cyt_mat[condt], na.rm = TRUE) / 10

# Create matrices for ANOVA and Tukey results
p_aov_mat <- matrix(NA, nrow = nCytokine, ncol = 3)
dimnames(p_aov_mat) <- list(cytokineNames, 
                            c("Group", "Treatment", "Interaction"))
p_groupComp_mat <- matrix(NA, nrow = nCytokine, ncol = 3)
dimnames(p_groupComp_mat) <- list(cytokineNames, 
                                  c("2-1", "3-1", "3-2"))
ssmd_groupComp_stm_mat <- mD_groupComp_stm_mat <- p_groupComp_stm_mat <- 
  p_groupComp_mat

for (i in 1:nCytokine) {
  Cytokine <- (cyt_mat[, i] + Cutoff)
  cytokine_aov <- aov(Cytokine ~ Group * Treatment, data = data_df)
  aov_table <- summary(cytokine_aov)[[1]]
  p_aov_mat[i, ] <- aov_table[1:3, 5]
  p_groupComp_mat[i, ] <- TukeyHSD(cytokine_aov)$Group[1:3, 4]
  p_groupComp_stm_mat[i, ] <- TukeyHSD(cytokine_aov)$`Group:Treatment`[1:3, 4]
  mD_groupComp_stm_mat[i, ] <- TukeyHSD(cytokine_aov)$`Group:Treatment`[1:3, 1]
  ssmd_groupComp_stm_mat[i, ] <- mD_groupComp_stm_mat[i, ] / 
    sqrt(2 * aov_table["Residuals", "Mean Sq"])
}

results <- cyt_skku(ExampleData1[, -c(3)], print_res_log = TRUE, 
                    group_cols = c("Group", "Treatment"))
```

<img src="output/EDA 5-1.png" width="100%" />

``` r
pdf("bar_error_plot_enriched.pdf")

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2,2), mar = c(3, 3, 2, 1))

par(mfrow = c(2,3), mar = c(8.1, 4.1, 4.1, 2.1))

for (k in 1:nCytokine) {
  result_mat <- results[1:9, , k]
  center_df <- data.frame(
    name = rownames(result_mat),
    result_mat[, c("center", "spread")],
    p.value = c(1, p_groupComp_stm_mat[k, 1:2]),
    effect.size = c(0, ssmd_groupComp_stm_mat[k, 1:2])
  )
  cyt_errbp(center_df, p_lab = TRUE, es_lab = TRUE, 
            class_symbol = TRUE,
            y_lab = "Concentration in log2 scale", 
            main = cytokineNames[k])
}
par(oldpar)
if (dev.cur() > 1) dev.off()
#> png 
#>   2
```

## 5. Univariate Analysis

### Two Sample T-test and Mann Whitney U Test

``` r
# Performing Two Sample T-test and Mann Whitney U Test
data_df <- ExampleData1[, -c(3)]
data_df <- dplyr::filter(data_df, Group != "ND", Treatment != "Unstimulated")
# Two sample T-test
cyt_ttest(data_df[, c(1:2, 5:6)], scale = "log2", verbose = TRUE, format_output = TRUE)
#>   Outcome Categorical      Comparison P_value
#> 1   IFN.G       Group   PreT2D vs T2D  0.0208
#> 2   IL.10       Group   PreT2D vs T2D  0.0248
#> 3   IFN.G   Treatment CD3/CD28 vs LPS  0.0000
#> 4   IL.10   Treatment CD3/CD28 vs LPS  0.0001
# Mann-Whitney U Test
cyt_ttest(data_df[, c(1:2, 5:6)], verbose = TRUE)
#> $IFN.G_Group
#> [1] 0.0085
#> 
#> $IL.10_Group
#> [1] 0.0119
#> 
#> $IFN.G_Treatment
#> [1] 0
#> 
#> $IL.10_Treatment
#> [1] 0
```

### ANOVA Comparisons Test

``` r
# Perform ANOVA comparisons test (example with 2 cytokines)
cyt_anova(data_df[, c(1:2, 5:6)], format_output = TRUE)
#> [1] Outcome     Categorical Comparison  P_adj      
#> <0 rows> (or 0-length row.names)
```

## 6. Multivariate Analysis

### Partial Least Squares Discriminant Analysis (PLS-DA)

``` r
# cyt_plsda function. 
data <- ExampleData1[, -c(3)]
data_df <- dplyr::filter(data, Group != "ND" & Treatment == "CD3/CD28")
cyt_splsda(data_df, 
           pdf_title = "example_spls_da_analysis.pdf", 
           colors = c("black", "purple"),
           bg = FALSE, scale = "log2", ellipse = TRUE,
           conf_mat = FALSE, var_num = 25, 
           cv_opt = "loocv", comp_num = 2, 
           pch_values = c(16, 4), group_col = "Group", group_col2 = "Treatment", 
           roc = TRUE)
#> png 
#>   2
```

## 7. Principal Component Analysis (PCA)

``` r
data <- ExampleData1[, -c(3,23)]
data_df <- filter(data, Group != "ND" & Treatment != "Unstimulated")
cyt_pca(data_df, 
        pdf_title = "example_pca_analysis.pdf", 
        colors = c("black", "red2"), 
        scale = "log2", 
        comp_num = 3, pch_values = c(16, 4), 
        style = "3D", group_col = "Group", group_col2 = "Treatment")
#> [1] "Results based on log2 transformation:"
#> png 
#>   2
cyt_pca(data_df, 
        pdf_title = "example_pca_analysis_2.pdf", 
        colors = c("black", "red2"), 
        scale = "log2", 
        comp_num = 2, pch_values = c(16, 4), 
        group_col = "Group")
#> [1] "Results based on log2 transformation:"
#> png 
#>   2
```

## 8. Volcano Plot

``` r
# Generating Volcano Plot
data_df <- ExampleData1[, -c(2:3)]
cyt_volc(data_df, group_col = "Group", 
                      cond1 = "T2D", cond2 = "ND", 
                      fold_change_thresh = 2.0, 
                      top_labels = 15)
#> $`T2D vs ND`
```

<img src="output/EDA 6-1.png" width="100%" />

## 9. Heatmap

``` r
# Generating Heat map
cyt_heatmap(data = data_df,
            scale = "log2",        # Optional scaling
            annotation_col_name = "Group",
            title = NULL)
```

<img src="output/EDA 7-1.png" width="100%" />

## 10. Dual Flashlight Plot

``` r
# Generating dual flashlights plot
data_df <- ExampleData1[, -c(2:3)]
dfp <- cyt_dualflashplot(data_df, group_var = "Group", 
                         group1 = "T2D", group2 = "ND", 
                         ssmd_thresh = -0.2, log2fc_thresh = 1, 
                         top_labels = 10)
# Print the plot
dfp
```

<img src="output/EDA 8-1.png" width="100%" />

``` r
# Print the table data used for plotting
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

## 11. Machine Learning Models

### Using XGBoost for classification

``` r
# Using XGBoost for classification
data_df0 <- ExampleData1
data_df <- data.frame(data_df0[, 1:3], log2(data_df0[, -c(1:3)]))
data_df <- data_df[, -c(2:3)]
data_df <- dplyr::filter(data_df, Group != "ND")

cyt_xgb(data = data_df, group_col = "Group",
                       nrounds = 500, max_depth = 4, eta = 0.05,
                       nfold = 5, cv = TRUE, eval_metric = "mlogloss",
                       early_stopping_rounds = NULL, top_n_features = 10,
                       verbose = 0, plot_roc = TRUE, print_results = FALSE)
```

<img src="output/ML1-1.png" width="100%" /><img src="output/ML1-2.png" width="100%" />

### Using Random Forest for classification

``` r
# Using Random Forest for classification
cyt_rf(data = data_df, group_col = "Group", k_folds = 5,
                     ntree = 1000, mtry = 4, run_rfcv = TRUE,
                     plot_roc = TRUE, verbose = FALSE)
```

<img src="output/ML2-1.png" width="100%" /><img src="output/ML2-2.png" width="100%" /><img src="output/ML2-3.png" width="100%" />
