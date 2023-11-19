
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
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.
