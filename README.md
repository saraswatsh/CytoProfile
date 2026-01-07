
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CytoProfile <a href="https://cytoprofile.cytokineprofile.org/"><img src="man/figures/logo.png" align="right" height="133" alt="CytoProfile website" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/saraswatsh/CytoProfile/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/saraswatsh/CytoProfile/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/CytoProfile)](https://CRAN.R-project.org/package=CytoProfile)
[![CRAN
checks](https://badges.cranchecks.info/summary/CytoProfile.svg)](https://cran.r-project.org/web/checks/check_results_CytoProfile.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/CytoProfile)](https://cran.r-project.org/package=CytoProfile)
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
[GitHub](https://github.com/saraswatsh/CytoProfile/tree/devel) with:

``` r
# install.packages("devtools")
devtools::install_github("saraswatsh/CytoProfile", ref = "devel")
```

Install CytoProfile from
[CRAN](https://cran.r-project.org/package=CytoProfile) with:

``` r
install.packages("CytoProfile")
```

See change log for the latest updates and changes at [News](NEWS.md)

To look at the vignettes included in the package, use:

``` r
browseVignettes("CytoProfile")
```

[Vignettes](https://cytoprofile.cytokineprofile.org/articles/index.html)
are also available on CytoProfile website to learn how to use the
package.

For more details on the package, please visit the [CytoProfile
website](https://cytoprofile.cytokineprofile.org/).
