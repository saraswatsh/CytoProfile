
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CytoProfile <a href="https://cytoprofile.cytokineprofile.org/"><img src="man/figures/logo.png" align="right" height="133" alt="CytoProfile website" /></a>

<!-- badges: start -->

**CRAN Build Status:** [![R-CMD-check Release
Build](https://github.com/saraswatsh/CytoProfile/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/saraswatsh/CytoProfile/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/CytoProfile)](https://CRAN.R-project.org/package=CytoProfile)
[![CRAN
checks](https://badges.cranchecks.info/summary/CytoProfile.svg)](https://cran.r-project.org/web/checks/check_results_CytoProfile.html)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/CytoProfile)](https://cran.r-project.org/package=CytoProfile)
[![Last Commit
Release](https://img.shields.io/github/last-commit/saraswatsh/CytoProfile/main)](https://github.com/saraswatsh/CytoProfile/commits/main/)

**Development Build Status:** [![R-CMD-check Dev
Build](https://github.com/saraswatsh/CytoProfile/actions/workflows/R-CMD-check.yaml/badge.svg?branch=devel)](https://github.com/saraswatsh/CytoProfile/actions/workflows/R-CMD-check.yaml)
[![Dev
Version](https://img.shields.io/badge/devel%20version-0.2.3.9000-red)](https://github.com/saraswatsh/CytoProfile/tree/devel)
[![Last Commit
Dev](https://img.shields.io/github/last-commit/saraswatsh/CytoProfile/devel)](https://github.com/saraswatsh/CytoProfile/commits/devel/)

[![CodeFactor](https://www.codefactor.io/repository/github/saraswatsh/cytoprofile/badge/devel)](https://www.codefactor.io/repository/github/saraswatsh/cytoprofile/overview/devel)

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

See change log for the latest updates and changes on release build at
[News](https://cytoprofile.cytokineprofile.org/news/index.html) and
development build at
[News](https://cytoprofile.cytokineprofile.org/dev/news/index.html).

To look at the vignettes included in the package, use:

``` r
browseVignettes("CytoProfile")
```

[Vignettes](https://cytoprofile.cytokineprofile.org//dev/articles/index.html)
are also available on CytoProfile website to learn how to use the
package.

For more details on the released build of the package, please visit the
[CytoProfile website](https://cytoprofile.cytokineprofile.org/). For the
development version of the package, please visit the [Development
CytoProfile website](https://cytoprofile.cytokineprofile.org/dev).
