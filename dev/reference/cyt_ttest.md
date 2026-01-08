# Two Sample T-test Comparisons. **\[deprecated\]**

This function performs pairwise comparisons between two groups for each
combination of a categorical predictor (with exactly two levels) and a
continuous outcome variable. It first converts any character variables
in `data` to factors and, if specified, applies a log2 transformation to
the continuous variables. Depending on the value of `scale`, the
function conducts either a two-sample t-test (if `scale = "log2"`) or a
Mann-Whitney U test (if `scale` is `NULL`). The resulting p-values are
printed and returned.

## Usage

``` r
cyt_ttest(data, scale = NULL, verbose = TRUE, format_output = FALSE)
```

## Arguments

- data:

  A matrix or data frame containing continuous and categorical
  variables.

- scale:

  A character specifying a transformation for continuous variables.
  Options are `NULL` (default) and `"log2"`. When `scale = "log2"`, a
  log2 transformation is applied and a two-sample t-test is used; when
  `scale` is `NULL`, a Mann-Whitney U test is performed.

- verbose:

  A logical indicating whether to print the p-values of the statistical
  tests. Default is `TRUE`.

- format_output:

  Logical. If TRUE, returns the results as a tidy data frame. Default is
  `FALSE`.

## Value

If `format_output` is FALSE, returns a list of p-values (named by
Outcome and Categorical variable). If TRUE, returns a data frame in a
tidy format.

## Author

Shubh Saraswat

## Examples

``` r
data_df <- ExampleData1[, -c(3)]
data_df <- dplyr::filter(data_df, Group != "ND", Treatment != "Unstimulated")
# Test example
cyt_ttest(
 data_df[, c(1:2, 5:6)],
 scale = "log2",
 verbose = TRUE,
 format_output = TRUE
)
#> Warning: `cyt_ttest()` was deprecated in CytoProfile 0.2.3.
#> ℹ Please use `cyt_univariate()` instead.
#> $results
#>   Outcome Categorical      Comparison
#> 1   IFN.G       Group   PreT2D vs T2D
#> 2   IL.10       Group   PreT2D vs T2D
#> 3   IFN.G   Treatment CD3/CD28 vs LPS
#> 4   IL.10   Treatment CD3/CD28 vs LPS
#>                                                Test Estimate Statistic P_value
#> 1 Wilcoxon rank sum test with continuity correction   -2.463    1599.0   0.008
#> 2 Wilcoxon rank sum test with continuity correction   -0.956    1625.0   0.012
#> 3 Wilcoxon rank sum test with continuity correction    9.024    4132.5   0.000
#> 4 Wilcoxon rank sum test with continuity correction    1.690    3091.0   0.000
#> 
```
