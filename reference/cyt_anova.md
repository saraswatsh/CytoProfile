# ANOVA Analysis on Continuous Variables.

This function performs an analysis of variance (ANOVA) for each
continuous variable against every categorical predictor in the input
data. Character columns are automatically converted to factors; all
factor columns are used as predictors while numeric columns are used as
continuous outcomes. For each valid predictor (i.e., with more than one
level and no more than 10 levels), Tukey's Honest Significant Difference
(HSD) test is conducted and the adjusted p-values for pairwise
comparisons are extracted.

## Usage

``` r
cyt_anova(data, format_output = FALSE)
```

## Arguments

- data:

  A data frame or matrix containing both categorical and continuous
  variables. Character columns will be converted to factors and used as
  predictors, while numeric columns will be used as continuous outcomes.

- format_output:

  Logical. If TRUE, returns the results as a tidy data frame instead of
  a list. Default is FALSE.

## Value

If `format_output` is FALSE (default), a list of adjusted p-values from
Tukey's HSD tests for each combination of continuous outcome and
categorical predictor. List elements are named in the format
"Outcome_Categorical". If `format_output` is TRUE, a data frame in a
tidy format.

## Author

Shubh Saraswat

## Examples

``` r
data("ExampleData1")
cyt_anova(ExampleData1[, c(1:2, 5:6)], format_output = TRUE)
#>                        Outcome Categorical            Comparison  P_adj
#> PreT2D-ND               GM.CSF       Group             PreT2D-ND 0.7731
#> T2D-ND                  GM.CSF       Group                T2D-ND 0.5373
#> T2D-PreT2D              GM.CSF       Group            T2D-PreT2D 0.1894
#> PreT2D-ND1               IFN.G       Group             PreT2D-ND 0.0883
#> T2D-ND1                  IFN.G       Group                T2D-ND 0.9779
#> T2D-PreT2D1              IFN.G       Group            T2D-PreT2D 0.0550
#> LPS-CD3/CD28            GM.CSF   Treatment          LPS-CD3/CD28 0.0000
#> Unstimulated-CD3/CD28   GM.CSF   Treatment Unstimulated-CD3/CD28 0.0000
#> Unstimulated-LPS        GM.CSF   Treatment      Unstimulated-LPS 0.3482
#> LPS-CD3/CD281            IFN.G   Treatment          LPS-CD3/CD28 0.0000
#> Unstimulated-CD3/CD281   IFN.G   Treatment Unstimulated-CD3/CD28 0.0000
#> Unstimulated-LPS1        IFN.G   Treatment      Unstimulated-LPS 0.9988
```
