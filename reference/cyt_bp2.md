# Boxplot Function Enhanced for Specific Group Comparisons.

This function generates a PDF file containing boxplots for each
combination of numeric and factor variables in the provided data. It
first converts any character columns to factors and checks that the data
contains at least one numeric and one factor column. If the scale
argument is set to "log2", all numeric columns are log2-transformed. The
function then creates boxplots using ggplot2 for each numeric variable
grouped by each factor variable.

## Usage

``` r
cyt_bp2(data, pdf_title, scale = NULL, y_lim = NULL)
```

## Arguments

- data:

  A matrix or data frame of raw data.

- pdf_title:

  A string representing the title (and filename) of the PDF file. If
  `NULL`, the boxplots are displayed on the current graphics device.
  Defaults to `NULL`.

- scale:

  Transformation option for continuous variables. Options are NULL
  (default) and "log2". When set to "log2", numeric columns are
  transformed using the log2 function.

- y_lim:

  An optional numeric vector defining the y-axis limits for the plots.

## Value

A PDF file containing the boxplots.

## Author

Shubh Saraswat

## Examples

``` r
# Loading data
data_df <- ExampleData1[, -c(3, 5:28)]
data_df <- dplyr::filter(data_df, Group == "T2D", Treatment == "Unstimulated")
cyt_bp2(data_df, pdf_title = NULL, scale = "log2")
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the CytoProfile package.
#>   Please report the issue at
#>   <https://github.com/saraswatsh/CytoProfile/issues>.


```
