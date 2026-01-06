# Boxplots for Continuous Variables with Optional Grouping

This function generates boxplots for numeric variables in a data frame
or matrix. It supports optional grouping by one or more categorical
variables. Numeric variables can be scaled using various transformations
before plotting. When grouping is not used, boxplots are arranged in
pages with a specified maximum number of plots per page. Plots can be
saved to a PDF file or displayed on the current graphics device.

## Usage

``` r
cyt_bp(
  data,
  pdf_title = NULL,
  group_by = NULL,
  bin_size = 25,
  y_lim = NULL,
  scale = c("none", "log2", "log10", "zscore", "custom"),
  custom_fn = NULL
)
```

## Arguments

- data:

  A matrix or data frame containing numeric and categorical variables.

- pdf_title:

  Optional string specifying the name of the PDF file to be created.
  When `NULL` (default), plots are drawn on the current graphics device.

- group_by:

  Optional character vector specifying one or more columns to use for
  grouping. If `NULL` (default) no grouping is applied.

- bin_size:

  Integer. Maximum number of boxplots per page when grouping is not
  used. Default is 25, as in the original `cyt_bp`.

- y_lim:

  Optional numeric vector giving y‑axis limits for the plots. Applies to
  all plots.

- scale:

  Character specifying a transformation for numeric variables. Accepts
  `"none"`, `"log2"`, `"log10"`, `"zscore"`, or `"custom"`. When
  `"custom"`, supply a function via `custom_fn`.

- custom_fn:

  A user supplied function to transform numeric columns when
  `scale = "custom"`.

## Value

Invisibly returns a list of `ggplot` objects. When `pdf_title` is
provided, plots are written to the PDF file.

## Examples

``` r
# Boxplots without grouping
cyt_bp(ExampleData1[, -c(1:3)], pdf_title = NULL, scale = "log2")

# Boxplots grouped by Group
cyt_bp(ExampleData1[, -c(3,5:28)], group_by = "Group", scale = "zscore")

```
