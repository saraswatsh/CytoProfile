# Violin Plots for Continuous Variables with Optional Grouping

`cyt_violin` produces violin plots for each numeric variable in `data`,
optionally grouped by one or more categorical variables. When grouping
is not specified, the function behaves similarly to `cyt_bp` but uses
violins instead of boxplots and supports pagination via the `bin_size`
argument. When grouping is provided, a separate violin is drawn for each
level (or interaction of levels) of the grouping variables. Users may
optionally overlay boxplots within each violin to visualize the median
and interquartile range.

## Usage

``` r
cyt_violin(
  data,
  pdf_title = NULL,
  group_by = NULL,
  bin_size = 25,
  y_lim = NULL,
  scale = c("none", "log2", "log10", "zscore", "custom"),
  custom_fn = NULL,
  boxplot_overlay = FALSE
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

  Integer. Maximum number of violins per page when grouping is not used.
  Default is 25, mirroring `cyt_bp`.

- y_lim:

  Optional numeric vector giving y‑axis limits for the plots. Applies to
  all plots.

- scale:

  Character specifying a transformation for numeric variables. Accepts
  "none", "log2", "log10", "zscore", or "custom". When "custom", supply
  a function via `custom_fn`.

- custom_fn:

  A user supplied function to transform numeric columns when
  `scale = "custom"`.

- boxplot_overlay:

  Logical. When `TRUE`, a narrow boxplot is drawn inside each violin to
  summarize the median and quartiles. Default is `FALSE`.

## Value

Invisibly returns a list of `ggplot` objects. When `pdf_title` is
provided, plots are written to the PDF file.

## Author

Shubh Saraswat

## Examples

``` r
# Violin plots without grouping
cyt_violin(ExampleData1[, -c(1:3)], pdf_title = NULL, scale = "zscore")

# Violin plots grouped by Group with boxplot overlay
cyt_violin(ExampleData1[, -c(3,5:28)], group_by = "Group",
                        boxplot_overlay = TRUE, scale = "log2")
```
