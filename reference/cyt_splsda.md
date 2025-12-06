# Analyze data with Sparse Partial Least Squares Discriminant Analysis (sPLS-DA).

This function conducts Sparse Partial Least Squares Discriminant
Analysis (sPLS-DA) on the provided data. It uses the specified
`group_col` (and optionally `group_col2`) to define class labels while
assuming the remaining columns contain continuous variables. The
function supports a log2 transformation via the `scale` parameter and
generates a series of plots, including classification plots, scree
plots, loadings plots, and VIP score plots. Optionally, ROC curves are
produced when `roc` is `TRUE`. Additionally, cross-validation is
supported via LOOCV or Mfold methods. When both `group_col` and
`group_col2` are provided and differ, the function analyzes each
treatment level separately.

## Usage

``` r
cyt_splsda(
  data,
  group_col = NULL,
  group_col2 = NULL,
  multilevel_col = NULL,
  batch_col = NULL,
  ind_names = FALSE,
  colors = NULL,
  pdf_title = NULL,
  ellipse = FALSE,
  bg = FALSE,
  conf_mat = FALSE,
  var_num,
  cv_opt = NULL,
  fold_num = 5,
  scale = NULL,
  comp_num = 2,
  pch_values,
  style = NULL,
  roc = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- data:

  A matrix or data frame containing the variables. Columns not specified
  by `group_col` or `group_col2` are assumed to be continuous variables
  for analysis.

- group_col:

  A string specifying the column name that contains the first group
  information. If `group_col2` is not provided, an overall analysis will
  be performed.

- group_col2:

  A string specifying the second grouping column. Default is `NULL`.

- multilevel_col:

  A string specifying the column name that identifies repeated
  measurements (e.g., patient or sample IDs). If provided, a multilevel
  analysis will be performed. Default is `NULL`.

- batch_col:

  A string specifying the column that identifies the batch or study for
  each sample.

- ind_names:

  If `TRUE`, the row names of the first (or second) data matrix is used
  as names. Default is `FALSE`. If a character vector is provided, these
  values will be used as names. If 'pch' is set this will overwrite the
  names as shapes. See ?mixOmics::plotIndiv for details.

- colors:

  A vector of colors for the groups or treatments. If `NULL`, a random
  palette (using `rainbow`) is generated based on the number of groups.

- pdf_title:

  A string specifying the file name for saving the PDF output. Default
  is `NULL` which generates figures in the current graphics device.

- ellipse:

  Logical. Whether to draw a 95\\ figures. Default is `FALSE`.

- bg:

  Logical. Whether to draw the prediction background in the figures.
  Default is `FALSE`.

- conf_mat:

  Logical. Whether to print the confusion matrix for the
  classifications. Default is `FALSE`.

- var_num:

  Numeric. The number of variables to be used in the PLS-DA model.

- cv_opt:

  Character. Option for cross-validation method: either "loocv" or
  "Mfold". Default is `NULL`.

- fold_num:

  Numeric. The number of folds to use if `cv_opt` is "Mfold". Default is
  5.

- scale:

  Character. Option for data transformation; if set to `"log2"`, a log2
  transformation is applied to the continuous variables. Default is
  `NULL`.

- comp_num:

  Numeric. The number of components to calculate in the sPLS-DA model.
  Default is 2.

- pch_values:

  A vector of integers specifying the plotting characters (pch values)
  to be used in the plots.

- style:

  Character. If set to `"3D"` or `"3d"` and `comp_num` equals 3, a 3D
  plot is generated using the `plot3D` package. Default is `NULL`.

- roc:

  Logical. Whether to compute and plot the ROC curve for the model.
  Default is `FALSE`.

- verbose:

  A logical value indicating whether to print additional informational
  output to the console. When `TRUE`, the function will display progress
  messages, and intermediate results when `FALSE` (the default), it runs
  quietly.

- seed:

  An integer specifying the seed for reproducibility (default is 123).

## Value

Plots consisting of the classification figures, component figures with
Variable of Importance in Projection (VIP) scores, and classifications
based on VIP scores greater than 1. ROC curves and confusion matrices
are also produced if requested.

## Details

When `verbose` is set to `TRUE`, additional information about the
analysis and confusion matrices are printed to the console. These can be
suppressed by keeping `verbose = FALSE`.

## References

Lê Cao, K.-A., Boitard, S. and Besse, P. (2011). Sparse PLS Discriminant
Analysis: biologically relevant feature selection and graphical displays
for multiclass problems. *BMC Bioinformatics* **12**:253.

## Author

Xiaohua Douglas Zhang and Shubh Saraswat

## Examples

``` r
# Loading Sample Data
data_df <- ExampleData1[,-c(3)]
data_df <- dplyr::filter(data_df, Group != "ND", Treatment != "Unstimulated")

cyt_splsda(data_df, pdf_title = NULL,
colors = c("black", "purple"), bg = FALSE, scale = "log2",
conf_mat = FALSE, var_num = 25, cv_opt = NULL, comp_num = 2,
pch_values = c(16, 4), style = NULL, ellipse = TRUE,
group_col = "Group", group_col2 = "Treatment", roc = FALSE, verbose = FALSE)
#> Warning: the standard deviation is zero








#> Warning: the standard deviation is zero








```
