# Run Random Forest Classification on Cytokine Data,

This function trains and evaluates a Random Forest classification model
on cytokine data. It includes feature importance visualization, cross-
validation for feature selection, and performance metrics such as
accuracy, sensitivity, and specificity. Optionally, for binary
classification, the function also plots the ROC curve and computes the
AUC.

## Usage

``` r
cyt_rf(
  data,
  group_col,
  ntree = 500,
  mtry = 5,
  train_fraction = 0.7,
  plot_roc = FALSE,
  k_folds = 5,
  step = 0.5,
  run_rfcv = TRUE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- data:

  A data frame containing the cytokine data, with one column as the
  grouping variable and the rest as numerical features.

- group_col:

  A string representing the name of the column with the grouping
  variable (the target variable for classification).

- ntree:

  An integer specifying the number of trees to grow in the forest
  (default is 500).

- mtry:

  An integer specifying the number of variables randomly selected at
  each split (default is 5).

- train_fraction:

  A numeric value between 0 and 1 representing the proportion of data to
  use for training (default is 0.7).

- plot_roc:

  A logical value indicating whether to plot the ROC curve and compute
  the AUC for binary classification (default is FALSE).

- k_folds:

  An integer specifying the number of folds for cross-validation
  (default is 5).

- step:

  A numeric value specifying the fraction of variables to remove at each
  step during cross-validation for feature selection (default is 0.5).

- run_rfcv:

  A logical value indicating whether to run Random Forest
  cross-validation for feature selection (default is TRUE).

- verbose:

  A logical value indicating whether to print additional informational
  output to the console. When `TRUE`, the function will display progress
  messages, and intermediate results when `FALSE` (the default), it runs
  quietly.

- seed:

  An integer specifying the seed for reproducibility (default is 123).

## Value

A list containing:

- model:

  The trained Random Forest model.

- confusion_matrix:

  The confusion matrix of the test set predictions.

- importance_plot:

  A ggplot object showing the variable importance plot based on Mean
  Decrease Gini.

- rfcv_result:

  Results from Random Forest cross-validation for feature selection (if
  `run_rfcv` is TRUE).

- importance_data:

  A data frame containing the variable importance based on the Gini
  index.

## Details

The function fits a Random Forest model to the provided data by
splitting it into training and test sets. It calculates performance
metrics such as accuracy, sensitivity, and specificity for both sets.
For binary classification, it can also plot the ROC curve and compute
the AUC. If `run_rfcv` is TRUE, cross-validation is performed to select
the optimal number of features. If `verbose` is TRUE, the function
prints additional information to the console, including training
results, test results, and plots.

## Author

Shubh Saraswat

## Examples

``` r
data.df0 <- ExampleData1
data.df <- data.frame(data.df0[, 1:3], log2(data.df0[, -c(1:3)]))
data.df <- data.df[, -c(2:3)]
data.df <- dplyr::filter(data.df, Group != "ND")

cyt_rf(
  data = data.df, group_col = "Group", k_folds = 5, ntree = 1000,
  mtry = 4, run_rfcv = TRUE, plot_roc = TRUE, verbose = FALSE
)



```
