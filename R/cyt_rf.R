#' Run Random Forest Classification on Cytokine Data,
#'
#' This function trains and evaluates a Random Forest classification model on
#' cytokine data. It includes feature importance visualization, cross-
#' validation for feature selection, and performance metrics such as accuracy,
#' sensitivity, and specificity. Optionally, for binary classification, the
#' function also plots the ROC curve and computes the AUC.
#'
#'
#' @param data A data frame containing the cytokine measurements.  One column
#'   should correspond to the grouping variable (the outcome) and the
#'   remaining columns should be numeric predictors.
#' @param group_col A string naming the column in `data` that contains the
#'   grouping variable.
#' @param ntree Integer specifying the number of trees to grow.  Default is
#'   500.
#' @param mtry Integer specifying the number of variables randomly sampled at
#'   each split.  Default is 5.
#' @param train_fraction Numeric between 0 and 1 giving the proportion of data
#'   used for training.  The remainder is used for testing.  Default is 0.7.
#' @param plot_roc Logical.  If `TRUE` and the problem is binary, an ROC curve
#'   and AUC will be computed and plotted for the test set.  Default is
#'   `FALSE`.
#' @param k_folds Integer specifying the number of folds for `rfcv` when
#'   `run_rfcv = TRUE`.  Default is 5.
#' @param step Numeric specifying the fraction of variables removed at each
#'   step during `rfcv`.  Default is 0.5.
#' @param run_rfcv Logical indicating whether to run Random Forest
#'   cross-validation for feature selection.  Default is `TRUE`.
#' @param verbose Logical indicating whether to print intermediate results.
#'   When `TRUE`, training and test performance metrics, confusion matrices
#'   and cross-validation details are printed.  Default is `FALSE`.
#' @param seed Optional integer seed for reproducibility.  Default is 123.
#' @param cv Logical indicating whether to perform a separate k-fold
#'   classification cross-validation using `caret`.  Default is `FALSE`.
#' @param cv_folds Integer specifying the number of folds for classification
#'   cross-validation when `cv = TRUE`.  Default is 5.
#' @param scale Character string specifying a transformation to apply to the
#'   numeric predictor columns prior to model fitting.  Options are
#'   "none" (no transformation), "log2", "log10", "zscore", or
#'   "custom".  When "custom" is selected a user defined function
#'   must be supplied via `custom_fn`.  Defaults to "none".
#' @param custom_fn A custom transformation function used when `scale = "custom"`.
#'   The function should take a numeric vector and return a numeric vector
#'   of the same length.  Ignored for other values of `scale`.
#'
#' @return An invisible list with components:
#'   \item{model}{The fitted `randomForest` model.}
#'   \item{confusion_matrix}{Confusion matrix on the test set.}
#'   \item{importance_plot}{A `ggplot2` object of the variable importance
#'     (mean decrease in Gini).}
#'   \item{importance_data}{A data frame of variable importance values.}
#'   \item{rfcv_result}{The `rfcv` object returned when `run_rfcv = TRUE`.}
#'   \item{rfcv_plot}{A `ggplot2` object of cross-validation error versus
#'     number of variables, returned when `run_rfcv = TRUE`.}
#'   \item{rfcv_data}{A data frame summarizing the `rfcv` error curve.}
#'   \item{roc_plot}{A `ggplot2` object of the ROC curve for binary
#'     classification when `plot_roc = TRUE`.}
#'   \item{cv_results}{A `caret` train object returned when `cv = TRUE` or
#'     `NULL` otherwise.}
#'
#' @details
#' The function first coerces the grouping variable to a factor and splits
#' the dataset into training and test subsets according to
#' `train_fraction`.  A Random Forest classifier is fit to the training
#' data using the specified `ntree` and `mtry` parameters.  The model
#' performance is assessed on both the training and test sets, and
#' results are printed when `verbose = TRUE`.  If `plot_roc = TRUE` and
#' the grouping variable has exactly two levels, an ROC curve is computed
#' on the test set and a plot is returned.  Variable importance is
#' extracted and visualized with a bar plot.  Optionally, cross-
#' validation for feature selection (`rfcv`) is performed and the error
#' curve is plotted.  A separate k-fold classification cross-
#' validation using `caret::train` can be requested via `cv = TRUE`.
#'
#' @examples
#' data.df0 <- ExampleData1
#' data.df <- data.frame(data.df0[, 1:3], log2(data.df0[, -c(1:3)]))
#' data.df <- data.df[, -c(2:3)]
#' data.df <- dplyr::filter(data.df, Group != "ND")
#'
#' cyt_rf(
#'   data = data.df, group_col = "Group", k_folds = 5, ntree = 1000,
#'   mtry = 4, run_rfcv = TRUE, plot_roc = TRUE, verbose = FALSE
#' )
#' @importFrom randomForest randomForest rfcv importance
#' @importFrom caret createDataPartition confusionMatrix
#' @import ggplot2
#' @importFrom pROC roc auc ggroc
#' @importFrom utils capture.output
#' @export
cyt_rf <- function(
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
  seed = 123,
  cv = FALSE,
  cv_folds = 5,
  scale = c("none", "log2", "log10", "zscore", "custom"),
  custom_fn = NULL
) {
  names(data) <- make.names(names(data), unique = TRUE)
  # Ensure grouping variable exists and is a factor
  if (!group_col %in% colnames(data)) {
    stop(sprintf("Column '%s' not found in data.", group_col))
  }
  # Apply optional scaling transformation to numeric predictors before modelling
  # Exclude the grouping column from the transformation even if it is numeric
  scale <- match.arg(scale)
  numeric_cols <- setdiff(names(data)[sapply(data, is.numeric)], group_col)
  if (length(numeric_cols) > 0) {
    data <- apply_scale(
      data,
      columns = numeric_cols,
      scale = scale,
      custom_fn = custom_fn
    )
  }
  data[[group_col]] <- as.factor(data[[group_col]])
  if (length(levels(data[[group_col]])) < 2) {
    stop("Grouping variable must have at least two levels.")
  }
  # Split into training and testing sets
  set.seed(seed)
  train_idx <- caret::createDataPartition(
    data[[group_col]],
    p = train_fraction,
    list = FALSE
  )
  train_data <- data[train_idx, , drop = FALSE]
  test_data <- data[-train_idx, , drop = FALSE]
  # Construct formula for random forest
  predictors <- setdiff(colnames(data), group_col)
  rf_formula <- as.formula(paste(
    group_col,
    "~",
    paste(predictors, collapse = "+")
  ))
  # Fit model
  rf_model <- randomForest::randomForest(
    rf_formula,
    data = train_data,
    ntree = ntree,
    mtry = mtry,
    importance = TRUE,
    do.trace = FALSE
  )
  # Training performance
  if (verbose) {
    cat("\n### RANDOM FOREST RESULTS ON TRAINING SET ###\n")
    cat(paste(utils::capture.output(rf_model), collapse = "\n"), "\n")
  }
  train_conf <- rf_model$confusion[, 1:(ncol(rf_model$confusion) - 1)]
  train_acc <- sum(diag(train_conf)) / sum(train_conf)
  if (verbose) {
    cat("\nAccuracy on training set:", round(train_acc, 3), "\n")
    # Sensitivity and specificity by class
    for (i in seq_len(nrow(train_conf))) {
      tp <- train_conf[i, i]
      fn <- sum(train_conf[i, ]) - tp
      fp <- sum(train_conf[, i]) - tp
      tn <- sum(train_conf) - (tp + fn + fp)
      sens <- tp / (tp + fn)
      spec <- tn / (tn + fp)
      cat(sprintf(
        "\nTraining set - Class '%s' metrics:\n",
        rownames(train_conf)[i]
      ))
      cat(sprintf("  Sensitivity: %.3f\n", sens))
      cat(sprintf("  Specificity: %.3f\n", spec))
    }
  }
  # Predictions on test set
  rf_pred <- predict(rf_model, newdata = test_data)
  test_conf <- caret::confusionMatrix(rf_pred, test_data[[group_col]])
  if (verbose) {
    cat("\n### PREDICTIONS ON TEST SET ###\n")
    cat(
      paste(utils::capture.output(print(test_conf$table)), collapse = "\n"),
      "\n"
    )
    cat(
      "\nAccuracy on test set:",
      round(test_conf$overall["Accuracy"], 3),
      "\n"
    )
    if (length(levels(data[[group_col]])) == 2) {
      sens_val <- test_conf$byClass["Sensitivity"]
      spec_val <- test_conf$byClass["Specificity"]
      cat("\nSensitivity by class:\n")
      cat(sprintf("Class: %s: %.3f\n", levels(data[[group_col]])[1], sens_val))
      cat(sprintf(
        "Class: %s: %.3f\n",
        levels(data[[group_col]])[2],
        1 - spec_val
      ))
      cat("\nSpecificity by class:\n")
      cat(sprintf("Class: %s: %.3f\n", levels(data[[group_col]])[2], spec_val))
      cat(sprintf(
        "Class: %s: %.3f\n",
        levels(data[[group_col]])[1],
        1 - sens_val
      ))
    } else {
      sens_vec <- test_conf$byClass[, "Sensitivity"]
      spec_vec <- test_conf$byClass[, "Specificity"]
      cat("\nSensitivity by class:\n")
      for (i in seq_along(sens_vec)) {
        cat(sprintf(
          "Class: %s: %.3f\n",
          levels(data[[group_col]])[i],
          sens_vec[i]
        ))
      }
      cat("\nSpecificity by class:\n")
      for (i in seq_along(spec_vec)) {
        cat(sprintf(
          "Class: %s: %.3f\n",
          levels(data[[group_col]])[i],
          spec_vec[i]
        ))
      }
    }
  }
  # ROC for binary classification
  roc_plot <- NULL
  if (plot_roc && length(levels(data[[group_col]])) == 2) {
    rf_prob <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    roc_obj <- pROC::roc(test_data[[group_col]], rf_prob, quiet = TRUE)
    auc_val <- pROC::auc(roc_obj)
    if (verbose) {
      cat("\nAUC:", round(auc_val, 3), "\n")
    }
    roc_plot <- pROC::ggroc(
      roc_obj,
      color = "blue",
      linewidth = 1.5,
      legacy.axes = TRUE
    ) +
      ggplot2::geom_abline(linetype = "dashed", color = "red", linewidth = 1) +
      ggplot2::labs(
        title = "ROC Curve (Test Set)",
        x = "1 - Specificity",
        y = "Sensitivity"
      ) +
      ggplot2::annotate(
        "text",
        x = 0.75,
        y = 0.25,
        label = paste("AUC =", round(auc_val, 3)),
        size = 5,
        color = "blue"
      ) +
      ggplot2::theme_minimal()
    # Do not print here; include in return list
  }
  # Variable importance data and plot
  importance_data <- data.frame(
    Variable = rownames(randomForest::importance(rf_model)),
    Gini = randomForest::importance(rf_model)[, "MeanDecreaseGini"]
  )
  importance_plot <- ggplot2::ggplot(
    importance_data,
    ggplot2::aes(x = reorder(Variable, Gini), y = Gini)
  ) +
    ggplot2::geom_bar(stat = "identity", fill = "red2") +
    ggplot2::coord_flip() +
    ggplot2::ggtitle("Variable Importance Plot (Mean Decrease in Gini)") +
    ggplot2::xlab("Features") +
    ggplot2::ylab("Importance (Gini Index)") +
    ggplot2::theme_minimal()
  # Feature selection cross-validation using rfcv
  rfcv_result <- NULL
  rfcv_plot <- NULL
  rfcv_data <- NULL
  if (run_rfcv) {
    if (verbose) {
      cat("\n### RANDOM FOREST CROSS-VALIDATION FOR FEATURE SELECTION ###\n")
    }
    x_train <- train_data[, predictors, drop = FALSE]
    y_train <- train_data[[group_col]]
    rfcv_result <- randomForest::rfcv(
      x_train,
      y_train,
      cv.fold = k_folds,
      step = step,
      do.trace = FALSE
    )
    rfcv_data <- data.frame(
      Variables = rfcv_result$n.var,
      Error = rfcv_result$error.cv
    )
    rfcv_plot <- ggplot2::ggplot(
      rfcv_data,
      ggplot2::aes(x = Variables, y = Error)
    ) +
      ggplot2::geom_line(color = "blue") +
      ggplot2::geom_point(color = "blue") +
      ggplot2::ggtitle("Cross-Validation Error vs. Number of Variables") +
      ggplot2::xlab("Number of Variables") +
      ggplot2::ylab("Cross-Validation Error") +
      ggplot2::theme_minimal()
    if (verbose) {
      cat("Random Forest CV completed for feature selection.\n")
    }
  }
  # Optional classification cross-validation via caret
  cv_results <- NULL
  if (cv) {
    if (verbose) {
      cat("\n### RANDOM FOREST k-FOLD CROSS-VALIDATION ###\n")
    }
    # caret requires predictors and outcome separately
    tr_ctrl <- caret::trainControl(
      method = "cv",
      number = cv_folds,
      classProbs = TRUE
    )
    cv_results <- caret::train(
      x = data[, predictors, drop = FALSE],
      y = data[[group_col]],
      method = "rf",
      trControl = tr_ctrl,
      tuneGrid = data.frame(mtry = mtry),
      ntree = ntree
    )
    if (verbose) {
      cat(
        "Cross-validation Accuracy:",
        round(max(cv_results$results$Accuracy), 3),
        "\n"
      )
    }
  }
  # Build and return results list invisibly
  results <- list(
    model = rf_model,
    confusion_matrix = test_conf,
    importance_plot = importance_plot,
    importance_data = importance_data,
    rfcv_result = rfcv_result,
    rfcv_plot = rfcv_plot,
    rfcv_data = rfcv_data,
    roc_plot = roc_plot,
    cv_results = cv_results
  )
  # Automatically display plots when available so the user does not need to assign the result
  # Importance plot
  if (!is.null(importance_plot)) {
    print(importance_plot)
  }
  # Feature selection cross-validation plot
  if (run_rfcv && !is.null(rfcv_plot)) {
    print(rfcv_plot)
  }
  # ROC curve for binary classification
  if (plot_roc && !is.null(roc_plot)) {
    print(roc_plot)
  }
  return(invisible(results))
}
