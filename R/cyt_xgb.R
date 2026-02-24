#' Run XGBoost Classification on Cytokine Data.
#'
#' This function trains and evaluates an XGBoost classification model on
#' cytokine data.
#' It allows for hyperparameter tuning, cross-validation, and visualizes
#' feature importance.
#'
#'
#' @param data A data frame containing numeric predictor variables and
#'   one grouping column.  Non-numeric predictor columns will be
#'   coerced to numeric if possible.
#' @param group_col Character string naming the column that contains
#'   the class labels (target variable).  Required.
#' @param train_fraction Numeric between 0 and 1 specifying the
#'   proportion of samples used for model training.  The remainder is
#'   reserved for testing.  Default is 0.7.
#' @param nrounds Integer specifying the number of boosting rounds.
#'   Default is 500.
#' @param max_depth Integer specifying the maximum depth of each tree.
#'   Default is 6.
#' @param learning_rate Numeric specifying the learning rate.  Default
#'   is 0.1.
#' @param nfold Integer specifying the number of folds for
#'   cross-validation when `cv = TRUE`.  Default is 5.
#' @param cv Logical indicating whether to perform cross-validation
#'   using `xgb.cv`.  Default is `FALSE`.
#' @param objective Character string specifying the objective
#'   function.  For multi-class classification use "multi:softprob";
#'   for binary classification use "binary:logistic".  Default is
#'   "multi:softprob".
#' @param early_stopping_rounds Integer specifying the number of rounds
#'   without improvement to trigger early stopping.  Default is NULL
#'   (no early stopping).
#' @param eval_metric Character specifying the evaluation metric used
#'   during training.  Default is "mlogloss".
#' @param min_split_loss Numeric specifying the minimum loss reduction
#'   required to make a further partition.  Default is 0.
#' @param colsample_bytree Numeric specifying the subsample ratio of
#'   columns when constructing each tree.  Default is 1.
#' @param subsample Numeric specifying the subsample ratio of the
#'   training instances.  Default is 1.
#' @param min_child_weight Numeric specifying the minimum sum of
#'   instance weight needed in a child.  Default is 1.
#' @param top_n_features Integer specifying the number of top features
#'   to display in the importance plot.  Default is 10.
#' @param verbose Integer (0, 1 or 2) controlling the verbosity of
#'   `xgb.train`.  Default is 1.  Larger values print more information.
#' @param plot_roc Logical indicating whether to plot the ROC curve and
#'   compute the AUC for binary classification.  Default is `FALSE`.
#' @param print_results Logical.  If `TRUE`, prints the confusion
#'   matrix, top features and cross-validation metrics to the
#'   console.  Default is `FALSE`.
#' @param seed Optional integer seed for reproducibility.  Default is
#'   123.
#' @param scale Character string specifying a transformation to apply
#'   to numeric predictors prior to model fitting.  Possible values
#'   are "none", "log2", "log10", "zscore" or "custom".  When set
#'   to "custom" a user defined function must be supplied via
#'   `custom_fn`.  Defaults to "none".
#' @param custom_fn A custom transformation function used when
#'   `scale = "custom"`.  Ignored otherwise.  Should take a numeric
#'   vector and return a numeric vector of the same length.
#'
#' @return An invisible list with elements: `model` (the trained
#'   xgboost object), `confusion_matrix` (test set confusion
#'   matrix), `importance` (variable importance matrix),
#'   `class_mapping` (mapping from class names to numeric labels),
#'   `cv_results` (cross-validation results when `cv=TRUE`),
#'   `importance_plot` (a `ggplot2` object of the top feature
#'   importance), and `roc_plot` (a ROC curve for binary
#'   classification when `plot_roc=TRUE`).  All plots are printed
#'   automatically.
#' @examples
#' data_df0 <- ExampleData1
#' data_df <- data.frame(data_df0[, 1:3], log2(data_df0[, -c(1:3)]))
#' data_df <- data_df[, -c(2:3)]
#' data_df <- dplyr::filter(data_df, Group != "ND")
#' cyt_xgb(
#'  data = data_df,
#'  group_col = "Group",
#'  nrounds = 250,
#'  max_depth = 4,
#'  min_split_loss = 0,
#'  learning_rate = 0.05,
#'  nfold = 5,
#'  cv = FALSE,
#'  objective = "multi:softprob",
#'  eval_metric = "auc",
#'  plot_roc = TRUE,
#'  print_results = FALSE)
#'
#' @importFrom xgboost xgb.DMatrix xgb.train xgb.importance xgb.ggplot.importance xgb.cv getinfo
#' @importFrom caret createDataPartition confusionMatrix
#' @import ggplot2
#' @importFrom pROC roc auc ggroc
#' @importFrom data.table copy
#' @export
cyt_xgb <- function(
  data,
  group_col,
  train_fraction = 0.7,
  nrounds = 500,
  max_depth = 6,
  learning_rate = 0.1,
  nfold = 5,
  cv = FALSE,
  objective = "multi:softprob",
  early_stopping_rounds = NULL,
  eval_metric = "mlogloss",
  min_split_loss = 0,
  colsample_bytree = 1,
  subsample = 1,
  min_child_weight = 1,
  top_n_features = 10,
  verbose = 0,
  plot_roc = FALSE,
  print_results = FALSE,
  seed = 123,
  scale = c("none", "log2", "log10", "zscore", "custom"),
  custom_fn = NULL
) {
  names(data) <- make.names(names(data), unique = TRUE)
  # Convert to data frame
  data <- as.data.frame(data)
  # Check group column
  if (!group_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data.", group_col))
  }
  # Apply optional scaling to numeric predictors (excluding the group column)
  scale <- match.arg(scale)
  id_cols <- group_col
  numeric_cols <- setdiff(names(data)[sapply(data, is.numeric)], id_cols)
  if (length(numeric_cols) > 0) {
    data <- apply_scale(
      data,
      columns = numeric_cols,
      scale = scale,
      custom_fn = custom_fn
    )
  }
  # Ensure grouping variable is a factor
  data[[group_col]] <- as.factor(data[[group_col]])
  # Create mapping from class labels to numeric values (0..n-1)
  class_labels <- levels(data[[group_col]])
  class_mapping <- setNames(seq_along(class_labels) - 1, class_labels)
  # Print mapping if requested
  if (print_results) {
    cat("\nGroup to Numeric Label Mapping\n")
    print(class_mapping)
  }
  # Replace group column with numeric labels
  data[[group_col]] <- as.numeric(data[[group_col]]) - 1
  # Prepare predictor matrix and label vector
  X <- as.matrix(data[, setdiff(names(data), group_col), drop = FALSE])
  y <- data[[group_col]]
  num_class <- length(unique(y))
  # Split into training and testing sets
  set.seed(seed)
  train_idx <- caret::createDataPartition(y, p = train_fraction, list = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  # Create DMatrix objects
  dtrain <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
  dtest <- xgboost::xgb.DMatrix(data = X_test, label = y_test)
  # Assemble parameters list
  params <- list(
    objective = objective,
    eval_metric = eval_metric,
    num_class = length(unique(y)),
    max_depth = max_depth,
    learning_rate = learning_rate,
    min_split_loss = min_split_loss,
    colsample_bytree = colsample_bytree,
    subsample = subsample,
    min_child_weight = min_child_weight
  )
  # Train XGBoost model with optional early stopping
  if (print_results) {
    cat("\nTRAINING XGBOOST MODEL\n")
  }
  xgb_model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    evals = list(train = dtrain, test = dtest),
    early_stopping_rounds = early_stopping_rounds,
    verbose = verbose
  )
  # Predictions on test set
  preds <- predict(xgb_model, X_test)

  # XGBoost >= 1.4: multi:softprob returns a matrix [nrow, num_class]
  # Older versions might return a flat vector, so handle both.
  if (is.matrix(preds)) {
    preds_matrix <- preds
  } else {
    preds_matrix <- matrix(preds, ncol = num_class, byrow = TRUE)
  }

  # Hard labels for confusion matrix (0 .. num_class-1)
  pred_labels <- max.col(preds_matrix) - 1

  # Confusion matrix on test set
  if (print_results) {
    cat("\nConfusion Matrix on Test Set\n")
  }
  # Ensure both factors share the same levels (0 .. num_class-1)
  test_levels <- sort(unique(y)) # global class set

  test_pred <- factor(pred_labels, levels = test_levels)
  test_ref <- factor(y_test, levels = test_levels)

  confusion_mat <- caret::confusionMatrix(test_pred, test_ref)

  if (print_results) {
    print(confusion_mat)
  }

  # Compute feature importance and create importance plot
  importance <- xgboost::xgb.importance(
    feature_names = colnames(X_train),
    model = xgb_model
  )
  top_features <- head(importance, top_n_features)
  if (print_results) {
    cat("\nTop", top_n_features, "Important Features\n")
    print(top_features)
  }
  if (!requireNamespace("Ckmeans.1d.dp", quietly = TRUE)) {
    warning(
      "Install 'Ckmeans.1d.dp' for a clustered importance plot. Falling back to a basic bar chart."
    )
    imp_plot <- ggplot2::ggplot(
      top_features,
      ggplot2::aes(x = reorder(Feature, Gain), y = Gain)
    ) +
      ggplot2::geom_bar(stat = "identity", fill = "red2") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Top Features by Gain",
        x = "Features",
        y = "Importance (Gain)"
      ) +
      ggplot2::theme_minimal()
  } else {
    imp_plot <- xgboost::xgb.ggplot.importance(
      importance_matrix = data.table::copy(top_features),
      top_n = top_n_features
    ) +
      ggplot2::geom_bar(stat = "identity", fill = "red2", show.legend = FALSE) +
      ggplot2::ggtitle("Top Features by Gain") +
      ggplot2::ylab("Importance (Gain)") +
      ggplot2::xlab("Features") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        legend.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(color = "black", size = 12, face = "bold"),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black")
      )
  }
  # Print importance plot automatically
  print(imp_plot)
  # Cross-validation using xgb.cv
  cv_results <- NULL
  if (cv) {
    if (print_results) {
      cat("\nCROSS-VALIDATION USING XGBOOST\n")
    }
    xgb_cv <- xgboost::xgb.cv(
      params = params,
      data = dtrain,
      nrounds = nrounds,
      nfold = nfold,
      early_stopping_rounds = early_stopping_rounds,
      verbose = verbose,
      prediction = TRUE
    )
    # Extract predictions and compute CV confusion matrix and accuracy
    cv_preds <- xgb_cv$cv_predict
    num_class <- length(unique(y))
    if (num_class == 2) {
      cv_pred_labels <- ifelse(cv_preds$pred[, 2] > cv_preds$pred[, 1], 1, 0)
    } else {
      cv_pred_labels <- max.col(cv_preds$pred) - 1
    }

    actual_labels <- xgboost::getinfo(dtrain, "label")

    # Use a common set of factor levels for CV predictions and labels
    class_levels <- 0:(num_class - 1)

    cv_pred <- factor(cv_pred_labels, levels = class_levels)
    cv_ref <- factor(actual_labels, levels = class_levels)

    if (length(cv_pred) != length(cv_ref)) {
      warning(
        "Length mismatch between CV predictions and labels; ",
        "skipping CV confusion matrix."
      )
      cv_confusion_mat <- NULL
    } else {
      cv_confusion_mat <- caret::confusionMatrix(cv_pred, cv_ref)
    }

    if (print_results && !is.null(cv_confusion_mat)) {
      cat("\nCross-Validation Confusion Matrix\n")
      print(cv_confusion_mat)
    }

    cv_accuracy <- sum(cv_pred_labels == actual_labels) / length(actual_labels)
    if (print_results) {
      cat("\nCross-Validation Accuracy: ", round(cv_accuracy, 3), "\n")
    }
  }
  # ROC curve for binary classification
  if (plot_roc) {
    if (num_class == 2) {
      # probability of class "1" (label 1)
      xgb_prob <- preds_matrix[, 2]

      if (length(xgb_prob) != length(y_test)) {
        warning(
          "Length mismatch between predicted probabilities and true labels; ",
          "skipping ROC/AUC."
        )
      } else {
        roc_obj <- pROC::roc(y_test, xgb_prob, quiet = TRUE)
        auc_value <- pROC::auc(roc_obj)

        if (print_results) {
          cat("\nAUC: ", round(auc_value, 3), "\n")
        }

        roc_plot <- pROC::ggroc(
          roc_obj,
          color = "blue",
          linewidth = 1.5,
          legacy.axes = TRUE
        ) +
          ggplot2::geom_abline(
            linetype = "dashed",
            color = "red",
            linewidth = 1
          ) +
          ggplot2::labs(
            title = "ROC Curve (Test Set)",
            x = "1 - Specificity",
            y = "Sensitivity"
          ) +
          ggplot2::annotate(
            "text",
            x = 0.75,
            y = 0.25,
            label = paste("AUC =", round(auc_value, 3)),
            size = 5,
            color = "blue"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_line(color = "grey95")
          )

        print(roc_plot)
      }
    } else {
      warning("ROC curve is only available for binary classification.")
    }
  }
  # Build result list
  res <- list(
    model = xgb_model,
    confusion_matrix = confusion_mat,
    importance = importance,
    class_mapping = class_mapping,
    cv_results = cv_results,
    importance_plot = imp_plot,
    roc_plot = roc_plot
  )
  return(invisible(res))
}
