#' Run XGBoost Classification on Cytokine Data.
#'
#' This function trains and evaluates an XGBoost classification model on
#' cytokine data.
#' It allows for hyperparameter tuning, cross-validation, and visualizes
#' feature importance.
#'
#' @param data A data frame containing the cytokine data, with one column as
#' the grouping variable and the rest as numerical features.
#' @param group_col A string representing the name of the column with the
#' grouping variable (i.e., the target variable for classification).
#' @param train_fraction A numeric value between 0 and 1 representing the
#' proportion of data to use for training (default is 0.7).
#' @param nrounds An integer specifying the number of boosting rounds
#' (default is 500).
#' @param max_depth An integer specifying the maximum depth of the trees
#' (default is 6).
#' @param learning_rate A numeric value representing the learning rate
#'   (default is 0.1).  This replaces the deprecated `eta` argument.
#' @param eta `r lifecycle::badge("deprecated")` Deprecated; use `learning_rate` instead.
#' @param nfold An integer specifying the number of folds for cross-validation
#' (default is 5).
#' @param cv A logical value indicating whether to perform cross-validation
#' (default is FALSE).
#' @param objective A string specifying the XGBoost objective function
#' (default is "multi:softprob" for multi-class classification).
#' @param early_stopping_rounds An integer specifying the number of rounds
#' with no improvement to stop training early (default is NULL).
#' @param eval_metric A string specifying the evaluation metric
#' (default is "mlogloss").
#' @param min_split_loss A numeric value for the minimum loss reduction
#'   required to make a further partition (default is 0).  This replaces
#'   the deprecated `gamma` argument.
#' @param gamma `r lifecycle::badge("deprecated")` Deprecated; use `min_split_loss` instead.
#' @param colsample_bytree A numeric value specifying the subsample ratio
#' of columns when constructing each tree (default is 1).
#' @param subsample A numeric value specifying the subsample ratio of the
#' training instances (default is 1).
#' @param min_child_weight A numeric value specifying the minimum sum of
#' instance weight needed in a child (default is 1).
#' @param top_n_features An integer specifying the number of top features to
#'  display in the importance plot (default is 10).
#' @param verbose An integer specifying the verbosity of the training
#' process (default is 1).
#' @param plot_roc A logical value indicating whether to plot the ROC curve
#' and calculate the AUC for binary classification (default is \code{FALSE}).
#' @param print_results A logical value indicating whether to print the results
#'  of the model training and evaluation (default is \code{FALSE}). If set to \code{TRUE},
#'  it will print the confusion matrix, and feature importance.
#' @param seed An integer specifying the seed for reproducibility (default is 123).
#'
#' @return A list containing:
#' \item{model}{The trained XGBoost model.}
#' \item{confusion_matrix}{The confusion matrix of the test set predictions.}
#' \item{importance}{The feature importance matrix for the top features.}
#' \item{class_mapping}{A named vector showing the mapping from class labels
#' to numeric values used for training.}
#' \item{cv_results}{Cross-validation results, if cross-validation was
#' performed (otherwise NULL).}
#' \item{plot}{A ggplot object showing the feature importance plot.}
#'
#' @author Shubh Saraswat
#'
#' @details
#' The function allows for training an XGBoost model on cytokine data,
#' splitting the data into training and test sets. If cross-validation is
#' enabled (`cv = TRUE`), it performs k-fold cross-validation and reports the
#' confusion matrix and accuracy.
#' The function also visualizes the top N important
#' features using `xgb.ggplot.importance()`.
#'
#' @examples
#' # Example usage:
#' data_df0 <- ExampleData1
#' data_df <- data.frame(data_df0[, 1:3], log2(data_df0[, -c(1:3)]))
#' data_df <- data_df[, -c(2:3)]
#' data_df <- dplyr::filter(data_df, Group != "ND")
#'
#' cyt_xgb(
#'  data = data_df,
#'  group_col = "Group",
#'  nrounds = 500,
#'  max_depth = 4,
#'  min_split_loss = 0,
#'  learning_rate = 0.05,
#'  nfold = 5,
#'  cv = FALSE,
#'  objective = "multi:softprob",
#'  eval_metric = "auc",
#'  early_stopping_rounds = NULL,
#'  top_n_features = 10,
#'  verbose = 0,
#'  plot_roc = TRUE,
#'  print_results = FALSE)
#'
#' @importFrom xgboost xgb.DMatrix xgb.train xgb.importance xgb.ggplot.importance xgb.cv getinfo
#' @importFrom caret createDataPartition confusionMatrix
#' @import ggplot2
#' @importFrom pROC roc auc ggroc
#' @export
#'
cyt_xgb <- function(
  data,
  group_col,
  train_fraction = 0.7,
  nrounds = 500,
  max_depth = 6,
  eta = lifecycle::deprecated(),
  learning_rate = 0.1,
  nfold = 5,
  cv = FALSE,
  objective = "multi:softprob",
  early_stopping_rounds = NULL,
  eval_metric = "mlogloss",
  gamma = lifecycle::deprecated(),
  min_split_loss = 0,
  colsample_bytree = 1,
  subsample = 1,
  min_child_weight = 1,
  top_n_features = 10,
  verbose = 1,
  plot_roc = FALSE,
  print_results = FALSE,
  seed = 123
) {
  # Handle deprecated eta -> learning_rate
  if (lifecycle::is_present(eta)) {
    lifecycle::deprecate_warn(
      when = "0.2.2",
      what = "cyt_xgb(eta)",
      with = "cyt_xgb(learning_rate)"
    )
    learning_rate <- eta
  }

  # Handle deprecated gamma -> min_split_loss
  if (lifecycle::is_present(gamma)) {
    lifecycle::deprecate_warn(
      when = "0.2.2",
      what = "cyt_xgb(gamma)",
      with = "cyt_xgb(min_split_loss)"
    )
    min_split_loss <- gamma
  }

  # Ensure the grouping variable is a factor
  data[[group_col]] <- as.factor(data[[group_col]])

  # Create a mapping from group names to numeric labels
  class_labels <- levels(data[[group_col]])
  class_mapping <- setNames(0:(length(class_labels) - 1), class_labels)
  if (print_results) {
    cat("\nGroup to Numeric Label Mapping\n")
    print(class_mapping)
  }

  # Convert group column to numeric values using the
  # mapping (starting from 0 for xgboost)
  data[[group_col]] <- as.numeric(data[[group_col]]) - 1

  # Prepare the dataset for xgboost (convert to matrix)
  X <- as.matrix(data[, -which(names(data) == group_col)])
  y <- data[[group_col]] # Numeric class values: 0, 1, 2, ..., n-1

  num_class <- length(unique(y))

  # Split the data into training and testing sets
  set.seed(seed)
  train_indices <- caret::createDataPartition(
    y,
    p = train_fraction,
    list = FALSE
  )
  X_train <- X[train_indices, ]
  y_train <- y[train_indices]
  X_test <- X[-train_indices, ]
  y_test <- y[-train_indices]

  # Prepare the DMatrix objects for xgboost
  dtrain <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
  dtest <- xgboost::xgb.DMatrix(data = X_test, label = y_test)

  # Set parameters for xgboost
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

  # Train the XGBoost model
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
  # Make predictions on the test set
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

  # For binary classification, reshape predictions and compute ROC/AUC
  roc_plot <- NULL

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
  # Feature importance - show only top_n features
  importance <- xgboost::xgb.importance(
    feature_names = colnames(X_train),
    model = xgb_model
  )
  top_features <- head(importance, top_n_features)

  if (print_results) {
    cat("\nTop", top_n_features, "Important Features\n")
    print(top_features)
  }

  ggplot_imp <- xgboost::xgb.ggplot.importance(
    importance_matrix = top_features,
    top_n = top_n_features
  )
  ggplot_imp <- ggplot_imp +
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
  print(ggplot_imp)

  # Cross-Validation (optional)
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
}
