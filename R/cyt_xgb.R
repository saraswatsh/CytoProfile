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
#' @param eta A numeric value representing the learning rate (default is 0.1).
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
#' @param gamma A numeric value for the minimum loss reduction required to
#'  make a further partition (default is 0).
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
#' and calculate the AUC for binary classification (default is FALSE).
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
#' @details
#' The function allows for training an XGBoost model on cytokine data,
#' splitting the data into training and test sets. If cross-validation is
#' enabled (`cv = TRUE`), it performs k-fold cross-validation and prints the
#'  best iteration based on the evaluation metric.
#' The function also visualizes the top N important
#' features using `xgb.ggplot.importance()`.
#'
#' @examples
#' # Example usage:
#' data_df0 <- ExampleData1
#' data_df <- data.frame(data_df0[, 1:3], log2(data_df0[, -c(1:3)]))
#' data_df <- data_df[, -c(2,3)]
#' data_df <- dplyr::filter(data_df, Group != "ND")
#'
#' cyt_xgb(
#'   data = data_df, group_col = "Group",
#'   nrounds = 500, max_depth = 4, eta = 0.05,
#'   nfold = 5, cv = TRUE, eval_metric = "mlogloss",
#'   early_stopping_rounds = NULL, top_n_features = 10,
#'   verbose = 0, plot_roc = TRUE
#' )
#'
#' @importFrom xgboost xgb.DMatrix xgb.train xgb.importance xgb.ggplot.importance xgb.cv getinfo
#' @importFrom caret createDataPartition confusionMatrix
#' @import ggplot2
#' @importFrom pROC roc auc ggroc
#' @export
#'
cyt_xgb <- function(data, group_col, train_fraction = 0.7,
                    nrounds = 500, max_depth = 6, eta = 0.1,
                    nfold = 5, cv = FALSE,
                    objective = "multi:softprob", early_stopping_rounds = NULL,
                    eval_metric = "mlogloss",
                    gamma = 0, colsample_bytree = 1, subsample = 1,
                    min_child_weight = 1,
                    top_n_features = 10, verbose = 1, plot_roc = FALSE) {
  # Ensure the grouping variable is a factor
  data[[group_col]] <- as.factor(data[[group_col]])

  # Create a mapping from group names to numeric labels
  class_labels <- levels(data[[group_col]])
  class_mapping <- setNames(0:(length(class_labels) - 1), class_labels)
  cat("\n### Group to Numeric Label Mapping ###\n")
  print(class_mapping)

  # Convert group column to numeric values using the
  # mapping (starting from 0 for xgboost)
  data[[group_col]] <- as.numeric(data[[group_col]]) - 1

  # Prepare the dataset for xgboost (convert to matrix)
  X <- as.matrix(data[, -which(names(data) == group_col)])
  y <- data[[group_col]] # Numeric class values: 0, 1, 2, ..., n-1

  # Split the data into training and testing sets
  set.seed(123)
  train_indices <- caret::createDataPartition(y, p = train_fraction, list = FALSE)
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
    eta = eta,
    gamma = gamma,
    colsample_bytree = colsample_bytree,
    subsample = subsample,
    min_child_weight = min_child_weight
  )

  # Train the XGBoost model
  cat("\n### TRAINING XGBOOST MODEL ###\n")
  xgb_model <- xgboost::xgb.train(
    params = params, data = dtrain, nrounds = nrounds,
    watchlist = list(train = dtrain, test = dtest),
    early_stopping_rounds = early_stopping_rounds, verbose = verbose
  )

  # Determine the evaluation column name dynamically (e.g., "test_mlogloss")
  eval_col_name <- paste0("test_", eval_metric)
  cat("\nBest iteration from training (based on", eval_metric, "):\n")
  print(xgb_model$evaluation_log[which.min(
    xgb_model$evaluation_log[[eval_col_name]]
  ), ])

  # Make predictions on the test set
  preds <- predict(xgb_model, X_test)
  pred_labels <- max.col(matrix(preds,
    ncol = length(unique(y_test)),
    byrow = TRUE
  )) - 1

  # For binary classification, reshape predictions and compute ROC/AUC
  if (length(unique(y_test)) == 2) {
    preds_matrix <- matrix(preds, ncol = 2, byrow = TRUE)
    xgb_prob <- preds_matrix[, 2]
    if (length(xgb_prob) != length(y_test)) {
      cat("The length of predicted probabilities does not match the length of
          true labels.")
    }
    roc_obj <- pROC::roc(y_test, xgb_prob)
    auc_value <- pROC::auc(roc_obj)
    cat("\nAUC: ", auc_value, "\n")

    roc_plot <- pROC::ggroc(roc_obj,
      color = "blue", linewidth = 1.5,
      legacy.axes = TRUE
    ) +
      geom_abline(linetype = "dashed", color = "red", linewidth = 1) +
      labs(
        title = "ROC Curve (Test Set)",
        x = "1 - Specificity",
        y = "Sensitivity"
      ) +
      annotate("text",
        x = 0.75, y = 0.25,
        label = paste("AUC =", round(auc_value, 3)),
        size = 5, color = "blue"
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95")
      )
    print(roc_plot)
  } else {
    cat("ROC curve is only available for binary classification.")
  }

  # Confusion matrix on test set
  cat("\n### Confusion Matrix on Test Set ###\n")
  confusion_mat <- caret::confusionMatrix(as.factor(pred_labels), as.factor(y_test))
  print(confusion_mat)

  # Feature importance - show only top_n features
  importance <- xgboost::xgb.importance(
    feature_names = colnames(X_train),
    model = xgb_model
  )
  top_features <- head(importance, top_n_features)
  cat("\n### Top", top_n_features, "Important Features ###\n")
  print(top_features)

  ggplot_imp <- xgboost::xgb.ggplot.importance(
    importance_matrix = top_features,
    top_n = top_n_features
  )
  ggplot_imp <- ggplot_imp +
    geom_bar(stat = "identity", fill = "red2", show.legend = FALSE) +
    ggtitle("Top Features by Gain") +
    ylab("Importance (Gain)") +
    xlab("Features") +
    theme_minimal() +
    theme(
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
    cat("\n### CROSS-VALIDATION USING XGBOOST ###\n")
    xgb_cv <- xgboost::xgb.cv(
      params = params, data = dtrain, nrounds = nrounds, nfold = nfold,
      early_stopping_rounds = early_stopping_rounds, verbose = verbose,
      prediction = TRUE
    )

    cat("\nBest iteration from cross-validation:\n")
    eval_col_name_cv <- paste0("test_", eval_metric, "_mean")
    print(xgb_cv$evaluation_log[which.min(
      xgb_cv$evaluation_log[[eval_col_name_cv]]
    ), ])

    cv_preds <- xgb_cv$pred
    num_class <- length(unique(y))
    if (num_class == 2) {
      cv_pred_labels <- ifelse(cv_preds[, 2] > cv_preds[, 1], 1, 0)
    } else {
      cv_pred_labels <- max.col(cv_preds) - 1
    }

    actual_labels <- xgboost::getinfo(dtrain, "label")
    cv_confusion_mat <- caret::confusionMatrix(
      as.factor(cv_pred_labels),
      as.factor(actual_labels)
    )
    print(cv_confusion_mat)

    cv_accuracy <- sum(cv_pred_labels == actual_labels) / length(actual_labels)
    cat("\nCross-Validation Accuracy: ", cv_accuracy, "\n")
  }

  return(list(
    model = xgb_model,
    confusion_matrix = confusion_mat,
    importance = top_features,
    class_mapping = class_mapping,
    cv_results = if (cv) xgb_cv else NULL,
    plot = ggplot_imp
  ))
}
