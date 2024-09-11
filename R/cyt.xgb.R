#' Run XGBoost Classification on Cytokine Data
#'
#' This function trains and evaluates an XGBoost classification model on cytokine data.
#' It allows for hyperparameter tuning, cross-validation, and visualizes feature importance.
#'
#' @param data A data frame containing the cytokine data, with one column as the grouping variable and the rest as numerical features.
#' @param group_col A string representing the name of the column with the grouping variable (i.e., the target variable for classification).
#' @param train_fraction A numeric value between 0 and 1 representing the proportion of data to use for training (default is 0.7).
#' @param nrounds An integer specifying the number of boosting rounds (default is 500).
#' @param max_depth An integer specifying the maximum depth of the trees (default is 6).
#' @param eta A numeric value representing the learning rate (default is 0.1).
#' @param nfold An integer specifying the number of folds for cross-validation (default is 5).
#' @param cv A logical value indicating whether to perform cross-validation (default is FALSE).
#' @param objective A string specifying the XGBoost objective function (default is "multi:softprob" for multi-class classification).
#' @param early_stopping_rounds An integer specifying the number of rounds with no improvement to stop training early (default is NULL).
#' @param eval_metric A string specifying the evaluation metric (default is "mlogloss").
#' @param gamma A numeric value for the minimum loss reduction required to make a further partition (default is 0).
#' @param colsample_bytree A numeric value specifying the subsample ratio of columns when constructing each tree (default is 1).
#' @param subsample A numeric value specifying the subsample ratio of the training instances (default is 1).
#' @param min_child_weight A numeric value specifying the minimum sum of instance weight needed in a child (default is 1).
#' @param top_n_features An integer specifying the number of top features to display in the importance plot (default is 10).
#' @param verbose An integer specifying the verbosity of the training process (default is 1).
#'
#' @return A list containing:
#' \item{model}{The trained XGBoost model.}
#' \item{confusion_matrix}{The confusion matrix of the test set predictions.}
#' \item{importance}{The feature importance matrix for the top features.}
#' \item{class_mapping}{A named vector showing the mapping from class labels to numeric values used for training.}
#' \item{cv_results}{Cross-validation results, if cross-validation was performed (otherwise NULL).}
#' \item{plot}{A ggplot object showing the feature importance plot.}
#'
#' @details
#' The function allows for training an XGBoost model on cytokine data, splitting the data into training and test sets.
#' If cross-validation is enabled (`cv = TRUE`), it performs k-fold cross-validation and prints the best iteration based on the evaluation metric.
#' The function also visualizes the top N important features using `xgb.ggplot.importance()`.
#'
#' @examples
#' # Example usage:
#' results <- cyt.xgb(data = cytodata, group_col = "Group", train_fraction = 0.8, cv = TRUE)
#' print(results$confusion_matrix)
#' print(results$plot)
#'
#' @import xgboost
#' @import caret
#' @import ggplot2
#' @export
#'
cyt.xgb <- function(data, group_col, train_fraction = 0.7,
                    nrounds = 500, max_depth = 6, eta = 0.1,
                    nfold = 5, cv = FALSE,
                    objective = "multi:softprob", early_stopping_rounds = NULL,
                    eval_metric = "mlogloss",
                    gamma = 0, colsample_bytree = 1, subsample = 1, min_child_weight = 1,
                    top_n_features = 10, verbose = 1) {
  # Ensure the grouping variable is a factor
  data[[group_col]] <- as.factor(data[[group_col]])

  # Create a mapping from group names to numeric labels
  class_labels <- levels(data[[group_col]])
  class_mapping <- setNames(0:(length(class_labels) - 1), class_labels)
  cat("\n### Group to Numeric Label Mapping ###\n")
  print(class_mapping)

  # Convert group column to numeric values using the mapping
  data[[group_col]] <- as.numeric(data[[group_col]]) - 1  # Start from 0 for xgboost

  # Prepare the dataset for xgboost (convert to matrix)
  X <- as.matrix(data[, -which(names(data) == group_col)])
  y <- data[[group_col]]  # Now the numeric class values (0, 1, 2, ..., n-1)

  # Split the data into training and testing sets
  set.seed(123)
  train_indices <- createDataPartition(y, p = train_fraction, list = FALSE)
  X_train <- X[train_indices, ]
  y_train <- y[train_indices]
  X_test <- X[-train_indices, ]
  y_test <- y[-train_indices]

  # Prepare the DMatrix objects for xgboost
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest <- xgb.DMatrix(data = X_test, label = y_test)

  # Set parameters for xgboost
  params <- list(
    objective = objective,  # Multi-class classification
    eval_metric = eval_metric,
    num_class = length(unique(y)),  # Number of classes
    max_depth = max_depth,
    eta = eta,
    gamma = gamma,
    colsample_bytree = colsample_bytree,
    subsample = subsample,
    min_child_weight = min_child_weight
  )

  # Train the XGBoost model
  cat("\n### TRAINING XGBOOST MODEL ###\n")
  xgb_model <- xgb.train(params = params, data = dtrain, nrounds = nrounds, watchlist = list(train = dtrain, test = dtest),
                         early_stopping_rounds = early_stopping_rounds, verbose = verbose)

  # Find the appropriate test evaluation metric dynamically
  eval_col_name <- paste0("test_", eval_metric)

  # Manually track the best iteration based on the minimum evaluation metric
  cat("\nBest iteration from training (based on", eval_metric, "):\n")
  print(xgb_model$evaluation_log[which.min(xgb_model$evaluation_log[[eval_col_name]]), ])

  # Make predictions on the test set
  preds <- predict(xgb_model, X_test)
  pred_labels <- max.col(matrix(preds, ncol = length(unique(y_test)), byrow = TRUE)) - 1  # Convert probability output to labels

  # Confusion matrix
  cat("\n### Confusion Matrix on Test Set ###\n")
  confusion_mat <- confusionMatrix(as.factor(pred_labels), as.factor(y_test))
  print(confusion_mat)

  # Feature importance - show only top_n features
  importance <- xgb.importance(feature_names = colnames(X_train), model = xgb_model)
  top_features <- head(importance, top_n_features)

  # Using xgb.ggplot.importance() for improved feature importance visualization
  cat("\n### Top", top_n_features, "Important Features ###\n")
  print(top_features)
  ggplot_imp <- xgb.ggplot.importance(importance_matrix = top_features, top_n = top_n_features)
  ggplot_imp <- ggplot_imp + ggtitle("Top Features by Gain") +
    ylab("Importance (Gain)") + xlab("Features") +
    theme_minimal()+
    theme(panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"),  # Ensure plot background is white
          legend.background = element_rect(fill = "white", colour = "white"),  # Ensure legend background is white
          axis.title = element_text(color = "black", size = 12, face = "bold"),  # Customize axis titles
          legend.title = element_text(color = "black", size = 10, face = "bold"),  # Customize legend title
          legend.text = element_text(color = "black"))

  print(ggplot_imp)

  # Cross-Validation (optional)
  if (cv) {
    cat("\n### CROSS-VALIDATION USING XGBOOST ###\n")
    xgb_cv <- xgb.cv(params = params, data = dtrain, nrounds = nrounds, nfold = nfold,
                     early_stopping_rounds = early_stopping_rounds, verbose = verbose)

    # Print only the best iteration of the CV results
    cat("\nBest iteration from cross-validation:\n")
    eval_col_name_cv <- paste0("test_", eval_metric,"_mean")
    print(xgb_cv$evaluation_log[which.min(xgb_cv$evaluation_log[[eval_col_name_cv]]), ])
  }

  return(list(model = xgb_model, confusion_matrix = confusion_mat, importance = top_features, class_mapping = class_mapping, cv_results = if (cv) xgb_cv else NULL,
              plot = ggplot_imp))
}
