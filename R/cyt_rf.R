#' Run Random Forest Classification on Cytokine Data,
#'
#' This function trains and evaluates a Random Forest classification model on
#' cytokine data. It includes feature importance visualization, cross-
#' validation for feature selection, and performance metrics such as accuracy,
#' sensitivity, and specificity. Optionally, for binary classification, the
#' function also plots the ROC curve and computes the AUC.
#'
#' @param data A data frame containing the cytokine data, with one column as
#'   the grouping variable and the rest as numerical features.
#' @param group_col A string representing the name of the column with the
#'   grouping variable (the target variable for classification).
#' @param ntree An integer specifying the number of trees to grow in the forest
#'   (default is 500).
#' @param mtry An integer specifying the number of variables randomly selected
#'   at each split (default is 5).
#' @param train_fraction A numeric value between 0 and 1 representing the
#'   proportion of data to use for training (default is 0.7).
#' @param plot_roc A logical value indicating whether to plot the ROC curve and
#'   compute the AUC for binary classification (default is FALSE).
#' @param k_folds An integer specifying the number of folds for
#'   cross-validation (default is 5).
#' @param step A numeric value specifying the fraction of variables to remove
#'   at each step during cross-validation for feature selection (default is
#'   0.5).
#' @param run_rfcv A logical value indicating whether to run Random Forest
#'   cross-validation for feature selection (default is TRUE).
#' @param verbose A logical value indicating whether to print additional
#'   informational output to the console. When \code{TRUE}, the function will
#'   display progress messages, and intermediate results when
#'   \code{FALSE} (the default), it runs quietly.
#' @param seed An integer specifying the seed for reproducibility (default is 123).
#'
#' @return A list containing:
#'   \item{model}{The trained Random Forest model.}
#'   \item{confusion_matrix}{The confusion matrix of the test set predictions.}
#'   \item{importance_plot}{A ggplot object showing the variable importance
#'     plot based on Mean Decrease Gini.}
#'   \item{rfcv_result}{Results from Random Forest cross-validation for feature
#'     selection (if \code{run_rfcv} is TRUE).}
#'   \item{importance_data}{A data frame containing the variable importance
#'     based on the Gini index.}
#'
#' @details
#' The function fits a Random Forest model to the provided data by splitting it
#' into training and test sets. It calculates performance metrics such as
#' accuracy, sensitivity, and specificity for both sets. For binary
#' classification, it can also plot the ROC curve and compute the AUC. If
#' \code{run_rfcv} is TRUE, cross-validation is performed to select the optimal
#' number of features.
#' If \code{verbose} is TRUE, the function prints additional information to the
#' console, including training results, test results, and plots.
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
#'
#' @importFrom randomForest randomForest rfcv importance
#' @importFrom caret createDataPartition confusionMatrix
#' @import ggplot2
#' @importFrom pROC roc auc ggroc
#' @importFrom utils capture.output
#' @export
#'

cyt_rf <- function(data, group_col, ntree = 500, mtry = 5,
                   train_fraction = 0.7, plot_roc = FALSE,
                   k_folds = 5, step = 0.5, run_rfcv = TRUE,
                   verbose = FALSE,
                   seed = 123) {
  # Ensure the grouping variable is a factor
  data[[group_col]] <- as.factor(data[[group_col]])

  # Split the data into training and testing sets
  set.seed(seed) # For reproducibility
  train_indices <- caret::createDataPartition(data[[group_col]], p = train_fraction, list = FALSE)
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]

  # Prepare the formula for the random forest model
  predictors <- setdiff(colnames(data), group_col)
  formula_rf <- as.formula(paste(group_col, "~", paste(predictors, collapse = "+")))

  # Fit the Random Forest model on training data
  rf_model <- randomForest::randomForest(formula_rf, data = train_data, ntree = ntree,
                                         mtry = mtry, importance = TRUE, do.trace = FALSE)

  # Print basic Random Forest results for training set if verbose is TRUE
  if (verbose) {
    cat("\n### RANDOM FOREST RESULTS ON TRAINING SET ###\n")
    cat(paste(capture.output(rf_model), collapse = "\n"), "\n")
  }

  # Calculate training metrics (from training confusion matrix)
  train_confusion <- rf_model$confusion[, 1:(ncol(rf_model$confusion) - 1)]
  accuracy_train <- sum(diag(train_confusion)) / sum(train_confusion)
  if (verbose) {
    cat("\nAccuracy on training set: ", round(accuracy_train, 3), "\n")
  }
  for (i in seq_len(nrow(train_confusion))) {
    true_positives <- train_confusion[i, i]
    false_negatives <- sum(train_confusion[i, ]) - true_positives
    false_positives <- sum(train_confusion[, i]) - true_positives
    true_negatives <- sum(train_confusion) - (true_positives + false_negatives + false_positives)
    sensitivity <- true_positives / (true_positives + false_negatives)
    specificity <- true_negatives / (true_negatives + false_positives)
    if (verbose) {
      cat("\nTraining set - Class '", rownames(train_confusion)[i], "' metrics:\n", sep = "")
      cat("  Sensitivity: ", round(sensitivity, 3), "\n")
      cat("  Specificity: ", round(specificity, 3), "\n")
    }
  }

  # Make predictions on the test set and calculate confusion matrix
  rf_pred <- predict(rf_model, newdata = test_data)
  if (verbose) {
    cat("\n### PREDICTIONS ON TEST SET ###\n")
  }
  confusion_mat <- caret::confusionMatrix(rf_pred, test_data[[group_col]])
  if (verbose) {
    cat(paste(capture.output(print(confusion_mat$table)), collapse = "\n"), "\n")
    cat("\nAccuracy on test set: ", round(confusion_mat$overall["Accuracy"], 3), "\n")
  }

  # Report sensitivity and specificity on the test set
  if (length(levels(data[[group_col]])) == 2) {
    # Two-class case: Sensitivity and Specificity are single values
    sensitivity_val <- confusion_mat$byClass["Sensitivity"]
    specificity_val <- confusion_mat$byClass["Specificity"]
    if (verbose) {
      cat("\nSensitivity by class:\n")
      cat(paste0("Class: ", levels(data[[group_col]])[1], ": ", round(sensitivity_val, 3), "\n"))
      cat(paste0("Class: ", levels(data[[group_col]])[2], ": ", round(1 - specificity_val, 3), "\n"))
      cat("\nSpecificity by class:\n")
      cat(paste0("Class: ", levels(data[[group_col]])[2], ": ", round(specificity_val, 3), "\n"))
      cat(paste0("Class: ", levels(data[[group_col]])[1], ": ", round(1 - sensitivity_val, 3), "\n"))
    }
  } else {
    # Multi-class case: Sensitivity and Specificity are vectors
    sensitivity_vec <- confusion_mat$byClass[, "Sensitivity"]
    specificity_vec <- confusion_mat$byClass[, "Specificity"]
    if (verbose) {
      cat("\nSensitivity by class:\n")
      for (i in seq_along(sensitivity_vec)) {
        cat(paste0("Class: ", levels(data[[group_col]])[i], ": ", round(sensitivity_vec[i], 3), "\n"))
      }
      cat("\nSpecificity by class:\n")
      for (i in seq_along(specificity_vec)) {
        cat(paste0("Class: ", levels(data[[group_col]])[i], ": ", round(specificity_vec[i], 3), "\n"))
      }
    }
  }

  # ROC curve and AUC for binary classification, if requested
  roc_plot <- NULL
  if (plot_roc && length(levels(data[[group_col]])) == 2) {
    rf_prob <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    roc_obj <- pROC::roc(test_data[[group_col]], rf_prob, quiet = TRUE)
    auc_value <- pROC::auc(roc_obj)
    if (verbose) {
      cat("\nAUC: ", round(auc_value, 3), "\n")
    }
    roc_plot <- pROC::ggroc(roc_obj, color = "blue", linewidth = 1.5, legacy.axes = TRUE) +
      ggplot2::geom_abline(linetype = "dashed", color = "red", linewidth = 1) +
      ggplot2::labs(
        title = "ROC Curve (Test Set)",
        x = "1 - Specificity",
        y = "Sensitivity"
      ) +
      ggplot2::annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_value, 3)),
                        size = 5, color = "blue") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background  = ggplot2::element_rect(fill = "white", color = NA),
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_line(color = "grey95")
      )
    print(roc_plot)
  }

  # Extract variable importance and generate plot
  importance_data <- data.frame(Variable = rownames(randomForest::importance(rf_model)),
                                Gini = randomForest::importance(rf_model)[, "MeanDecreaseGini"])
  vip_plot <- ggplot2::ggplot(importance_data, ggplot2::aes(x = reorder(Variable, Gini), y = Gini)) +
    ggplot2::geom_bar(stat = "identity", fill = "red2") +
    ggplot2::coord_flip() +
    ggplot2::ggtitle("Variable Importance Plot (Mean Decrease in Gini)") +
    ggplot2::xlab("Features") +
    ggplot2::ylab("Importance (Gini Index)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "white", colour = "white"),
      plot.background = ggplot2::element_rect(fill = "white", colour = "white"),
      legend.background = ggplot2::element_rect(fill = "white", colour = "white"),
      axis.title = ggplot2::element_text(color = "black", size = 12, face = "bold"),
      legend.title = ggplot2::element_text(color = "black", size = 10, face = "bold"),
      legend.text = ggplot2::element_text(color = "black")
    )
  print(vip_plot)

  # Optional: Random Forest Cross-Validation for Feature Selection
  rfcv_data <- NULL
  rfcv_result <- NULL
  rfcv_plot <- NULL
  if (run_rfcv) {
    if (verbose) {
      cat("\n### RANDOM FOREST CROSS-VALIDATION FOR FEATURE SELECTION ###\n")
    }
    x_train <- train_data[, predictors]
    y_train <- train_data[[group_col]]
    rfcv_result <- randomForest::rfcv(x_train, y_train, cv.fold = k_folds, step = step, do.trace = FALSE)
    rfcv_data <- data.frame(Variables = rfcv_result$n.var, Error = rfcv_result$error.cv)
    rfcv_plot <- ggplot2::ggplot(rfcv_data, ggplot2::aes(x = Variables, y = Error)) +
      ggplot2::geom_line(color = "blue") +
      ggplot2::geom_point(color = "blue") +
      ggplot2::ggtitle("Cross-Validation Error vs. Number of Variables") +
      ggplot2::xlab("Number of Variables") +
      ggplot2::ylab("Cross-Validation Error") +
      ggplot2::scale_x_continuous(breaks = 1:(ncol(train_data) - 1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background  = ggplot2::element_rect(fill = "white", color = NA),
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_line(color = "grey95")
      )
    print(rfcv_plot)
    if (verbose) {
      cat("Random Forest CV completed for feature selection. Check the plot for error vs. number of variables.\n")
    }
  }

  # return(list(
  #   model = rf_model,
  #   confusion_matrix = confusion_mat,
  #   importance_plot = vip_plot,
  #   rfcv_result = rfcv_result,
  #   rfcv_plot = rfcv_plot,
  #   importance_data = importance_data,
  #   rfcv_data = rfcv_data,
  #   roc_plot = roc_plot
  # ))
}

