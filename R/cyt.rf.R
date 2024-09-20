#' Run Random Forest Classification on Cytokine Data
#'
#' This function trains and evaluates a Random Forest classification model on cytokine data.
#' It includes feature importance visualization, cross-validation for feature selection, and performance metrics like accuracy, sensitivity, and specificity.
#' Optionally, for binary classification, the function also plots the ROC curve and calculates the AUC.
#'
#' @param data A data frame containing the cytokine data, with one column as the grouping variable and the rest as numerical features.
#' @param group_col A string representing the name of the column with the grouping variable (i.e., the target variable for classification).
#' @param ntree An integer specifying the number of trees to grow in the forest (default is 500).
#' @param mtry An integer specifying the number of variables randomly selected at each split (default is 5).
#' @param train_fraction A numeric value between 0 and 1 representing the proportion of data to use for training (default is 0.7).
#' @param plot_roc A logical value indicating whether to plot the ROC curve and calculate the AUC for binary classification (default is FALSE).
#' @param k_folds An integer specifying the number of folds for cross-validation (default is 5).
#' @param step A numeric value specifying the fraction of variables to remove at each step during cross-validation for feature selection (default is 0.5).
#' @param run_rfcv A logical value indicating whether to run Random Forest cross-validation for feature selection (default is TRUE).
#'
#' @return A list containing:
#' \item{model}{The trained Random Forest model.}
#' \item{confusion_matrix}{The confusion matrix of the test set predictions.}
#' \item{importance_plot}{A ggplot object showing the variable importance plot based on Mean Decrease Gini.}
#' \item{rfcv_result}{Results from Random Forest cross-validation for feature selection (if `run_rfcv` is TRUE).}
#' \item{importance_data}{A data frame containing the variable importance based on the Gini index.}
#'
#' @details
#' The function fits a Random Forest model to the provided data, splitting it into training and test sets.
#' It calculates key performance metrics such as accuracy, sensitivity, and specificity for both the training and test sets.
#' For binary classification tasks, it can also plot the ROC curve and calculate the AUC.
#' If `run_rfcv` is set to TRUE, Random Forest cross-validation is performed to identify the optimal number of features for classification.
#'
#' @examples
#' # Example usage:
#' results <- cyt.rf(data = cytodata, group_col = "Group", ntree = 500, plot_roc = TRUE)
#' print(results$confusion_matrix)
#' print(results$importance_plot)
#'
#' @import randomForest
#' @import caret
#' @import ggplot2
#' @import pROC
#' @export

cyt.rf <- function(data, group_col, ntree = 500, mtry = 5, train_fraction = 0.7, plot_roc = FALSE,
                   k_folds = 5, step = 0.5, run_rfcv = TRUE) {

  # Ensure the grouping variable is a factor
  data[[group_col]] <- as.factor(data[[group_col]])

  # Split the data into training and testing sets
  set.seed(123) # For reproducibility
  train_indices <- createDataPartition(data[[group_col]], p = train_fraction, list = FALSE)
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]

  # Prepare the formula for the random forest model
  predictors <- setdiff(colnames(data), group_col)
  formula <- as.formula(paste(group_col, "~", paste(predictors, collapse = "+")))

  # Fit the Random Forest model on training data
  rf_model <- randomForest(formula, data = train_data, ntree = ntree, mtry = mtry, importance = TRUE)

  # Print basic Random Forest results with descriptive text
  cat("\n### RANDOM FOREST RESULTS ON TRAINING SET ###\n")
  print(rf_model)

  # Additional metrics from training model
  # Extract confusion matrix from the model
  train_confusion <- rf_model$confusion[, 1:(ncol(rf_model$confusion) - 1)]

  # Calculate overall accuracy
  accuracy <- sum(diag(train_confusion)) / sum(train_confusion)
  cat("\nAccuracy on training set: ", accuracy, "\n")

  # Calculate sensitivity and specificity for each class
  for (i in 1:nrow(train_confusion)) {
    true_positives <- train_confusion[i, i]
    false_negatives <- sum(train_confusion[i, ]) - true_positives
    false_positives <- sum(train_confusion[, i]) - true_positives
    true_negatives <- sum(train_confusion) - (true_positives + false_negatives + false_positives)

    # Sensitivity
    sensitivity <- true_positives / (true_positives + false_negatives)

    # Specificity
    specificity <- true_negatives / (true_negatives + false_positives)

    cat("\nClass '", rownames(train_confusion)[i], "' metrics:\n", sep = "")
    cat("  Sensitivity: ", round(sensitivity, 3), "\n")
    cat("  Specificity: ", round(specificity, 3), "\n")
  }

  # Make predictions on the test set
  rf_pred <- predict(rf_model, newdata = test_data)

  # Confusion matrix on test set
  cat("\n### PREDICTIONS ON TEST SET ###\n")
  confusion_mat <- confusionMatrix(rf_pred, test_data[[group_col]])
  print(confusion_mat$table)

  # Additional metrics from confusion matrix
  cat("\nAccuracy on test set: ", confusion_mat$overall['Accuracy'], "\n")

  # Check if there are two or more classes
  if (length(levels(data[[group_col]])) == 2) {
    # Two-class case: Sensitivity and Specificity are single values
    sensitivity <- confusion_mat$byClass["Sensitivity"]
    specificity <- confusion_mat$byClass["Specificity"]

    cat("\nSensitivity by class:\n")
    cat(paste0("Class: ", levels(data[[group_col]])[1], ": ", round(sensitivity, 3), "\n"))
    cat(paste0("Class: ", levels(data[[group_col]])[2], ": ", round(1 - specificity, 3), "\n"))

    cat("\nSpecificity by class:\n")
    cat(paste0("Class: ", levels(data[[group_col]])[2], ": ", round(specificity, 3), "\n"))
    cat(paste0("Class: ", levels(data[[group_col]])[1], ": ", round(1 - sensitivity, 3), "\n"))

  } else {
    # Multi-class case: Sensitivity and Specificity are vectors
    sensitivity <- confusion_mat$byClass[,"Sensitivity"]
    specificity <- confusion_mat$byClass[,"Specificity"]

    # Print sensitivity and specificity, labeled by class
    cat("\nSensitivity by class:\n")
    for (i in seq_along(sensitivity)) {
      cat(paste0("Class: ", levels(data[[group_col]])[i], ": ", round(sensitivity[i], 4), "\n"))
    }

    cat("\nSpecificity by class:\n")
    for (i in seq_along(specificity)) {
      cat(paste0("Class: ", levels(data[[group_col]])[i], ": ", round(specificity[i], 4), "\n"))
    }
  }
  if (plot_roc && length(unique(test_data[[group_col]])) == 2) {  # Only for binary classification
    rf_prob <- predict(rf_model, newdata = test_data, type = "prob")[, 2]

    # Compute ROC and AUC using pROC package
    roc_obj <- roc(test_data[[group_col]], rf_prob)
    auc_value <- auc(roc_obj)
    cat("\nAUC: ", auc_value, "\n")

    # Create the ROC plot using ggroc()
    roc_plot <- ggroc(roc_obj, color = "blue", linewidth = 1.5, legacy.axes = TRUE) +
      geom_abline(linetype = "dashed", color = "red", linewidth = 1) +  # Random classifier line
      labs(title = "ROC Curve (Test Set)", x = "1 - Specificity", y = "Sensitivity") +
      annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_value, 3)), size = 5, color = "blue") +
      theme_minimal()

    # Display the ROC plot
    print(roc_plot)
  }
  # Extract variable importance based on Mean Decrease Gini and plot with ggplot2
  importance_data <- data.frame(Variable = rownames(importance(rf_model)), Gini = importance(rf_model)[, "MeanDecreaseGini"])

  # Plot Gini-based variable importance (VIP) using ggplot2
  vip_plot <- ggplot(importance_data, aes(x = reorder(Variable, Gini), y = Gini)) +
    geom_bar(stat = "identity", fill = "red2") +
    coord_flip() +
    ggtitle("Variable Importance Plot (Mean Decrease in Gini)") +
    xlab("Features") +
    ylab("Importance (Gini Index)") +
    theme_minimal()+
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"),  # Ensure plot background is white
          legend.background = element_rect(fill = "white", colour = "white"),  # Ensure legend background is white
          axis.title = element_text(color = "black", size = 12, face = "bold"),  # Customize axis titles
          legend.title = element_text(color = "black", size = 10, face = "bold"),  # Customize legend title
          legend.text = element_text(color = "black"))

  print(vip_plot)

  # Optional: Random Forest Cross-Validation for Feature Selection

  if (run_rfcv) {
    cat("\n### RANDOM FOREST CROSS-VALIDATION FOR FEATURE SELECTION ###\n")
    X <- train_data[, predictors]
    y <- train_data[[group_col]]
    rfcv_result <- rfcv(X, y, cv.fold = k_folds, step = step)

    # Convert rfcv results to a data frame for ggplot
    rfcv_data <- data.frame(Variables = rfcv_result$n.var, Error = rfcv_result$error.cv)

    # Plot Cross-Validation Error vs Number of Variables using ggplot2
    rfcv_plot <- ggplot(rfcv_data, aes(x = Variables, y = Error)) +
      geom_line(color = "blue") +
      geom_point(color = "blue") +
      ggtitle("Cross-Validation Error vs. Number of Variables") +
      xlab("Number of Variables") +
      ylab("Cross-Validation Error") +
      scale_x_continuous(breaks = 1:(ncol(train_data)-1) )+  # Ensure all variable counts are showns
      theme_minimal()

    print(rfcv_plot)

    cat("Random Forest CV completed for feature selection. Check the plot for error vs. number of variables.\n")
  } else {
    rfcv_result <- NULL
  }

  return(list(model = rf_model, confusion_matrix = confusion_mat, importance_plot = vip_plot,
              rfcv_result = rfcv_result, importance_data = importance_data, rfcv_data = rfcv_data))
}
