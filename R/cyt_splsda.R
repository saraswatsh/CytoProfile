#' Analyze data with Sparse Partial Least Squares Discriminant Analysis
#' (sPLS-DA).
#'
#' @param data A matrix or data frame containing the variables. Columns not
#'   specified by \code{group_col} or \code{group_col2} are assumed to be continuous
#'   variables for analysis.
#' @param group_col A string specifying the column name that contains the first group
#'   information. If \code{group_col2} is not provided, an overall analysis will
#'   be performed.
#' @param group_col2 A string specifying the second grouping column. Default is
#'   \code{NULL}.
#' @param colors A vector of colors for the groups or treatments. If
#'   \code{NULL}, a random palette (using \code{rainbow}) is generated based on
#'   the number of groups.
#' @param pdf_title A string specifying the file name for saving the PDF output.
#' @param ellipse Logical. Whether to draw a 95\% confidence ellipse on the
#'   figures. Default is \code{FALSE}.
#' @param bg Logical. Whether to draw the prediction background in the figures.
#'   Default is \code{FALSE}.
#' @param conf_mat Logical. Whether to print the confusion matrix for the
#'   classifications. Default is \code{FALSE}.
#' @param var_num Numeric. The number of variables to be used in the PLS-DA model.
#' @param cv_opt Character. Option for cross-validation method: either
#'   "loocv" or "Mfold". Default is \code{NULL}.
#' @param fold_num Numeric. The number of folds to use if \code{cv_opt} is
#'   "Mfold". Default is 5.
#' @param scale Character. Option for data transformation; if set to
#'   \code{"log2"}, a log2 transformation is applied to the continuous
#'   variables. Default is \code{NULL}.
#' @param comp_num Numeric. The number of components to calculate in the sPLS-DA
#'   model. Default is 2.
#' @param pch_values A vector of integers specifying the plotting characters
#'   (pch values) to be used in the plots.
#' @param style Character. If set to \code{"3D"} or \code{"3d"} and
#'   \code{comp_num} equals 3, a 3D plot is generated using the
#'   \code{plot3D} package. Default is \code{NULL}.
#' @param roc Logical. Whether to compute and plot the ROC curve for the model.
#'   Default is \code{FALSE}.
#'
#' @description
#' This function conducts Sparse Partial Least Squares Discriminant Analysis
#' (sPLS-DA) on the provided data. It uses the specified \code{group_col} (and
#' optionally \code{group_col2}) to define class labels while assuming the remaining
#' columns contain continuous variables. The function supports a log2
#' transformation via the \code{scale} parameter and generates a series of plots,
#' including classification plots, scree plots, loadings plots, and VIP score
#' plots. Optionally, ROC curves are produced when \code{roc} is \code{TRUE}.
#' Additionally, cross-validation is supported via LOOCV or Mfold methods. When
#' both \code{group_col} and \code{group_col2} are provided and differ, the function
#' analyzes each treatment level separately.
#'
#' @return A PDF file containing the classification figures, component figures
#'   with Variable of Importance in Projection (VIP) scores, and classifications
#'   based on VIP scores greater than 1. ROC curves and confusion matrices are also
#'   produced if requested.
#'
#' @examples
#' # Loading Sample Data
#' data_df <- ExampleData1[,-c(3)]
#' data_df <- dplyr::filter(data_df, Group != "ND", Treatment != "Unstimulated")
#'
#' cyt_splsda(data_df, pdf_title = "Example sPLS-DA Analysis.pdf",
#' colors = c("black", "purple"), bg = FALSE, scale = "log2",
#' conf_mat = FALSE, var_num = 25, cv_opt = NULL, comp_num = 2,
#' pch_values = c(16, 4), style = "3d", ellipse = TRUE,
#' group_col = "Group", group_col2 = "Treatment", roc = FALSE)
#'
#' @export
#' @importFrom mixOmics splsda background.predict perf vip auroc plotIndiv plotLoadings
#' @import ggplot2
#' @importFrom plot3D scatter3D
#' @importFrom reshape2 melt
#' @importFrom caret confusionMatrix

cyt_splsda <- function(data, group_col = NULL, group_col2 = NULL, colors = NULL,
                      pdf_title, ellipse = FALSE, bg = FALSE, conf_mat = FALSE,
                      var_num, cv_opt = NULL, fold_num = 5, scale = NULL,
                      comp_num = 2, pch_values, style = NULL, roc = FALSE) {
  # If one factor is missing, use the provided column for
  # both grouping and treatment.
  if (!is.null(group_col) && is.null(group_col2)) {
    message("No second grouping column provided; performing overall analysis.")
    group_col2 <- group_col
  }
  if(is.null(group_col) && !is.null(group_col2)) {
    stop("No first grouping column provided; must provide the first grouping column.")
  }
  if (is.null(group_col) && is.null(group_col2)) {
    stop("At least one grouping column must be provided.")
  }

  # Optionally apply log2 transformation
  if (!is.null(scale) && scale == "log2") {
    data <- data.frame(
      data[, c(group_col, group_col2)],
      log2(data[, !(names(data) %in% c(group_col, group_col2))])
    )
    message("Results based on log2 transformation:")
  } else if (is.null(scale)) {
    message("Results based on no transformation:")
  }

  # Extract the grouping variable from your data (using group_col or group_col2)
  # Extract grouping variable(s)
  if (group_col == group_col2) {
    group_vec <- data[[group_col]]
  } else {
    # Combine the two grouping columns into a composite factor
    group_vec <- data[[group_col2]]
  }

  # Now perform the check for pch_values:
  if (is.null(pch_values)) {
    stop("Please enter a vector of pch values, e.g. c(16, 4).")
  }
  if (group_col == group_col2) {
    if (length(pch_values) < length(unique(data[[group_col]]))) {
      stop("Please ensure the number of pch values provided (", length(pch_values),
           ") is at least equal to the number of unique groups (", length(unique(data[[group_col]])),
           ") from the grouping column.")
    }
  } else {
    # When group_col and group_col2 differ, use the levels of group_col for pch
    if (length(pch_values) < length(unique(data[[group_col]]))) {
      stop("Please ensure the number of pch values provided (", length(pch_values),
           ") is at least equal to the number of unique groups (", length(unique(data[[group_col]])),
           ") from the first grouping column.")
    }
  }

  # Generate a color palette if not provided (based on the
  # grouping variable levels in the entire dataset)
  if (is.null(colors)) {
    num_groups <- length(unique(data[[group_col]]))
    colors <- rainbow(num_groups)
  }

  pdf(file = pdf_title, width = 8.5, height = 8)

  # Case 1: Only one factor provided (both columns are the same)
  if (group_col == group_col2) {
    overall_analysis <- "Overall Analysis"

    # Remove the factor column from predictors and keep only numeric columns
    the_data_df <- data[, !(names(data) %in% c(group_col))]
    the_data_df <- the_data_df[, sapply(the_data_df, is.numeric)]

    the_groups <- as.vector(data[[group_col]])
    if (length(unique(the_groups)) < 2) {
      stop("The grouping variable must have at least two levels for PLS-DA.
           Please provide an appropriate grouping column.")
    }

    cytokine_splsda <- mixOmics::splsda(the_data_df, the_groups,
      scale = TRUE, ncomp = comp_num,
      keepX = rep(var_num, comp_num)
    )

    splsda_predict <- predict(cytokine_splsda, the_data_df, dist = "max.dist")
    prediction1 <- cbind(original = the_groups, splsda_predict$class$max.dist)
    accuracy1 <- (sum(prediction1[, 1] == prediction1[, 3]) /
                    length(prediction1[, 1]))
    acc1 <- 100 * signif(accuracy1, digits = 2)

    bg_maxdist <- mixOmics::background.predict(cytokine_splsda,
                                                comp.predicted = 2,
                                                dist = "max.dist",
                                                xlim = c(-15,15),
                                                ylim = c(-15,15)
    )
    group_factors <- seq_len(length(levels(factor(the_groups))))

    plot_args <- list(cytokine_splsda,
      ind.names = NA, legend = TRUE, col = colors,
      pch = pch_values, pch.levels = group_factors,
      title = paste(
        overall_analysis, "With Accuracy:",
        acc1, "%"
      ),
      legend.title = group_col
    )
    if (ellipse) plot_args$ellipse <- TRUE
    if (bg) plot_args$background <- bg_maxdist
    do.call(mixOmics::plotIndiv, plot_args)

    if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
      cytokine_scores <- cytokine_splsda$variates$X
      plot3D::scatter3D(cytokine_scores[, 1], cytokine_scores[, 2],
        cytokine_scores[, 3],
        pch = pch_values, col = colors,
        xlab = "Component 1", ylab = "Component 2",
        zlab = "Component 3",
        main = paste("3D Plot:", overall_analysis),
        theta = 20, phi = 30, bty = "g", colkey = FALSE
      )
    }

    # If roc = TRUE, compute and plot ROC curve for the overall model
    if (roc) {
      roc_obj <- mixOmics::auroc(
        object = cytokine_splsda, newdata = the_data_df,
        outcome.test = the_groups,
        plot = TRUE, roc.comp = comp_num,
        title = paste0("ROC Curve:", overall_analysis), print = FALSE
      )
    }

    # Cross-validation methods
    if (!is.null(cv_opt)) {
      if (cv_opt == "loocv") {
        set.seed(123)
        loocv_results <- mixOmics::perf(cytokine_splsda, validation = "loo")
        loocv_error_rate <- loocv_results$error.rate$overall[
          "comp2",
          "max.dist"
        ]
        loocv_acc <- 100 * signif(1 - loocv_error_rate, digits = 2)
        print(paste0("LOOCV Accuracy: ", loocv_acc, "%"))

        error_rates <- loocv_results$error.rate$overall[, "max.dist"]
        error_df <- as.data.frame(error_rates)
        error_df$Component <- rownames(error_df)
        error_df <- reshape2::melt(error_df,
          id.vars = "Component",
          variable.name = "Distance",
          value.name = "ErrorRate"
        )

        a <- ggplot2::ggplot(error_df, aes(
          x = Component, y = ErrorRate,
          color = Distance, group = 1
        )) +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 3) +
          ggplot2::labs(
            title = paste("LOOCV Error Rate:", overall_analysis),
            x = "Number of Components",
            y = "Error Rate"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggplot2::scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      } else if (cv_opt == "Mfold") {
        set.seed(123)
        fold_results <- mixOmics::perf(cytokine_splsda,
          validation = "Mfold",
          folds = fold_num, nrepeat = 1000
        )
        fold_error_rate <- fold_results$error.rate$overall[
          "comp2",
          "max.dist"
        ]
        fold_acc <- 100 * signif(1 - fold_error_rate, digits = 2)
        print(paste0("Mfold Accuracy: ", fold_acc, "%"))

        error_rates <- fold_results$error.rate$overall[, "max.dist"]
        error_df <- as.data.frame(error_rates)
        error_df$Component <- rownames(error_df)
        error_df <- reshape2::melt(error_df,
          id.vars = "Component",
          variable.name = "Distance",
          value.name = "ErrorRate"
        )

        a <- ggplot2::ggplot(error_df, aes(
          x = Component, y = ErrorRate,
          color = Distance, group = 1
        )) +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 3) +
          ggplot2::labs(
            title = paste("Mfold Error Rate:", overall_analysis),
            x = "Number of Components",
            y = "Error Rate"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggplot2::scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      }
    }

    # Loadings plot for each component
    for (comp in 1:comp_num) {
      mixOmics::plotLoadings(cytokine_splsda,
        comp = comp, contrib = "max", method = "mean",
        size.name = 1, size.legend = 1, legend.color = colors,
        title = paste("Component", comp, ":", overall_analysis),
        size.title = 1, legend.title = group_col
      )
    }

    # VIP scores and plot for PLS-DA with VIP > 1
    all_vip_scores <- mixOmics::vip(cytokine_splsda)
    for (comp in 1:comp_num) {
      vscore <- as.data.frame(all_vip_scores[, comp, drop = FALSE])
      vscore$metabo <- rownames(vscore)
      vscore$comp <- vscore[, 1]
      bar <- vscore[, c("metabo", "comp")]
      bar <- bar[order(bar$comp, decreasing = TRUE), ]

      a <- ggplot2::ggplot(bar, aes(x = metabo, y = comp)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::scale_y_continuous(limits = c(0, max(bar$comp))) +
        ggplot2::geom_hline(yintercept = 1, color = "grey") +
        ggplot2::scale_x_discrete(limits = factor(bar$metabo)) +
        ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
        ggplot2::labs(x = "", y = "VIP score") +
        ggplot2::ggtitle(paste("Component", comp)) +
        ggplot2::theme(
          panel.grid = element_blank(),
          panel.background = element_rect(
            color = "black",
            fill = "transparent"
          )
        )
      print(a)
    }

    # PLS-DA on VIP > 1: Subset predictors with VIP > 1
    condt_variable <- all_vip_scores[, 1] > 1
    keep_x <- sum(condt_variable)
    the_data_mat <- the_data_df[, condt_variable, drop = FALSE]
    cytokine_splsda2 <- mixOmics::splsda(the_data_mat, the_groups,
      scale = TRUE, ncomp = comp_num,
      keepX = rep(keep_x, comp_num)
    )

    splsda_predict2 <- predict(cytokine_splsda2, the_data_mat,
      dist = "max.dist"
    )
    prediction2 <- cbind(original = the_groups, splsda_predict2$class$max.dist)
    accuracy2 <- (sum(prediction2[, 1] == prediction2[, 3]) /
                    length(prediction2[, 1]))
    acc2 <- 100 * signif(accuracy2, digits = 2)

    # Create a grid of values
    bg_maxdist2 <- mixOmics::background.predict(cytokine_splsda2,
                                                comp.predicted = 2,
                                                dist = "max.dist",
                                                xlim = c(-15,15),
                                                ylim = c(-15,15)
    )

    plot_args2 <- list(cytokine_splsda2,
      ind.names = NA, legend = TRUE, col = colors,
      pch = pch_values, pch.levels = group_factors,
      title = paste(
        overall_analysis, "(VIP>1)",
        "With Accuracy:", acc2, "%"
      ),
      legend.title = group_col
    )
    if (ellipse) plot_args2$ellipse <- TRUE
    if (bg) plot_args2$background <- bg_maxdist2
    do.call(mixOmics::plotIndiv, plot_args2)

    if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
      cytokine_scores2 <- cytokine_splsda2$variates$X
      plot3D::scatter3D(cytokine_scores2[, 1], cytokine_scores2[, 2],
        cytokine_scores2[, 3],
        pch = pch_values, col = colors,
        xlab = "Component 1", ylab = "Component 2",
        zlab = "Component 3",
        main = paste("3D Plot:", overall_analysis, "(VIP>1)"),
        theta = 20, phi = 30, bty = "g", colkey = FALSE
      )
    }
    # Loadings plot for each component with VIP > 1
    for (comp in 1:comp_num) {
      mixOmics::plotLoadings(cytokine_splsda2,
                             comp = comp, contrib = "max", method = "mean",
                             size.name = 1, size.legend = 1, legend.color = colors,
                             title = paste("Component", comp, "(VIP > 1):", overall_analysis),
                             size.title = 1, legend.title = group_col
      )
    }
    if (!is.null(cv_opt)) {
      if (cv_opt == "loocv") {
        set.seed(123)
        loocv_results2 <- mixOmics::perf(cytokine_splsda2, validation = "loo")
        loocv_error_rate2 <- loocv_results2$error.rate$overall[
          "comp2",
          "max.dist"
        ]
        loocv_acc2 <- 100 * signif(1 - loocv_error_rate2, digits = 2)
        print(paste0("LOOCV Accuracy (VIP>1): ", loocv_acc2, "%"))

        error_rates2 <- loocv_results2$error.rate$overall[, "max.dist"]
        error_df2 <- as.data.frame(error_rates2)
        error_df2$Component <- rownames(error_df2)
        error_df2 <- reshape2::melt(error_df2,
          id.vars = "Component",
          variable.name = "Distance",
          value.name = "ErrorRate"
        )

        a <- ggplot2::ggplot(error_df2, aes(
          x = Component, y = ErrorRate, color =
            Distance, group = 1
        )) +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 3) +
          ggplot2::labs(
            title = paste("LOOCV Error Rate (VIP>1):", overall_analysis),
            x = "Number of Components",
            y = "Error Rate"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggplot2::scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      } else if (cv_opt == "Mfold") {
        set.seed(123)
        fold_results2 <- mixOmics::perf(cytokine_splsda2,
          validation = "Mfold",
          folds = fold_num, nrepeat = 1000
        )
        fold_error_rate2 <- fold_results2$error.rate$overall[
          "comp2",
          "max.dist"
        ]
        fold_acc2 <- 100 * signif(1 - fold_error_rate2, digits = 2)
        print(paste0("Mfold Accuracy (VIP>1): ", fold_acc2, "%"))

        error_rates2 <- fold_results2$error.rate$overall[, "max.dist"]
        error_df2 <- as.data.frame(error_rates2)
        error_df2$Component <- rownames(error_df2)
        error_df2 <- reshape2::melt(error_df2,
          id.vars = "Component",
          variable.name = "Distance",
          value.name = "ErrorRate"
        )

        a <- ggplot2::ggplot(error_df2, aes(
          x = Component, y = ErrorRate,
          color = Distance, group = 1
        )) +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 3) +
          ggplot2::labs(
            title = paste("Mfold Error Rate (VIP>1):", overall_analysis),
            x = "Number of Components",
            y = "Error Rate"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggplot2::scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      }
    }

    if (conf_mat == TRUE) {
      cat("Overall Confusion Matrix for PLS-DA Comparison\n")

      # Confusion Matrix for main model
      cm <- caret::confusionMatrix(
        data = as.factor(prediction1[, 3]), # predicted
        reference = as.factor(prediction1[, 1]) # actual
      )
      print(cm$table)
      cat("Accuracy:", signif(cm$overall["Accuracy"], 2), "\n")

      # Check if binary or multi-class
      if (nlevels(as.factor(prediction1[, 1])) == 2) {
        cat("Sensitivity:", signif(cm$byClass["Sensitivity"], 2), "\n")
        cat("Specificity:", signif(cm$byClass["Specificity"], 2), "\n")
      } else {
        cat("\nPer-Class Sensitivity:\n")
        print(signif(cm$byClass[, "Sensitivity"], 2))
        cat("\nPer-Class Specificity:\n")
        print(signif(cm$byClass[, "Specificity"], 2))
        macro_sens <- mean(cm$byClass[, "Sensitivity"], na.rm = TRUE)
        macro_spec <- mean(cm$byClass[, "Specificity"], na.rm = TRUE)
        cat("\nMacro-Averaged Sensitivity:", signif(macro_sens, 2), "\n")
        cat("Macro-Averaged Specificity:", signif(macro_spec, 2), "\n")
      }

      cat("Overall Confusion Matrix for PLS-DA Comparison with VIP Score > 1\n")

      # Confusion Matrix for VIP>1 model
      cm_vip <- caret::confusionMatrix(
        data = as.factor(prediction2[, 3]), # predicted
        reference = as.factor(prediction2[, 1]) # actual
      )
      print(cm_vip$table)
      cat("Accuracy:", signif(cm_vip$overall["Accuracy"], 2), "\n")
      if (nlevels(as.factor(prediction2[, 1])) == 2) {
        cat("Sensitivity:", signif(cm_vip$byClass["Sensitivity"], 2), "\n")
        cat("Specificity:", signif(cm_vip$byClass["Specificity"], 2), "\n")
      } else {
        cat("\nPer-Class Sensitivity:\n")
        print(signif(cm_vip$byClass[, "Sensitivity"], 2))
        cat("\nPer-Class Specificity:\n")
        print(signif(cm_vip$byClass[, "Specificity"], 2))
        macro_sens_vip <- mean(cm_vip$byClass[, "Sensitivity"], na.rm = TRUE)
        macro_spec_vip <- mean(cm_vip$byClass[, "Specificity"], na.rm = TRUE)
        cat("\nMacro-Averaged Sensitivity:", signif(macro_sens_vip, 2), "\n")
        cat("Macro-Averaged Specificity:", signif(macro_spec_vip, 2), "\n")
      }
    }
    # If roc = TRUE, compute and plot ROC curve for the overall model of VIP > 1
    if (roc) {
      roc_obj2 <- mixOmics::auroc(
        object = cytokine_splsda2, newdata = the_data_mat,
        outcome.test = the_groups,
        plot = TRUE, roc.comp = comp_num,
        title = paste0("ROC Curve (VIP>1):", overall_analysis), print = FALSE
      )
    }
  } else {
    # Case 2: Both group and treatment columns are provided and they differ.
    levels_vec <- unique(data[[group_col2]])
    for (i in seq_along(levels_vec)) {
      current_level <- levels_vec[i]
      overall_analysis <- current_level
      condt <- data[[group_col2]] == current_level

      the_data_df <- data[condt, -which(names(data) %in% c(
        group_col,
        group_col2
      ))]
      the_data_df <- the_data_df[, sapply(the_data_df, is.numeric)]
      the_groups <- as.vector(data[condt, group_col])

      if (length(unique(the_groups)) < 2) {
        stop("The grouping variable must have at least two levels for PLS-DA.
             Please provide an appropriate grouping column.")
      }

      cytokine_splsda <- mixOmics::splsda(the_data_df, the_groups,
        scale = TRUE, ncomp = comp_num,
        keepX = rep(var_num, comp_num)
      )

      splsda_predict <- predict(cytokine_splsda, the_data_df,
        dist = "max.dist"
      )
      prediction1 <- cbind(
        original = the_groups,
        splsda_predict$class$max.dist
      )
      accuracy1 <- (sum(prediction1[, 1] == prediction1[, 3]) /
                      length(prediction1[, 1]))
      acc1 <- 100 * signif(accuracy1, digits = 2)

      # Create a grid of values
      bg_maxdist <- mixOmics::background.predict(cytokine_splsda,
                                                 comp.predicted = 2,
                                                 dist = "max.dist",
                                                 xlim = c(-15,15),
                                                 ylim = c(-15,15)
      )

      group_factors <- seq_len(length(levels(factor(the_groups))))

      plot_args <- list(cytokine_splsda,
        ind.names = NA, legend = TRUE, col = colors,
        pch = pch_values[group_factors], pch.levels = pch_values[group_factors],
        title = paste(
          overall_analysis, "With Accuracy:",
          acc1, "%"
        ),
        legend.title = group_col
      )
      if (ellipse) plot_args$ellipse <- TRUE
      if (bg) plot_args$background <- bg_maxdist
      do.call(mixOmics::plotIndiv, plot_args)

      if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
        cytokine_scores <- cytokine_splsda$variates$X
        plot3D::scatter3D(cytokine_scores[, 1], cytokine_scores[, 2],
          cytokine_scores[, 3],
          pch = pch_values, col = colors,
          xlab = "Component 1", ylab = "Component 2",
          zlab = "Component 3",
          main = paste("3D Plot:", overall_analysis),
          theta = 20, phi = 30, bty = "g", colkey = FALSE
        )
      }

      # If roc = TRUE, compute and plot ROC curve for the overall model
      if (roc) {
        roc_obj <- mixOmics::auroc(
          object = cytokine_splsda, newdata = the_data_df, outcome.test =
            the_groups,
          plot = TRUE, roc.comp = comp_num,
          title = paste0("ROC Curve:", overall_analysis), print = FALSE
        )
      }

      if (!is.null(cv_opt)) {
        if (cv_opt == "loocv") {
          set.seed(123)
          loocv_results <- mixOmics::perf(cytokine_splsda, validation = "loo")
          loocv_error_rate <- loocv_results$error.rate$overall[
            "comp2",
            "max.dist"
          ]
          loocv_acc <- 100 * signif(1 - loocv_error_rate, digits = 2)
          print(paste0(current_level, " LOOCV Accuracy: ", loocv_acc, "%"))

          error_rates <- loocv_results$error.rate$overall[, "max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df,
            id.vars = "Component",
            variable.name = "Distance",
            value.name = "ErrorRate"
          )

          a <- ggplot2::ggplot(error_df, aes(
            x = Component, y = ErrorRate,
            color = Distance, group = 1
          )) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 3) +
            ggplot2::labs(
              title = paste("LOOCV Error Rate:", overall_analysis),
              x = "Number of Components",
              y = "Error Rate"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ggplot2::scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        } else if (cv_opt == "Mfold") {
          set.seed(123)
          fold_results <- mixOmics::perf(cytokine_splsda,
            validation = "Mfold",
            folds = fold_num, nrepeat = 1000
          )
          fold_error_rate <- fold_results$error.rate$overall[
            "comp2",
            "max.dist"
          ]
          fold_acc <- 100 * signif(1 - fold_error_rate, digits = 2)
          print(paste0(current_level, " Mfold Accuracy: ", fold_acc, "%"))

          error_rates <- fold_results$error.rate$overall[, "max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df,
            id.vars = "Component",
            variable.name = "Distance",
            value.name = "ErrorRate"
          )

          a <- ggplot2::ggplot(error_df, aes(
            x = Component, y = ErrorRate,
            color = Distance, group = 1
          )) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 3) +
            ggplot2::labs(
              title = paste("Mfold Error Rate:", overall_analysis),
              x = "Number of Components",
              y = "Error Rate"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ggplot2::scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
      }

      # Loadings plot for each component
      for (comp in 1:comp_num) {
        mixOmics::plotLoadings(cytokine_splsda,
          comp = comp, contrib = "max", method = "mean",
          size.name = 1, size.legend = 1, legend.color = colors,
          title = paste("Component", comp, ":", overall_analysis),
          size.title = 1, legend.title = group_col
        )
      }

      all_vip_scores <- mixOmics::vip(cytokine_splsda)
      for (comp in 1:comp_num) {
        vscore <- as.data.frame(all_vip_scores[, comp, drop = FALSE])
        vscore$metabo <- rownames(vscore)
        vscore$comp <- vscore[, 1]
        bar <- vscore[, c("metabo", "comp")]
        bar <- bar[order(bar$comp, decreasing = TRUE), ]

        a <- ggplot2::ggplot(bar, aes(x = metabo, y = comp)) +
          ggplot2::geom_bar(stat = "identity", position = "dodge") +
          ggplot2::scale_y_continuous(limits = c(0, max(bar$comp))) +
          ggplot2::geom_hline(yintercept = 1, color = "grey") +
          ggplot2::scale_x_discrete(limits = factor(bar$metabo)) +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
          ggplot2::labs(x = "", y = "VIP score") +
          ggplot2::ggtitle(paste("Component", comp)) +
          ggplot2::theme(
            panel.grid = element_blank(),
            panel.background = element_rect(
              color = "black",
              fill = "transparent"
            )
          )
        print(a)
      }

      condt_variable <- all_vip_scores[, 1] > 1
      keep_x <- sum(condt_variable)
      the_data_mat <- the_data_df[, condt_variable, drop = FALSE]
      cytokine_splsda2 <- mixOmics::splsda(the_data_mat, the_groups,
        scale = TRUE, ncomp = comp_num,
        keepX = rep(keep_x, comp_num)
      )

      splsda_predict2 <- predict(cytokine_splsda2, the_data_mat,
        dist = "max.dist"
      )
      prediction2 <- cbind(
        original = the_groups,
        splsda_predict2$class$max.dist
      )
      accuracy2 <- (sum(prediction2[, 1] == prediction2[, 3]) /
                      length(prediction2[, 1]))
      acc2 <- 100 * signif(accuracy2, digits = 2)

      # Create a grid of values
      bg_maxdist2 <- mixOmics::background.predict(cytokine_splsda2,
                                                 comp.predicted = 2,
                                                 dist = "max.dist",
                                                 xlim = c(-15,15),
                                                 ylim = c(-15,15)
      )

      plot_args2 <- list(cytokine_splsda2,
        ind.names = NA, legend = TRUE, col = colors,
        pch = pch_values[group_factors], pch.levels = pch_values[group_factors],
        title = paste(
          overall_analysis, "(VIP>1)",
          "With Accuracy:", acc2, "%"
        ),
        legend.title = group_col
      )
      if (ellipse) plot_args2$ellipse <- TRUE
      if (bg) plot_args2$background <- bg_maxdist2
      do.call(mixOmics::plotIndiv, plot_args2)

      if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
        cytokine_scores2 <- cytokine_splsda2$variates$X
        plot3D::scatter3D(cytokine_scores2[, 1], cytokine_scores2[, 2],
          cytokine_scores2[, 3],
          pch = pch_values, col = colors,
          xlab = "Component 1", ylab = "Component 2",
          zlab = "Component 3",
          main = paste("3D Plot:", overall_analysis, "(VIP>1)"),
          theta = 20, phi = 30, bty = "g", colkey = FALSE
        )
      }
      # Loadings plot for each component with VIP > 1
      # Loadings plot for each component
      for (comp in 1:comp_num) {
        mixOmics::plotLoadings(cytokine_splsda2,
                               comp = comp, contrib = "max", method = "mean",
                               size.name = 1, size.legend = 1, legend.color = colors,
                               title = paste("Component", comp, "(VIP > 1):", overall_analysis),
                               size.title = 1, legend.title = group_col
        )
      }
      if (roc) {
        roc_obj2 <- mixOmics::auroc(
          object = cytokine_splsda2, newdata = the_data_mat,
          outcome.test = the_groups,
          plot = TRUE, roc.comp = comp_num,
          title = paste0("ROC Curve (VIP>1):", overall_analysis), print = FALSE
        )
      }

      if (!is.null(cv_opt)) {
        if (cv_opt == "loocv") {
          set.seed(123)
          loocv_results2 <- mixOmics::perf(cytokine_splsda2, validation = "loo")
          loocv_error_rate2 <- loocv_results2$error.rate$overall[
            "comp2",
            "max.dist"
          ]
          loocv_acc2 <- 100 * signif(1 - loocv_error_rate2, digits = 2)
          print(paste0(
            current_level, " LOOCV Accuracy (VIP>1): ",
            loocv_acc2, "%"
          ))

          error_rates2 <- loocv_results2$error.rate$overall[, "max.dist"]
          error_df2 <- as.data.frame(error_rates2)
          error_df2$Component <- rownames(error_df2)
          error_df2 <- reshape2::melt(error_df2,
            id.vars = "Component",
            variable.name = "Distance",
            value.name = "ErrorRate"
          )

          a <- ggplot2::ggplot(error_df2, aes(
            x = Component, y = ErrorRate,
            color = Distance, group = 1
          )) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 3) +
            ggplot2::labs(
              title = paste("LOOCV Error Rate (VIP>1):", overall_analysis),
              x = "Number of Components",
              y = "Error Rate"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ggplot2::scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        } else if (cv_opt == "Mfold") {
          set.seed(123)
          fold_results2 <- mixOmics::perf(cytokine_splsda2,
            validation = "Mfold",
            folds = fold_num, nrepeat = 1000
          )
          fold_error_rate2 <- fold_results2$error.rate$overall[
            "comp2",
            "max.dist"
          ]
          fold_acc2 <- 100 * signif(1 - fold_error_rate2, digits = 2)
          print(paste0(
            current_level, " Mfold Accuracy (VIP>1): ",
            fold_acc2, "%"
          ))

          error_rates2 <- fold_results2$error.rate$overall[, "max.dist"]
          error_df2 <- as.data.frame(error_rates2)
          error_df2$Component <- rownames(error_df2)
          error_df2 <- reshape2::melt(error_df2,
            id.vars = "Component",
            variable.name = "Distance",
            value.name = "ErrorRate"
          )

          a <- ggplot2::ggplot(error_df2, aes(
            x = Component, y = ErrorRate,
            color = Distance, group = 1
          )) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 3) +
            ggplot2::labs(
              title = paste("Mfold Error Rate (VIP>1):", overall_analysis),
              x = "Number of Components",
              y = "Error Rate"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ggplot2::scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
      }
      if (conf_mat == TRUE) {
        cat(paste0("Confusion Matrix for PLS-DA Comparison: ", current_level, "\n"))

        # Confusion Matrix for main model
        cm <- caret::confusionMatrix(
          data = as.factor(prediction1[, 3]), # predicted
          reference = as.factor(prediction1[, 1]) # actual
        )
        print(cm$table)
        cat("Accuracy:", signif(cm$overall["Accuracy"], 2), "\n")

        # Check if binary or multi-class
        if (nlevels(as.factor(prediction1[, 1])) == 2) {
          cat("Sensitivity:", signif(cm$byClass["Sensitivity"], 2), "\n")
          cat("Specificity:", signif(cm$byClass["Specificity"], 2), "\n")
        } else {
          cat("\nPer-Class Sensitivity:\n")
          print(signif(cm$byClass[, "Sensitivity"], 2))
          cat("\nPer-Class Specificity:\n")
          print(signif(cm$byClass[, "Specificity"], 2))
          macro_sens <- mean(cm$byClass[, "Sensitivity"], na.rm = TRUE)
          macro_spec <- mean(cm$byClass[, "Specificity"], na.rm = TRUE)
          cat("\nMacro-Averaged Sensitivity:", signif(macro_sens, 2), "\n")
          cat("Macro-Averaged Specificity:", signif(macro_spec, 2), "\n")
        }

        cat(paste0("Confusion Matrix for PLS-DA Comparison with VIP > 1: ", current_level, "\n"))

        # Confusion Matrix for VIP>1 model
        cm_vip <- caret::confusionMatrix(
          data = as.factor(prediction2[, 3]), # predicted
          reference = as.factor(prediction2[, 1]) # actual
        )
        print(cm_vip$table)
        cat("Accuracy:", signif(cm_vip$overall["Accuracy"], 2), "\n")
        if (nlevels(as.factor(prediction2[, 1])) == 2) {
          cat("Sensitivity:", signif(cm_vip$byClass["Sensitivity"], 2), "\n")
          cat("Specificity:", signif(cm_vip$byClass["Specificity"], 2), "\n")
        } else {
          cat("\nPer-Class Sensitivity:\n")
          print(signif(cm_vip$byClass[, "Sensitivity"], 2))
          cat("\nPer-Class Specificity:\n")
          print(signif(cm_vip$byClass[, "Specificity"], 2))
          macro_sens_vip <- mean(cm_vip$byClass[, "Sensitivity"], na.rm = TRUE)
          macro_spec_vip <- mean(cm_vip$byClass[, "Specificity"], na.rm = TRUE)
          cat("\nMacro-Averaged Sensitivity:", signif(macro_sens_vip, 2), "\n")
          cat("Macro-Averaged Specificity:", signif(macro_spec_vip, 2), "\n")
        }
      }
    }
  }
  dev.off()
}
