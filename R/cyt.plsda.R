
#' Analyze data with Sparse Partial Least Squares Discriminant Analysis (sPLS-DA).
#'
#' @param x.df A matrix or data frame containing the variables. Columns not specified by \code{group.col} or \code{trt.col} are assumed to be continuous variables for analysis.
#' @param group.col A string specifying the column name that contains group information. If \code{trt.col} is not provided, it will be used for both grouping and treatment.
#' @param trt.col A string specifying the column name for treatments. Default is \code{NULL}.
#' @param colors A vector of colors for the groups or treatments. If \code{NULL}, a random palette (using \code{rainbow}) is generated based on the number of groups.
#' @param title A string specifying the file name for saving the PDF output.
#' @param ellipse Logical. Whether to draw a 95% confidence ellipse on the figures. Default is \code{FALSE}.
#' @param bg Logical. Whether to draw the prediction background in the figures. Default is \code{FALSE}.
#' @param conf.mat Logical. Whether to print the confusion matrix for the classifications. Default is \code{FALSE}.
#' @param var.num Numeric. The number of variables to be used in the PLS-DA model.
#' @param cv.opt Character. Option for cross-validation method: either "loocv" or "Mfold". Default is \code{NULL}.
#' @param fold.num Numeric. The number of folds to use if \code{cv.opt} is "Mfold". Default is 5.
#' @param scale Character. Option for data transformation; if set to \code{"log2"}, a log2 transformation is applied to the continuous variables. Default is \code{NULL}.
#' @param comp.num Numeric. The number of components to calculate in the sPLS-DA model. Default is 2.
#' @param pch.values A vector of integers specifying the plotting characters to be used in the plots.
#' @param style Character. If set to \code{"3D"} or \code{"3d"} and \code{comp.num} equals 3, a 3D plot is generated using the \code{plot3D} package. Default is \code{NULL}.
#' @param roc Logical. Whether to compute and plot the ROC curve for the model. Default is \code{FALSE}.
#'
#' @description
#' This function conducts Sparse Partial Least Squares Discriminant Analysis (sPLS-DA) on the provided data.
#' It uses the specified \code{group.col} (and optionally \code{trt.col}) to define class labels while assuming the remaining columns
#' contain continuous variables. The function supports a log2 transformation via the \code{scale} parameter and generates a series of plots,
#' including classification plots, scree plots, loadings plots, and VIP score plots. Optionally, ROC curves are produced when \code{roc} is \code{TRUE}.
#' Additionally, cross-validation is supported via LOOCV or Mfold methods. When both \code{group.col} and \code{trt.col} are provided and differ,
#' the function analyzes each treatment level separately.
#'
#' @return A PDF file containing the classification figures, component figures with Variable of Importance in Projection (VIP) scores,
#' and classifications based on VIP scores greater than 1. ROC curves and confusion matrices are also produced if requested.
#'
#' @examples
#' \dontrun{
#' # Example using overall analysis (single factor)
#' cyt.plsda(x.df = my_data, group.col = "Group", title = "PLSDA_overall.pdf",
#'           var.num = 25, pch.values = c(16, 4, 3))
#'
#' # Example using separate group and treatment columns
#' cyt.plsda(x.df = my_data, group.col = "Group", trt.col = "Treatment", title = "PLSDA_byTreatment.pdf",
#'           colors = c("black", "purple", "red2"), ellipse = TRUE, bg = TRUE, conf.mat = TRUE,
#'           var.num = 25, cv.opt = "loocv", comp.num = 2, pch.values = c(16, 4, 3))
#'
#' # Example with ROC curve and 3D plot (comp.num == 3)
#' cyt.plsda(x.df = my_data, group.col = "Group", trt.col = "Treatment", title = "PLSDA_ROC_3D.pdf",
#'           colors = c("black", "purple", "red2"), ellipse = TRUE, bg = TRUE, conf.mat = TRUE,
#'           var.num = 25, cv.opt = "Mfold", fold.num = 5, scale = "log2", comp.num = 3,
#'           pch.values = c(16, 4, 3), style = "3D", roc = TRUE)
#' }
#'
#' @export
#' @import mixOmics
#' @import ggplot2
#' @import plot3D
#' @import reshape2
#' @import caret

cyt.plsda <- function(x.df, group.col = NULL, trt.col = NULL, colors = NULL, title,
                      ellipse = FALSE, bg = FALSE, conf.mat = FALSE, var.num, cv.opt = NULL,
                      fold.num = 5, scale = NULL, comp.num = 2, pch.values, style = NULL, roc = FALSE){

  # If one factor is missing, use the provided column for both grouping and treatment.
  if(is.null(group.col) && !is.null(trt.col)){
    message("No group column provided; using the treatment column as the grouping variable.")
    group.col <- trt.col
  }
  if(is.null(trt.col) && !is.null(group.col)){
    message("No treatment column provided; using the group column as the treatment variable.")
    trt.col <- group.col
  }
  if(is.null(group.col) && is.null(trt.col)){
    stop("At least one factor column must be provided.")
  }

  # Optionally apply log2 transformation
  if(!is.null(scale) && scale == "log2"){
    x.df <- data.frame(x.df[, c(group.col, trt.col)],
                       log2(x.df[, !(names(x.df) %in% c(group.col, trt.col))]))
    print("Results based on log2 transformation:")
  } else if (is.null(scale)){
    print("Results based on no transformation:")
  }

  # Convert factor column names to lowercase for consistency
  names(x.df)[names(x.df) %in% c(group.col, trt.col)] <-
    tolower(names(x.df)[names(x.df) %in% c(group.col, trt.col)])
  group.col <- tolower(group.col)
  trt.col <- tolower(trt.col)

  # Generate a color palette if not provided (based on the grouping variable levels in the entire dataset)
  if (is.null(colors)) {
    num_groups <- length(unique(x.df[[group.col]]))
    colors <- rainbow(num_groups)
  }

  pdf(file = title, width = 8.5, height = 8)

  # Case 1: Only one factor provided (both columns are the same)
  if(group.col == trt.col){
    Title <- "Overall Analysis"

    # Remove the factor column from predictors and keep only numeric columns
    theData.df <- x.df[, !(names(x.df) %in% c(group.col))]
    theData.df <- theData.df[, sapply(theData.df, is.numeric)]

    theGroups <- as.vector(x.df[[group.col]])
    if(length(unique(theGroups)) < 2){
      stop("The grouping variable must have at least two levels for PLS-DA. Please provide an appropriate grouping column.")
    }

    cytokine.splsda <- mixOmics::splsda(theData.df, theGroups, scale = TRUE, ncomp = comp.num,
                              keepX = rep(var.num, comp.num))

    splsda.predict <- predict(cytokine.splsda, theData.df, dist = "max.dist")
    Prediction1 <- cbind(original = theGroups, splsda.predict$class$max.dist)
    accuracy1 <- sum(Prediction1[,1] == Prediction1[,3]) / length(Prediction1[,1])
    acc1 <- 100 * signif(accuracy1, digits = 2)

    bg.maxdist <- background.predict(cytokine.splsda, comp.predicted = 2, dist = 'max.dist')
    group_factors <- sort(unique(theGroups))

    plot_args <- list(cytokine.splsda, ind.names = NA, legend = TRUE, col = colors,
                      pch = pch.values, pch.levels = group_factors,
                      title = paste(Title, "With Accuracy:", acc1, "%"), legend.title = group.col)
    if(ellipse) plot_args$ellipse <- TRUE
    if(bg) plot_args$background <- bg.maxdist
    do.call(plotIndiv, plot_args)

    if(!is.null(style) && comp.num == 3 && (tolower(style) == "3d")){
      cytokine.scores <- cytokine.splsda$variates$X
      plot3D::scatter3D(cytokine.scores[,1], cytokine.scores[,2], cytokine.scores[,3],
                        pch = pch.values, col = colors,
                        xlab = "Component 1", ylab = "Component 2", zlab = "Component 3",
                        main = paste("3D Plot:", Title),
                        theta = 20, phi = 30, bty = "g", colkey = FALSE)
    }

    # If roc = TRUE, compute and plot ROC curve for the overall model
    if(roc){
      # Using auroc() from mixOmics: this will use the overall model, training data, and outcome vector.
      roc_obj <- auroc(object = cytokine.splsda, newdata = theData.df, outcome.test = theGroups,
                       plot = TRUE, roc.comp = comp.num,
                       title = paste0("ROC Curve:", Title), print = FALSE)
    }

    ## Cross-validation methods
    if(!is.null(cv.opt)) {
      if(cv.opt == "loocv"){
        set.seed(123)
        loocv_results <- perf(cytokine.splsda, validation = "loo")
        loocv_error_rate <- loocv_results$error.rate$overall["comp2", "max.dist"]
        loocv_acc <- 100 * signif(1 - loocv_error_rate, digits = 2)
        print(paste0("LOOCV Accuracy: ", loocv_acc, "%"))

        error_rates <- loocv_results$error.rate$overall[,"max.dist"]
        error_df <- as.data.frame(error_rates)
        error_df$Component <- rownames(error_df)
        error_df <- reshape2::melt(error_df, id.vars = "Component",
                                   variable.name = "Distance", value.name = "ErrorRate")

        a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
          geom_line() +
          geom_point(size = 3) +
          labs(title = paste("LOOCV Error Rate:", Title),
               x = "Number of Components",
               y = "Error Rate") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      } else if(cv.opt == "Mfold"){
        set.seed(123)
        fold_results <- perf(cytokine.splsda, validation = "Mfold", folds = fold.num, nrepeat = 1000)
        fold_error_rate <- fold_results$error.rate$overall["comp2", "max.dist"]
        fold_acc <- 100 * signif(1 - fold_error_rate, digits = 2)
        print(paste0("Mfold Accuracy: ", fold_acc, "%"))

        error_rates <- fold_results$error.rate$overall[,"max.dist"]
        error_df <- as.data.frame(error_rates)
        error_df$Component <- rownames(error_df)
        error_df <- reshape2::melt(error_df, id.vars = "Component",
                                   variable.name = "Distance", value.name = "ErrorRate")

        a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
          geom_line() +
          geom_point(size = 3) +
          labs(title = paste("Mfold Error Rate:", Title),
               x = "Number of Components",
               y = "Error Rate") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      }
    }

    # Loadings plot for each component
    for(comp in 1:comp.num) {
      plotLoadings(cytokine.splsda, comp = comp, contrib = 'max', method = 'mean',
                   size.name = 1, size.legend = 1, legend.color = colors,
                   title = paste("Component", comp, ":", Title), size.title = 1, legend.title = group.col)
    }

    # VIP scores and plot for PLS-DA with VIP > 1
    all_VIP_scores <- vip(cytokine.splsda)
    for (comp in 1:comp.num) {
      Vscore <- as.data.frame(all_VIP_scores[, comp, drop = FALSE])
      Vscore$metabo <- rownames(Vscore)
      Vscore$comp <- Vscore[,1]
      bar <- Vscore[, c('metabo','comp')]
      bar <- bar[order(bar$comp, decreasing = TRUE), ]

      a <- ggplot(bar, aes(x = metabo, y = comp)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_y_continuous(limits = c(0, max(bar$comp))) +
        geom_hline(yintercept = 1, color = "grey") +
        scale_x_discrete(limits = factor(bar$metabo)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
        labs(x = "", y = "VIP score") +
        ggtitle(paste("Component", comp)) +
        theme(panel.grid = element_blank(),
              panel.background = element_rect(color = 'black', fill = 'transparent'))
      print(a)
    }

    # PLS-DA on VIP > 1: Subset predictors with VIP > 1
    condtVariable <- all_VIP_scores[, 1] > 1
    KeepX <- sum(condtVariable)
    theData.mat <- theData.df[, condtVariable, drop = FALSE]
    cytokine.splsda2 <- mixOmics::splsda(theData.mat, theGroups, scale = TRUE, ncomp = comp.num,
                               keepX = rep(KeepX, comp.num))

    splsda.predict2 <- predict(cytokine.splsda2, theData.mat, dist = "max.dist")
    Prediction2 <- cbind(original = theGroups, splsda.predict2$class$max.dist)
    accuracy2 <- sum(Prediction2[,1] == Prediction2[,3]) / length(Prediction2[,1])
    acc2 <- 100 * signif(accuracy2, digits = 2)

    bg.maxdist2 <- background.predict(cytokine.splsda2, comp.predicted = 2, dist = 'max.dist')

    plot_args2 <- list(cytokine.splsda2, ind.names = NA, legend = TRUE, col = colors,
                       pch = pch.values, pch.levels = group_factors,
                       title = paste(Title, "(VIP>1)", "With Accuracy:", acc2, "%"), legend.title = group.col)
    if(ellipse) plot_args2$ellipse <- TRUE
    if(bg) plot_args2$background <- bg.maxdist2
    do.call(plotIndiv, plot_args2)

    if(!is.null(style) && comp.num == 3 && (tolower(style) == "3d")){
      cytokine.scores2 <- cytokine.splsda2$variates$X
      plot3D::scatter3D(cytokine.scores2[,1], cytokine.scores2[,2], cytokine.scores2[,3],
                        pch = pch.values, col = colors,
                        xlab = "Component 1", ylab = "Component 2", zlab = "Component 3",
                        main = paste("3D Plot:", Title, "(VIP>1)"),
                        theta = 20, phi = 30, bty = "g", colkey = FALSE)
    }

    if(!is.null(cv.opt)) {
      if(cv.opt == "loocv"){
        set.seed(123)
        loocv_results2 <- perf(cytokine.splsda2, validation = "loo")
        loocv_error_rate2 <- loocv_results2$error.rate$overall["comp2", "max.dist"]
        loocv_acc2 <- 100 * signif(1 - loocv_error_rate2, digits = 2)
        print(paste0("LOOCV Accuracy (VIP>1): ", loocv_acc2, "%"))

        error_rates2 <- loocv_results2$error.rate$overall[,"max.dist"]
        error_df2 <- as.data.frame(error_rates2)
        error_df2$Component <- rownames(error_df2)
        error_df2 <- reshape2::melt(error_df2, id.vars = "Component",
                                    variable.name = "Distance", value.name = "ErrorRate")

        a <- ggplot(error_df2, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
          geom_line() +
          geom_point(size = 3) +
          labs(title = paste("LOOCV Error Rate (VIP>1):", Title),
               x = "Number of Components",
               y = "Error Rate") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      } else if(cv.opt == "Mfold"){
        set.seed(123)
        fold_results2 <- perf(cytokine.splsda2, validation = "Mfold", folds = fold.num, nrepeat = 1000)
        fold_error_rate2 <- fold_results2$error.rate$overall["comp2", "max.dist"]
        fold_acc2 <- 100 * signif(1 - fold_error_rate2, digits = 2)
        print(paste0("Mfold Accuracy (VIP>1): ", fold_acc2, "%"))

        error_rates2 <- fold_results2$error.rate$overall[,"max.dist"]
        error_df2 <- as.data.frame(error_rates2)
        error_df2$Component <- rownames(error_df2)
        error_df2 <- reshape2::melt(error_df2, id.vars = "Component",
                                    variable.name = "Distance", value.name = "ErrorRate")

        a <- ggplot(error_df2, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
          geom_line() +
          geom_point(size = 3) +
          labs(title = paste("Mfold Error Rate (VIP>1):", Title),
               x = "Number of Components",
               y = "Error Rate") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      }
    }

    if(conf.mat == TRUE){
      print("Confusion Matrix for PLS-DA Comparison")
      print(get.confusion_matrix(truth = Prediction1[,1], predicted = Prediction1[,2]))
      print("Confusion Matrix for PLS-DA Comparison with VIP Score > 1")
      print(get.confusion_matrix(truth = Prediction2[,1], predicted = Prediction2[,2]))
    }
    # If roc = TRUE, compute and plot ROC curve for the overall model of VIP > 1
    if(roc){
      roc_obj2 <- auroc(object = cytokine.splsda2, newdata = theData.mat, outcome.test = theGroups, plot = TRUE,
                        roc.comp = comp.num, title = paste0("ROC Curve (VIP>1):", Title), print = FALSE)
    }


  } else {
    # Case 2: Both group and treatment columns are provided and they differ.
    levels.vec <- unique(x.df[[trt.col]])
    for(i in seq_along(levels.vec)) {
      current.level <- levels.vec[i]
      Title <- current.level
      condt <- x.df[[trt.col]] == current.level

      theData.df <- x.df[condt, -which(names(x.df) %in% c(group.col, trt.col))]
      theData.df <- theData.df[, sapply(theData.df, is.numeric)]
      theGroups <- as.vector(x.df[condt, group.col])

      if(length(unique(theGroups)) < 2){
        stop("The grouping variable must have at least two levels for PLS-DA. Please provide an appropriate grouping column.")
      }

      cytokine.splsda <- mixOmics::splsda(theData.df, theGroups, scale = TRUE, ncomp = comp.num,
                                keepX = rep(var.num, comp.num))

      splsda.predict <- predict(cytokine.splsda, theData.df, dist = "max.dist")
      Prediction1 <- cbind(original = theGroups, splsda.predict$class$max.dist)
      accuracy1 <- sum(Prediction1[,1] == Prediction1[,3]) / length(Prediction1[,1])
      acc1 <- 100 * signif(accuracy1, digits = 2)

      bg.maxdist <- background.predict(cytokine.splsda, comp.predicted = 2, dist = 'max.dist')
      group_factors <- sort(unique(theGroups))

      plot_args <- list(cytokine.splsda, ind.names = NA, legend = TRUE, col = colors,
                        pch = pch.values, pch.levels = group_factors,
                        title = paste(Title, "With Accuracy:", acc1, "%"), legend.title = group.col)
      if(ellipse) plot_args$ellipse <- TRUE
      if(bg) plot_args$background <- bg.maxdist
      do.call(plotIndiv, plot_args)

      if(!is.null(style) && comp.num == 3 && (tolower(style) == "3d")){
        cytokine.scores <- cytokine.splsda$variates$X
        plot3D::scatter3D(cytokine.scores[,1], cytokine.scores[,2], cytokine.scores[,3],
                          pch = pch.values, col = colors,
                          xlab = "Component 1", ylab = "Component 2", zlab = "Component 3",
                          main = paste("3D Plot:", Title),
                          theta = 20, phi = 30, bty = "g", colkey = FALSE)
      }

      # If roc = TRUE, compute and plot ROC curve for the overall model
      if(roc){
        # Using auroc() from mixOmics: this will use the overall model, training data, and outcome vector.
        roc_obj <- auroc(object = cytokine.splsda, newdata = theData.df, outcome.test = theGroups,
                         plot = TRUE, roc.comp = comp.num,
                         title = paste0("ROC Curve:", Title), print = FALSE)
      }

      if(!is.null(cv.opt)) {
        if(cv.opt == "loocv"){
          set.seed(123)
          loocv_results <- perf(cytokine.splsda, validation = "loo")
          loocv_error_rate <- loocv_results$error.rate$overall["comp2", "max.dist"]
          loocv_acc <- 100 * signif(1 - loocv_error_rate, digits = 2)
          print(paste0(current.level, " LOOCV Accuracy: ", loocv_acc, "%"))

          error_rates <- loocv_results$error.rate$overall[,"max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df, id.vars = "Component",
                                     variable.name = "Distance", value.name = "ErrorRate")

          a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("LOOCV Error Rate:", Title),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        } else if(cv.opt == "Mfold"){
          set.seed(123)
          fold_results <- perf(cytokine.splsda, validation = "Mfold", folds = fold.num, nrepeat = 1000)
          fold_error_rate <- fold_results$error.rate$overall["comp2", "max.dist"]
          fold_acc <- 100 * signif(1 - fold_error_rate, digits = 2)
          print(paste0(current.level, " Mfold Accuracy: ", fold_acc, "%"))

          error_rates <- fold_results$error.rate$overall[,"max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df, id.vars = "Component",
                                     variable.name = "Distance", value.name = "ErrorRate")

          a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("Mfold Error Rate:", Title),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
      }


      for(comp in 1:comp.num) {
        plotLoadings(cytokine.splsda, comp = comp, contrib = 'max', method = 'mean',
                     size.name = 1, size.legend = 1, legend.color = colors,
                     title = paste("Component", comp, ":", Title), size.title = 1, legend.title = group.col)
      }

      all_VIP_scores <- vip(cytokine.splsda)
      for (comp in 1:comp.num) {
        Vscore <- as.data.frame(all_VIP_scores[, comp, drop = FALSE])
        Vscore$metabo <- rownames(Vscore)
        Vscore$comp <- Vscore[,1]
        bar <- Vscore[, c('metabo','comp')]
        bar <- bar[order(bar$comp, decreasing = TRUE), ]

        a <- ggplot(bar, aes(x = metabo, y = comp)) +
          geom_bar(stat = "identity", position = "dodge") +
          scale_y_continuous(limits = c(0, max(bar$comp))) +
          geom_hline(yintercept = 1, color = "grey") +
          scale_x_discrete(limits = factor(bar$metabo)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
          labs(x = "", y = "VIP score") +
          ggtitle(paste("Component", comp)) +
          theme(panel.grid = element_blank(),
                panel.background = element_rect(color = 'black', fill = 'transparent'))
        print(a)
      }

      condtVariable <- all_VIP_scores[, 1] > 1
      KeepX <- sum(condtVariable)
      theData.mat <- theData.df[, condtVariable, drop = FALSE]
      cytokine.splsda2 <- mixOmics::splsda(theData.mat, theGroups, scale = TRUE, ncomp = comp.num,
                                 keepX = rep(KeepX, comp.num))

      splsda.predict2 <- predict(cytokine.splsda2, theData.mat, dist = "max.dist")
      Prediction2 <- cbind(original = theGroups, splsda.predict2$class$max.dist)
      accuracy2 <- sum(Prediction2[,1] == Prediction2[,3]) / length(Prediction2[,1])
      acc2 <- 100 * signif(accuracy2, digits = 2)


      bg.maxdist2 <- background.predict(cytokine.splsda2, comp.predicted = 2, dist = 'max.dist')

      plot_args2 <- list(cytokine.splsda2, ind.names = NA, legend = TRUE, col = colors,
                         pch = pch.values, pch.levels = group_factors,
                         title = paste(Title, "(VIP>1)", "With Accuracy:", acc2, "%"), legend.title = group.col)
      if(ellipse) plot_args2$ellipse <- TRUE
      if(bg) plot_args2$background <- bg.maxdist
      do.call(plotIndiv, plot_args2)

      if(!is.null(style) && comp.num == 3 && (tolower(style) == "3d")){
        cytokine.scores2 <- cytokine.splsda2$variates$X
        plot3D::scatter3D(cytokine.scores2[,1], cytokine.scores2[,2], cytokine.scores2[,3],
                          pch = pch.values, col = colors,
                          xlab = "Component 1", ylab = "Component 2", zlab = "Component 3",
                          main = paste("3D Plot:", Title, "(VIP>1)"),
                          theta = 20, phi = 30, bty = "g", colkey = FALSE)
      }

      if(roc){
        roc_obj2 <- auroc(object = cytokine.splsda2, newdata = theData.mat, outcome.test = theGroups, plot = TRUE,
                          roc.comp = comp.num, title = paste0("ROC Curve (VIP>1):", Title), print = FALSE)
      }

      if(!is.null(cv.opt)) {
        if(cv.opt == "loocv"){
          set.seed(123)
          loocv_results2 <- perf(cytokine.splsda2, validation = "loo")
          loocv_error_rate2 <- loocv_results2$error.rate$overall["comp2", "max.dist"]
          loocv_acc2 <- 100 * signif(1 - loocv_error_rate2, digits = 2)
          print(paste0(current.level, " LOOCV Accuracy (VIP>1): ", loocv_acc2, "%"))

          error_rates2 <- loocv_results2$error.rate$overall[,"max.dist"]
          error_df2 <- as.data.frame(error_rates2)
          error_df2$Component <- rownames(error_df2)
          error_df2 <- reshape2::melt(error_df2, id.vars = "Component",
                                      variable.name = "Distance", value.name = "ErrorRate")

          a <- ggplot(error_df2, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("LOOCV Error Rate (VIP>1):", Title),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        } else if(cv.opt == "Mfold"){
          set.seed(123)
          fold_results2 <- perf(cytokine.splsda2, validation = "Mfold", folds = fold.num, nrepeat = 1000)
          fold_error_rate2 <- fold_results2$error.rate$overall["comp2", "max.dist"]
          fold_acc2 <- 100 * signif(1 - fold_error_rate2, digits = 2)
          print(paste0(current.level, " Mfold Accuracy (VIP>1): ", fold_acc2, "%"))

          error_rates2 <- fold_results2$error.rate$overall[,"max.dist"]
          error_df2 <- as.data.frame(error_rates2)
          error_df2$Component <- rownames(error_df2)
          error_df2 <- reshape2::melt(error_df2, id.vars = "Component",
                                      variable.name = "Distance", value.name = "ErrorRate")

          a <- ggplot(error_df2, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("Mfold Error Rate (VIP>1):", Title),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
      }

      if(conf.mat == TRUE){
        print(paste0(current.level, " Confusion Matrix for PLS-DA Comparison"))
        #print(get.confusion_matrix(truth = Prediction1[,1], predicted = Prediction1[,3]))
        cm <- caret::confusionMatrix(data = as.factor(Prediction1[,3]), reference = as.factor(Prediction1[,1]))
        print(cm$table)
        print(signif(cm$overall["Accuracy"]),2)
        print(signif(cm$byClass["Sensitivity"]),2)
        print(signif(cm$byClass["Specificity"]),2)
        print(paste0(current.level, " Confusion Matrix for PLS-DA Comparison with VIP Score > 1"))
        # print(get.confusion_matrix(truth = Prediction2[,1], predicted = Prediction2[,3]))
        cm_vip <- caret::confusionMatrix(data = as.factor(Prediction2[,3]), reference = as.factor(Prediction2[,1]))
        print(cm_vip$table)
        print(signif(cm_vip$overall["Accuracy"]),2)
        print(signif(cm_vip$byClass["Sensitivity"]),2)
        print(signif(cm_vip$byClass["Specificity"]),2)
      }
    }
  }

  dev.off()
}

