#########################################################################################################
# Function to generate a pdf file of plots of pls-da, components, and VIP score of components.
# Author: Shubh Saraswat
# Arguments:
#   x.df: a matrix or data frame with groups, stimulation or treatments, and continuous variables
#   colors:  option to have specific colors representing different groups
#   ellipse: option to have ellipse drawn in the classification figure. Default set to false.
#   bg: option to have a prediction background drawn in the classification figure. Default set to false.
#   conf.mat: option to print confusion matrix to show classification from model. Default set to false.
#########################################################################################################

#' Analyze data with Sparse Partial Least Squares Discriminant Analysis (sPLS-DA).
#'
#' @param x.df A matrix or data frame of variables.
#' @param group.col String value of the column name for groups. Default set to NULL
#' @param trt.col String value of the column name for treatments.Default set to NULL
#' @param colors Vector of colors to be set, list of colors to be set equal to the number of groups or treatments. Default set to NULL and will result in random colors.
#' @param title Title of the PDF file to be saved with figures.
#' @param ellipse If ellipse should be drawn in the figures or not. Default set to FALSE.
#' @param bg If prediction background should be drawn in the figures or not. Default set to FALSE.
#' @param conf.mat If confusion matrix of the classifications be printed at the end or not. Default set to FALSE.
#' @param var.num Number of variables to be used in PLS-DA.
#' @param cv.opt Option for "loocv" or "Mfold" cross-validation methods.
#' @param fold.num Number of folds if cv.opt equals "Mfold".
#' @param scale Option to transform data using log2. Default set to NULL.
#' @param comp.num Number of components to be calculated. Default set to 2.
#' @param pch.values A vector of integers specifying the plotting characters to be used in the plot.
#' @param style A character input to have a 3D plot using '3D' or '3d'. Default set to NULL.
#' @description
#' This function conducts Sparse Partial Least Squares Discriminant Analysis based on the groups and treatment inputed. The rest of columns will be assumed as continuous variables to be used in analysis.
#' The function supports log2 transformation of the data and raw data analysis. The function also supports cross-validation methods such as Leave-One-Out Cross-Validation (LOOCV) and Mfold Cross-Validation.
#' @return A PDF file consisting of classification figures, component figures including Variable of importance in projection (VIP) scores, and
#' classification of groups with VIP scores greater than 1.
#' @examples
#' \dontrun{
#' cyt.plsda(x.df, title = "Example PLS-DA Analysis.pdf", bg = TRUE, conf.mat = TRUE, scale = "log2",
#' var.num = 25, cv.opt = "loocv", comp.num = 2, colors = c("black", "purple", "red2"),
#' pch.values = c(16,4,3))
#' }
#' @import mixOmics
#' @export

cyt.plsda <- function(x.df, group.col, trt.col = NULL, colors = NULL, title, ellipse = FALSE, bg = FALSE,
                      conf.mat = FALSE, var.num, cv.opt = NULL, fold.num = 5,
                      scale = NULL, comp.num = 2, pch.values, style = NULL){

  if(!is.null(scale) && scale == "log2"){
    # Log2 transforming cytokines
    x.df <- data.frame(x.df[, c(group.col, trt.col)], log2(x.df[, !(names(x.df) %in% c(group.col, trt.col))]))
    print("Results based on log2 transformation:")
  } else if (is.null(scale)){
    print("Results based on no transformation:")
  }

  # Making the specified columns lowercase
  names(x.df)[names(x.df) %in% c(group.col, trt.col)] <- tolower(names(x.df)[names(x.df) %in% c(group.col, trt.col)])
  group.col <- tolower(group.col)
  if(!is.null(trt.col)) trt.col <- tolower(trt.col)
  # Create a vector for treatment
  levels.vec <- unique(x.df[[trt.col]])

  # Generate a color palette based on the number of groups
  if (is.null(colors)) {
    num_groups <- length(unique(x.df[[group.col]]))
    colors <- rainbow(num_groups)
  }

  pdf(file = title)
  for(i in seq_along(levels.vec)) {
    #i = 1
    current.level <- levels.vec[i]
    condt <- x.df[[trt.col]] == current.level
    Title <- current.level
    theData.df <- x.df[condt, -which(names(x.df) %in% c(group.col, trt.col))]
    theGroups <- x.df[condt, group.col]

    cytokine.splsda <- splsda(theData.df, theGroups, scale = TRUE, ncomp = comp.num, keepX = c(var.num, var.num))

    splsda.predict <- predict(cytokine.splsda, theData.df, dist = "max.dist")
    Prediction1 <- cbind(original = theGroups, splsda.predict$class$max.dist)
    accuracy1 <- sum(Prediction1[,1] == Prediction1[,2]) / length(Prediction1[,1])
    acc1 <- 100 * signif(accuracy1, digits = 2)

    # Creating a shaded predicted background using max.dist
    bg.maxdist <- background.predict(cytokine.splsda, comp.predicted = 1, dist = 'max.dist')

    group_factors <- sort(unique(theGroups))

    # PLS-DA plot using plotIndiv
    plot_args <- list(cytokine.splsda, ind.names = NA, legend = TRUE, col = colors, pch = pch.values, pch.levels = group_factors, title = paste(Title, "With Accuracy:", acc1, "%"))
    if(ellipse) plot_args$ellipse <- TRUE
    if(bg) plot_args$background <- bg.maxdist
    do.call(plotIndiv, plot_args)

    # 3D Plot for initial PLS-DA plot
    if(!is.null(style) && comp.num == 3 && (style == "3D" || style == "3d")){
      cytokine.scores <- cytokine.splsda$variates$X
      plot3D::scatter3D(cytokine.scores[,1], cytokine.scores[,2], cytokine.scores[,3], pch = pch.values, col = colors,
                        xlab = "Component 1", ylab = "Component 2", zlab = "Component 3", main = paste("3D Plot", ":", Title),
                        theta = 20, phi = 30, bty = "g", colkey = FALSE)
    }
    ## Cross-validation methods
    if(!is.null(cv.opt)) {
      if(cv.opt == "loocv"){
        set.seed(123)
        loocv_results <- perf(cytokine.splsda, validation = "loo")
        loocv_error_rate <- loocv_results$error.rate$overall["comp1", "max.dist"]
        loocv_acc <- 1 - loocv_error_rate
        loocv_acc <- 100 * signif(loocv_acc, digits = 2)
        print(paste0(current.level, " LOOCV Accuracy: ", loocv_acc, "%"))

        error_rates <- loocv_results$error.rate$overall[,"max.dist"]
        error_df <- as.data.frame(error_rates)
        error_df$Component <- rownames(error_df)
        error_df <- reshape2::melt(error_df, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

        a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
          geom_line() +
          geom_point(size = 3) +
          labs(title = paste("LOOCV Error Rate", ":", Title),
               x = "Number of Components",
               y = "Error Rate") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      } else if(cv.opt == "Mfold"){
        set.seed(123)
        fold_results <- perf(cytokine.splsda, validation = "Mfold", folds = fold.num, nrepeat = 1000)
        fold_error_rate <- fold_results$error.rate$overall["comp1", "max.dist"]
        fold_acc <- 1 - fold_error_rate
        fold_acc <- 100 * signif(fold_acc, digits = 2)
        print(paste0(current.level, " Mfold Accuracy: ", fold_acc, "%"))

        error_rates <- fold_results$error.rate$overall[,"max.dist"]
        error_df <- as.data.frame(error_rates)
        error_df$Component <- rownames(error_df)
        error_df <- reshape2::melt(error_df, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

        a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
          geom_line() +
          geom_point(size = 3) +
          labs(title = paste("Mfold Error Rate", ":", Title),
               x = "Number of Components",
               y = "Error Rate") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      }
    }
    # Loadings plot
    for(comp in 1:comp.num) {
      plotLoadings(cytokine.splsda, comp = comp, contrib = 'max', method = 'mean', size.name = 1, size.legend = 1, legend.color = colors, title = paste("Component", comp, ":", Title), size.title = 1, legend.title = "Group")
    }

    # VIP scores and plot for PLS-DA with VIP > 1
    all_VIP_scores <- vip(cytokine.splsda)

    for (comp in 1:comp.num) {
      Vscore <- as.data.frame(all_VIP_scores[, comp, drop = FALSE])
      Vscore$metabo <- rownames(Vscore)
      Vscore$comp <- Vscore[,1]
      bar <- Vscore[,c('metabo','comp')]
      bar <- bar[order(bar$comp, decreasing = TRUE), ]

      a <- ggplot(bar, aes(x = metabo, y = comp)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_y_continuous(limits = c(0, max(bar$comp))) +
        geom_hline(yintercept = 1, color = "grey") +
        scale_x_discrete(limits = factor(bar$metabo)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
        labs(x = "", y = "VIP score") +
        ggtitle(paste("Component", comp)) +
        theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'))
      print(a)
    }

    # PLS-DA on VIP > 1
    condtVariable <- all_VIP_scores[, 1] > 1
    KeepX <- sum(condtVariable)
    theData.mat <- theData.df[, condtVariable, drop = FALSE]
    cytokine.splsda2 <- splsda(theData.mat, theGroups, scale = TRUE, ncomp = comp.num, keepX = c(KeepX, KeepX))

    splsda.predict2 <- predict(cytokine.splsda2, theData.mat, dist = "max.dist")
    Prediction2 <- cbind(original = theGroups, splsda.predict2$class$max.dist)
    accuracy2 <- sum(Prediction2[,1] == Prediction2[,2]) / length(Prediction2[,1])
    acc2 <- 100 * signif(accuracy2, digits = 2)

    # PLS-DA plot using plotIndiv for VIP > 1
    plot_args2 <- list(cytokine.splsda2, ind.names = NA, legend = TRUE, col = colors, pch = pch.values, pch.levels = group_factors, title = paste(Title, "(VIP>1)", "With Accuracy:", acc2, "%"))
    if(ellipse) plot_args2$ellipse <- TRUE
    if(bg) plot_args2$background <- bg.maxdist
    do.call(plotIndiv, plot_args2)

    # 3D Plot for VIP > 1
    if(!is.null(style) && comp.num == 3 && (style == "3D" || style == "3d")){
      cytokine.scores2 <- cytokine.splsda2$variates$X
      plot3D::scatter3D(cytokine.scores2[,1], cytokine.scores2[,2], cytokine.scores2[,3], pch = pch.values, col = colors,
                        xlab = "Component 1", ylab = "Component 2", zlab = "Component 3", main = paste("3D Plot", ":", Title, "(VIP>1)"),
                        theta = 20, phi = 30, bty = "g", colkey = FALSE)
    }

    ## Cross-validation methods
    if(!is.null(cv.opt)) {
      if(cv.opt == "loocv"){
        set.seed(123)
        loocv_results2 <- perf(cytokine.splsda2, validation = "loo")
        loocv_error_rate2 <- loocv_results2$error.rate$overall["comp1", "max.dist"]
        loocv_acc2 <- 1 - loocv_error_rate2
        loocv_acc2 <- 100 * signif(loocv_acc2, digits = 2)
        print(paste0(current.level, " LOOCV Accuracy (VIP>1): ", loocv_acc2, "%"))

        error_rates2 <- loocv_results2$error.rate$overall[,"max.dist"]
        error_df2 <- as.data.frame(error_rates2)
        error_df2$Component <- rownames(error_df2)
        error_df2 <- reshape2::melt(error_df2, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

        a <- ggplot(error_df2, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
          geom_line() +
          geom_point(size = 3) +
          labs(title = paste("LOOCV Error Rate (VIP>1)", ":", Title),
               x = "Number of Components",
               y = "Error Rate") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      } else if(cv.opt == "Mfold"){
        set.seed(123)
        fold_results2 <- perf(cytokine.splsda2, validation = "Mfold", folds = fold.num, nrepeat = 1000)
        fold_error_rate2 <- fold_results2$error.rate$overall["comp1", "max.dist"]
        fold_acc2 <- 1 - fold_error_rate2
        fold_acc2 <- 100 * signif(fold_acc2, digits = 2)
        print(paste0(current.level, " Mfold Accuracy (VIP>1): ", fold_acc2, "%"))

        error_rates2 <- fold_results2$error.rate$overall[,"max.dist"]
        error_df2 <- as.data.frame(error_rates2)
        error_df2$Component <- rownames(error_df2)
        error_df2 <- reshape2::melt(error_df2, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

        a <- ggplot(error_df2, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
          geom_line() +
          geom_point(size = 3) +
          labs(title = paste("Mfold Error Rate (VIP>1)", ":", Title),
               x = "Number of Components",
               y = "Error Rate") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          scale_color_manual(values = "red", labels = "max.dist")
        print(a)
      }
    }

    # Loadings plot for VIP > 1
    for(comp in 1:comp.num) {
      plotLoadings(cytokine.splsda2, comp = comp, contrib = 'max', method = 'mean', size.name = 1, size.legend = 1, legend.color = colors, title = paste("Component", comp, ":", Title, "(VIP>1)"), size.title = 1, legend.title = "Group")
    }

    # Confusion matrix printing
    if(conf.mat == TRUE){
      print(paste0(current.level, " Confusion Matrix for PLS-DA Comparison"))
      print(get.confusion_matrix(truth = Prediction1[,1], predicted = Prediction1[,2]))
      print(paste0(current.level, " Confusion Matrix for PLS-DA Comparison with VIP Score > 1"))
      print(get.confusion_matrix(truth = Prediction2[,1], predicted = Prediction2[,2]))
    }
  }
  dev.off()
}

