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
#' This function takes a data frame assuming that the first two columns consists
#' of groups and stimulation/treatment of patients and using the rest of the columns as
#' components for the PLS-DA analysis.
#'
#' @param x.df A matrix or data frame of variables.
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
#' This function takes a matrix or data frame of variables where the first two columns must consist of "Group"
#' and "Stimulation" or "Group" and "Treatment". The rest of columns will be cytokines that will be analyzed.
#' The function also uses log2 transformation of the data, so the data that is inputed, must be raw data without
#' any transformation.
#'
#' @return A PDF file consisting of classification figures, component figures including Variable of importance in projection (VIP) scores, and
#' classification of groups with VIP scores greater than 1.
#' @examples
#' \dontrun{
#' cyt.plsda(x.df, title = "Example PLS-DA Analysis.pdf", bg = TRUE, conf.mat = TRUE, scale = "log2",
#' var.num = 25, cv.opt = "loocv", comp.num = 2, colors = c("black", "purple", "red2"),
#' pch.values = c(16,4,3))
#' }
#' @export

cyt.plsda = function(data.df, colors = NULL, title, ellipse = FALSE, bg = FALSE,
                     conf.mat = FALSE, var.num, cv.opt = NULL, fold.num = 5,
                     scale = NULL, comp.num = 2, pch.values, style = NULL){
  if(!is.null(scale) && scale == "log2"){
    # Log2 transforming cytokines
    data.df = data.frame(data.df[,c(1:2)], log2(data.df[, -c(1:2)]))
    print("Results based on log2 transformation:")
  }
  else if (is.null(scale)){
    print("Results based on no transformation:")
  }

  # Making the first two columns to be lowercase
  names(data.df)[1:2] <- tolower(names(data.df)[1:2])
  # Creating a table to have two separate vectors for group and stimulation
  a = table( data.df[, c(1,2)] ); a

  if("treatment" %in% names(data.df)[1:2]){
    Treatment.vec = dimnames(a)$treatment; Treatment.vec
  }else{
    Stimulation.vec = dimnames(a)$stimulation; Stimulation.vec
  }
  Group.vec = dimnames(a)$group; Group.vec

  # Generate a color palette based on the number of groups
  if (is.null(colors)) {
    num_groups <- length(unique(data.df$group))
    colors <- rainbow(num_groups)
  }

  if("treatment" %in% names(data.df)[1:2]){
    pdf( file= title)
    for(i in 1:length(Treatment.vec) ) {
      #i = 1
      theTrt = Treatment.vec[i]
      condt  = data.df[, "treatment"] == theTrt
      Title  = theTrt
      theData.df = data.df[condt,-c(1:2)]
      theGroups  = data.df[condt, "group"]

      cytokine.splsda =
        splsda(theData.df, theGroups, scale=TRUE, ncomp=comp.num, keepX = c(var.num, var.num))

      splsda.predict = predict( cytokine.splsda, theData.df, dist="max.dist")
      Prediction1 = cbind( original = theGroups, splsda.predict$class$max.dist )
      accuracy1 = sum(Prediction1[,1] == Prediction1[,2])/length(Prediction1[,1]); accuracy1
      acc1 = 100*signif(accuracy1, digits = 2)

      # Creating a shaded predicted background using max.dist
      bg.maxdist <-  background.predict(cytokine.splsda,
                                        comp.predicted = 1,
                                        dist = 'max.dist')

      group_factors = sort(unique(theGroups))

      # Conditions to have either ellipses or prediction background or both or neither on graphs.
      if(ellipse == TRUE & bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=TRUE, title= paste (Title, "With Accuracy:",acc1,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=TRUE, title= paste (Title, "With Accuracy:",acc1,"%"))
      }
      else if(bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=FALSE, title= paste (Title, "With Accuracy:",acc1,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == FALSE & bg == FALSE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=FALSE, title= paste (Title, "With Accuracy:",acc1,"%"))
      }
      if(!is.null(style) && comp.num == 3 && (style == "3D" || style == "3d")){
        # 3D scatter plot using plot3D package
        cytokine.scores <- cytokine.splsda$variates$X
        plot3D::scatter3D(cytokine.scores[,1], cytokine.scores[,2], cytokine.scores[,3], pch=pch.values, col=colors,
                          xlab="Component 1", ylab="Component 2", zlab="Component 3", main=paste("3D Plot", ":", Title),
                          theta=20, phi=30, bty="g", colkey=FALSE)
      }
      else if(!is.null(style)){
        stop("Please enter a valid style for 3D plot:'3d' or '3D' or enter valid number of components.")
      }
      ## Cross-validation methods
      if(!is.null(cv.opt)) {
        # LOOCV Different method
        if(cv.opt == "loocv"){
          set.seed(123) # For reproducibility
          loocv_results = perf(cytokine.splsda, validation = "loo")
          loocv_error_rate = loocv_results$error.rate$overall["comp1", "max.dist"]
          loocv_acc = 1 - loocv_error_rate; loocv_acc
          loocv_acc = 100*signif(loocv_acc, digits = 2)
          print(paste0(theTrt, " ",unique(theGroups)[1]," vs ",unique(theGroups)[2]," LOOCV Accuracy: ", loocv_acc))

          # Extracting data for plotting
          error_rates <- loocv_results$error.rate$overall[,"max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

          # Plotting with ggplot2
          a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("LOOCV Error Rate", ":", theTrt),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
        #Mfold
        else if(cv.opt == "Mfold"){
          set.seed(123) # For reproducibility
          fold_results = perf(cytokine.splsda, validation = "Mfold", folds = fold.num, nrepeat = 1000)
          fold_error_rate = fold_results$error.rate$overall["comp1", "max.dist"]
          fold_acc = 1 - fold_error_rate; fold_acc
          fold_acc = 100*signif(fold_acc, digits = 2)
          print(paste0(theTrt, " ",unique(theGroups)[1]," vs ",unique(theGroups)[2]," Mfold Accuracy: ", fold_acc))

          # Extracting data for plotting
          error_rates <- fold_results$error.rate$overall[,"max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

          # Plotting with ggplot2
          a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("Mfold Error Rate", ":", theTrt),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
      }

      num_rows <- nrow(cytokine.splsda$loadings$X)
      loading.vip.df <- data.frame(matrix(ncol = 0, nrow = num_rows))

      # Loop over the number of components specified by 'comp.num'
      for (comp in 1:comp.num) {
        # Plot loadings for the current component
        plotLoadings(cytokine.splsda, comp=comp, contrib='max', method='mean', size.name=1,
                     size.legend=1, legend.color=colors, title=paste("Component", comp, ":", Title), size.title=1,
                     legend.title = "Group")

        # Add columns to the data frame for the loadings and VIP for the current component
        loading.vip.df[[paste("loading.comp", comp, sep = "")]] <- cytokine.splsda$loadings$X[, comp]
        loading.vip.df[[paste("vip.comp", comp, sep = "")]] <- vip(cytokine.splsda)
      }

      # Update the column names for loadings and VIP for all components
      colnames(loading.vip.df) <- unlist(lapply(1:comp.num, function(comp) {
        c(paste("loading.comp", comp, sep = ""), paste("vip.comp", comp, sep = ""))
      }))


      # Calculate VIP scores once
      all_VIP_scores <- vip(cytokine.splsda)

      # Use all_VIP_scores for plotting
      for (comp in 1:comp.num) {
        Vscore <- as.data.frame(all_VIP_scores[, comp, drop = FALSE])
        Vscore$metabo <- rownames(Vscore)
        Vscore$comp <- Vscore[,1]  # assuming the first column contains the VIP scores
        bar <- Vscore[,c('metabo','comp')]
        rownames(bar) <- 1:nrow(bar)
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

        print(a)  # This prints the plot for each component
      }

      # Conduct further analysis based on a specific component's VIP score
      # For example, use the VIP scores from the first component
      Vscore <- as.data.frame(all_VIP_scores[, 1, drop = FALSE])
      Vscore$metabo <- rownames(Vscore)
      Vscore$comp <- Vscore[,1]
      condtVariable <- Vscore$comp > 1
      KeepX <- sum(condtVariable)
      theData.mat <- theData.df[, condtVariable, drop = FALSE]
      cytokine.splsda2 <- splsda(theData.mat, theGroups, scale=TRUE, ncomp=comp.num, keepX = c(KeepX, KeepX))


      splsda.predict = predict( cytokine.splsda2, theData.mat, dist="max.dist")
      Prediction2 = cbind( original = theGroups, splsda.predict$class$max.dist )
      accuracy2 = sum(Prediction2[,1] == Prediction2[,2])/length(Prediction2[,1]); accuracy2
      acc2 = 100*signif(accuracy2, digits = 2)
      # Creating a shaded predicted background using max.dist
      bg.maxdist <-  background.predict(cytokine.splsda2,
                                        comp.predicted = 1,
                                        dist = 'max.dist')


      # Conditions to have either ellipses or prediction background or both or neither on graphs.
      if(ellipse == TRUE & bg == TRUE){
        plotIndiv(cytokine.splsda2, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=TRUE, title= paste (Title, "(VIP>1)", "With Accuracy:",acc2,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == TRUE) {
        plotIndiv(cytokine.splsda2, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=TRUE, title=paste(Title, "(VIP>1)","With Accuracy", acc2, "%" ))
      }
      else if(bg == TRUE){
        plotIndiv(cytokine.splsda2, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=FALSE, title=paste(Title, "(VIP>1)","With Accuracy", acc2, "%" ),
                  background = bg.maxdist)
      }
      else if(ellipse == FALSE & bg == FALSE){
        plotIndiv(cytokine.splsda2, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=FALSE, title= paste (Title, "(VIP>1)", "With Accuracy:",acc2,"%"))
      }
      if(!is.null(style) && comp.num == 3 && (style == "3D" || style == "3d")){
        # 3D scatter plot using plot3D package
        cytokine.scores2 <- cytokine.splsda2$variates$X
        plot3D::scatter3D(cytokine.scores2[,1], cytokine.scores2[,2], cytokine.scores2[,3], pch=pch.values, col=colors,
                          xlab="Component 1", ylab="Component 2", zlab="Component 3", main=paste("3D Plot", ":", Title, "(VIP>1)"),
                          theta=20, phi=30, bty="g", colkey=FALSE)
      }
      else if(!is.null(style)){
        stop("Please enter a valid style for 3D plot:'3d' or '3D' or enter valid number of components.")
      }
      ## Cross-validation methods
      if(!is.null(cv.opt)){
        # LOOCV Different method
        if(cv.opt == "loocv"){
          set.seed(123) # For reproducibility
          loocv_results2 = perf(cytokine.splsda2, validation = "loo")
          loocv_error_rate2 = loocv_results2$error.rate$overall["comp1", "max.dist"]
          loocv_acc2 = 1 - loocv_error_rate2; loocv_acc2
          loocv_acc2 = 100*signif(loocv_acc2, digits = 2)
          print(paste0(theTrt," ",unique(theGroups)[1],"vs",unique(theGroups)[2]," LOOCV Accuracy (VIP>1) Cytokines: ", loocv_acc2))

          # Extracting data for plotting
          error_rates <- loocv_results2$error.rate$overall[,"max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

          # Plotting with ggplot2
          a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("LOOCV Error Rate VIP > 1", ":", theTrt),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
        #Mfold
        else if(cv.opt == "Mfold"){
          set.seed(123) # For reproducibility
          fold_results2 = perf(cytokine.splsda2, validation = "Mfold", folds = fold.num, nrepeat = 1000)
          fold_error_rate2 = fold_results2$error.rate$overall["comp1", "max.dist"]
          fold_acc2 = 1 - fold_error_rate2; fold_acc2
          fold_acc2 = 100*signif(fold_acc2, digits = 2)
          print(paste0(theTrt," ",unique(theGroups)[1],"vs",unique(theGroups)[2]," Mfold Accuracy (VIP>1) Cytokines: ", fold_acc2))

          # Extracting data for plotting
          error_rates <- fold_results2$error.rate$overall[,"max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

          # Plotting with ggplot2
          a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("Mfold Error Rate VIP > 1", ":", theTrt),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
      }

      num_rows2 <- nrow(cytokine.splsda2$loadings$X)
      loading.vip.df2 <- data.frame(matrix(ncol = 0, nrow = num_rows2))

      # Loop over the number of components specified by 'comp.num'
      for (comp in 1:comp.num) {
        # Plot loadings for the current component
        plotLoadings(cytokine.splsda2, comp=comp, contrib='max', method='mean', size.name=1,
                     size.legend=1, legend.color=colors, title=paste("Component", comp, ":", Title), size.title=1,
                     legend.title = "Group")

        # Add columns to the data frame for the loadings and VIP for the current component
        loading.vip.df2[[paste("loading.comp", comp, sep = "")]] <- cytokine.splsda2$loadings$X[, comp]
        loading.vip.df2[[paste("vip.comp", comp, sep = "")]] <- vip(cytokine.splsda2)
      }

      # Update the column names for loadings and VIP for all components
      colnames(loading.vip.df2) <- unlist(lapply(1:comp.num, function(comp) {
        c(paste("loading.comp", comp, sep = ""), paste("vip.comp", comp, sep = ""))
      }))

      # Prints confusion matrix
      if(conf.mat == TRUE){
        print(paste0(theTrt," ",unique(theGroups)[1]," vs ",unique(theGroups)[2]," Confusion Matrix for PLS-DA Comparison"))
        print(get.confusion_matrix(truth = Prediction1[,1], predicted = Prediction1[,2])) # Confusion matrix for all variables in model
        print(paste0(theTrt," ",unique(theGroups)[1]," vs ",unique(theGroups)[2]," Confusion Matrix for PLS-DA Comparison with VIP Score > 1"))
        print(get.confusion_matrix(truth = Prediction2[,1], predicted = Prediction2[,2])) # Confusion matrix for all variables with VIP Score > 1
      }
    }
    dev.off()
  }else{
    pdf( file= title)
    for(i in 1:length(Stimulation.vec) ) {
      theTrts = Stimulation.vec[i]
      condt  = data.df[, "stimulation"] == theTrts
      Title  = theTrts
      theData.df = data.df[condt,-c(1:2)]
      theGroups  = data.df[condt, "group"]

      cytokine.splsda =
        splsda(theData.df, theGroups, scale=TRUE, ncomp=comp.num, keepX = c(var.num, var.num))

      splsda.predict = predict( cytokine.splsda, theData.df, dist="max.dist")
      Prediction1 = cbind( original = theGroups, splsda.predict$class$max.dist )
      accuracy1 = sum(Prediction1[,1] == Prediction1[,2])/length(Prediction1[,1]); accuracy1
      acc1 = 100*signif(accuracy1, digits = 2)

      # Creating a shaded predicted background using max.dist
      bg.maxdist <-  background.predict(cytokine.splsda,
                                        comp.predicted = 1,
                                        dist = 'max.dist')


      group_factors = sort(unique(theGroups))

      # Conditions to have either ellipses or prediction background or both or neither on graphs.
      if(ellipse == TRUE & bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=TRUE, title= paste (Title, "With Accuracy:",acc1,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=TRUE, title= paste (Title, "With Accuracy:",acc1,"%"))
      }
      else if(bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=FALSE, title= paste (Title, "With Accuracy:",acc1,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == FALSE & bg == FALSE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=FALSE, title= paste (Title, "With Accuracy:",acc1,"%"))
      }

      if(!is.null(style) && comp.num == 3 && (style == "3D" || style == "3d")){
        # 3D scatter plot using plot3D package
        cytokine.scores <- cytokine.splsda$variates$X
        plot3D::scatter3D(cytokine.scores[,1], cytokine.scores[,2], cytokine.scores[,3], pch=pch.values, col=colors,
                          xlab="Component 1", ylab="Component 2", zlab="Component 3", main=paste("3D Plot", ":", Title),
                          theta=20, phi=30, bty="g", colkey=FALSE)
      }
      else if(!is.null(style)){
        stop("Please enter a valid style for 3D plot:'3d' or '3D' or enter valid number of components.")
      }
      ## Cross-validation methods
      # LOOCV Different method
      if(!is.null(cv.opt)) {
        # LOOCV Different method
        if(cv.opt == "loocv"){
          set.seed(123) # For reproducibility
          loocv_results = perf(cytokine.splsda, validation = "loo")
          loocv_error_rate = loocv_results$error.rate$overall["comp1", "max.dist"]
          loocv_acc = 1 - loocv_error_rate; loocv_acc
          loocv_acc = 100*signif(loocv_acc, digits = 2)
          print(paste0(theTrts," ",unique(theGroups)[1],"vs",unique(theGroups)[2]," LOOCV Accuracy: ", loocv_acc))

          # Extracting data for plotting
          error_rates <- loocv_results$error.rate$overall[,"max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

          # Plotting with ggplot2
          a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("LOOCV Error Rate", ":", theTrt),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
        #Mfold
        else if(cv.opt == "Mfold"){
          set.seed(123) # For reproducibility
          fold_results = perf(cytokine.splsda, validation = "Mfold", folds = fold.num, nrepeat = 1000)
          fold_error_rate = fold_results$error.rate$overall["comp1", "max.dist"]
          fold_acc = 1 - fold_error_rate; fold_acc
          fold_acc = 100*signif(fold_acc, digits = 2)
          print(paste0(theTrts," ",unique(theGroups)[1],"vs",unique(theGroups)[2]," Mfold Accuracy: ", fold_acc))

          # Extracting data for plotting
          error_rates <- fold_results$error.rate$overall[,"max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

          # Plotting with ggplot2
          a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("Mfold Error Rate", ":", theTrt),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
      }

      num_rows <- nrow(cytokine.splsda$loadings$X)
      loading.vip.df <- data.frame(matrix(ncol = 0, nrow = num_rows))

      # Loop over the number of components specified by 'comp.num'
      for (comp in 1:comp.num) {
        # Plot loadings for the current component
        plotLoadings(cytokine.splsda, comp=comp, contrib='max', method='mean', size.name=1,
                     size.legend=1, legend.color=colors, title=paste("Component", comp, ":", Title), size.title=1,
                     legend.title = "Group")

        # Add columns to the data frame for the loadings and VIP for the current component
        loading.vip.df[[paste("loading.comp", comp, sep = "")]] <- cytokine.splsda$loadings$X[, comp]
        loading.vip.df[[paste("vip.comp", comp, sep = "")]] <- vip(cytokine.splsda)
      }

      # Update the column names for loadings and VIP for all components
      colnames(loading.vip.df) <- unlist(lapply(1:comp.num, function(comp) {
        c(paste("loading.comp", comp, sep = ""), paste("vip.comp", comp, sep = ""))
      }))


      # Calculate VIP scores once
      all_VIP_scores <- vip(cytokine.splsda)

      # Use all_VIP_scores for plotting
      for (comp in 1:comp.num) {
        Vscore <- as.data.frame(all_VIP_scores[, comp, drop = FALSE])
        Vscore$metabo <- rownames(Vscore)
        Vscore$comp <- Vscore[,1]  # assuming the first column contains the VIP scores
        bar <- Vscore[,c('metabo','comp')]
        rownames(bar) <- 1:nrow(bar)
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

        print(a)  # This prints the plot for each component
      }

      # Conduct further analysis based on a specific component's VIP score
      # For example, use the VIP scores from the first component
      Vscore <- as.data.frame(all_VIP_scores[, 1, drop = FALSE])
      Vscore$metabo <- rownames(Vscore)
      Vscore$comp <- Vscore[,1]
      condtVariable <- Vscore$comp > 1
      KeepX <- sum(condtVariable)
      theData.mat <- theData.df[, condtVariable, drop = FALSE]
      cytokine.splsda2 <- splsda(theData.mat, theGroups, scale=TRUE, ncomp=comp.num, keepX = c(KeepX, KeepX))


      splsda.predict = predict( cytokine.splsda2, theData.mat, dist="max.dist")
      Prediction2 = cbind( original = theGroups, splsda.predict$class$max.dist )
      accuracy2 = sum(Prediction2[,1] == Prediction2[,2])/length(Prediction2[,1]); accuracy2
      acc2 = 100*signif(accuracy2, digits = 2)

      # Creating a shaded predicted background using max.dist
      bg.maxdist <-  background.predict(cytokine.splsda2,
                                        comp.predicted = 1,
                                        dist = 'max.dist')



      # Conditions to have either ellipses or prediction background or both or neither on graphs.
      if(ellipse == TRUE & bg == TRUE){
        plotIndiv(cytokine.splsda2, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=TRUE, title= paste (Title, "(VIP>1)", "With Accuracy:",acc2,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == TRUE) {
        plotIndiv(cytokine.splsda2, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=TRUE, title=paste(Title, "(VIP>1)", "With Accuracy:",acc2,"%"))
      }
      else if(bg == TRUE){
        plotIndiv(cytokine.splsda2, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=FALSE, title=paste(Title, "(VIP>1)", "With Accuracy:",acc2,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == FALSE & bg == FALSE){
        plotIndiv(cytokine.splsda2, ind.names=NA, legend=TRUE, col=colors, pch=pch.values,
                  pch.levels = group_factors,
                  ellipse=FALSE, title= paste (Title, "(VIP>1)", "With Accuracy:",acc2,"%"))
      }
      if(!is.null(style) && comp.num == 3 && (style == "3D" || style == "3d")){
        # 3D scatter plot using plot3D package
        cytokine.scores2 <- cytokine.splsda2$variates$X
        plot3D::scatter3D(cytokine.scores2[,1], cytokine.scores2[,2], cytokine.scores2[,3], pch=pch.values, col=colors,
                          xlab="Component 1", ylab="Component 2", zlab="Component 3", main=paste("3D Plot", ":", Title, "(VIP>1)"),
                          theta=20, phi=30, bty="g", colkey=FALSE)
      }
      else if(!is.null(style)){
        stop("Please enter a valid style for 3D plot:'3d' or '3D' or enter valid number of components.")
      }
      ## Cross-validation methods
      if(!is.null(cv.opt)){
        # LOOCV Different method
        if(cv.opt == "loocv"){
          set.seed(123) # For reproducibility
          loocv_results2 = perf(cytokine.splsda2, validation = "loo")
          loocv_error_rate2 = loocv_results2$error.rate$overall["comp1", "max.dist"]
          loocv_acc2 = 1 - loocv_error_rate2; loocv_acc2
          loocv_acc2 = 100*signif(loocv_acc2, digits = 2)
          print(paste0(theTrts," ",unique(theGroups)[1],"vs",unique(theGroups)[2]," LOOCV Accuracy (VIP>1) Cytokines: ", loocv_acc2))

          # Extracting data for plotting
          error_rates <- loocv_results2$error.rate$overall[,"max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

          # Plotting with ggplot2
          a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("LOOCV Error Rate VIP > 1", ":", theTrt),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
        #Mfold
        else if(cv.opt == "Mfold"){
          set.seed(123) # For reproducibility
          fold_results2 = perf(cytokine.splsda2, validation = "Mfold", folds = fold.num, nrepeat = 1000)
          fold_error_rate2 = fold_results2$error.rate$overall["comp1", "max.dist"]
          fold_acc2 = 1 - fold_error_rate2; fold_acc2
          fold_acc2 = 100*signif(fold_acc2, digits = 2)
          print(paste0(theTrts," ",unique(theGroups)[1],"vs",unique(theGroups)[2]," Mfold Accuracy (VIP>1) Cytokines: ", fold_acc2))

          # Extracting data for plotting
          error_rates <- fold_results2$error.rate$overall[,"max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(error_df, id.vars = "Component", variable.name = "Distance", value.name = "ErrorRate")

          # Plotting with ggplot2
          a <- ggplot(error_df, aes(x = Component, y = ErrorRate, color = Distance, group = 1)) +
            geom_line() +
            geom_point(size = 3) +
            labs(title = paste("LOOCV Error Rate VIP > 1", ":", theTrt),
                 x = "Number of Components",
                 y = "Error Rate") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_color_manual(values = "red", labels = "max.dist")
          print(a)
        }
      }

      num_rows2 <- nrow(cytokine.splsda2$loadings$X)
      loading.vip.df2 <- data.frame(matrix(ncol = 0, nrow = num_rows2))

      # Loop over the number of components specified by 'comp.num'
      for (comp in 1:comp.num) {
        # Plot loadings for the current component
        plotLoadings(cytokine.splsda2, comp=comp, contrib='max', method='mean', size.name=1,
                     size.legend=1, legend.color=colors, title=paste("Component", comp, ":", Title), size.title=1,
                     legend.title = "Group")

        # Add columns to the data frame for the loadings and VIP for the current component
        loading.vip.df2[[paste("loading.comp", comp, sep = "")]] <- cytokine.splsda2$loadings$X[, comp]
        loading.vip.df2[[paste("vip.comp", comp, sep = "")]] <- vip(cytokine.splsda2)
      }

      # Update the column names for loadings and VIP for all components
      colnames(loading.vip.df2) <- unlist(lapply(1:comp.num, function(comp) {
        c(paste("loading.comp", comp, sep = ""), paste("vip.comp", comp, sep = ""))
      }))
      # Prints confusion matrix
      if(conf.mat == TRUE){
        print(paste0(theTrts," ",unique(theGroups)[1]," vs ",unique(theGroups)[2]," Confusion Matrix for PLS-DA Comparison"))
        print(get.confusion_matrix(truth = Prediction1[,1], predicted = Prediction1[,2])) # Confusion matrix for all variables in model
        print(paste0(theTrts," ",unique(theGroups)[1]," vs ",unique(theGroups)[2]," Confusion Matrix for PLS-DA Comparison with VIP Score > 1"))
        print(get.confusion_matrix(truth = Prediction2[,1], predicted = Prediction2[,2])) # Confusion matrix for all variables with VIP Score > 1
      }
    }
    dev.off()
  }
}

