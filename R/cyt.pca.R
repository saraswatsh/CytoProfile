#' Analyze data with Principal Component Analysis (PCA) for cytokines.
#'
#' @param x.df A data frame with cytokine data.
#' @param colors Vector of colors to be set, list of colors to be set equal to the number of groups or treatments. Default set to NULL and will result in random colors.
#' @param title Title of pdf file to be saved with figures.
#' @param ellipse 95% confidence region. Default set to FALSE.
#' @param comp.num Number of components to be used in PCA. Default set to 2.
#' @param scale Scale the data. Default set to NULL. Uses log2 transformation.
#' @param pch.values Vector of plotting characters to be used. Default set to NULL.
#' @param style Use "3d" or "3D" to generate a 3D plot. Default set to NULL.
#'
#' @return A pdf file with PCA plots (2D and 3D), loadings plots, biplots, and correlation circle plots.
#' @export
#'
#' @examples
#' data <- cytodata[,-c(1,4)]
#' data.df <- filter(data, Group != "ND" & Treatment != "Unstimulated")
#' data.df <- data.df[,-22]
#' cyt.pca(data.df, title = "PCA_Example_Analysis.pdf" ,colors = c("black", "red2"), scale = "log2", comp.num = 3, pch.values = c(16,4), style = "3D")
#' cyt.pca(data.df, title = "PCA_Example_Analysis2.pdf" ,colors = c("black", "red2"), scale = "log2", comp.num = 2, ellipse = TRUE, pch.values = c(16,4))

cyt.pca <- function(data.df, colors = NULL, title, ellipse = FALSE, comp.num = 2, scale = NULL, pch.values = NULL, style = NULL) {
  if(!is.null(scale) && scale == "log2"){
    # Log2 transforming cytokines
    data.df <- data.frame(data.df[,c(1:2)], log2(data.df[, -c(1:2)]))
    print("Results based on log2 transformation:")
  } else {
    print("Results based on no transformation:")
  }

  # Making the first two columns to be lowercase
  names(data.df)[1:2] <- tolower(names(data.df)[1:2])
  # Creating a table to have two separate vectors for group and stimulation
  a <- table(data.df[, c(1, 2)])

  Group.vec <- dimnames(a)$group
  if("treatment" %in% names(data.df)[1:2]){
    Treatment.vec <- dimnames(a)$treatment
  } else {
    Stimulation.vec <- dimnames(a)$stimulation
  }
  # Generate a color palette based on the number of groups
  if (is.null(colors)) {
    num_groups <- length(unique(data.df$group))
    colors <- rainbow(num_groups)
  }

  if("treatment" %in% names(data.df)[1:2]){
    pdf(file=title)
    for(i in 1:length(Treatment.vec)) {
      theTrt <- Treatment.vec[i]
      condt  <- data.df[, "treatment"] == theTrt
      theData.df <- data.df[condt, -c(1:2)]
      theGroups  <- data.df[condt, "group"]

      cytokine.pca <- pca(theData.df, ncomp = comp.num, center = TRUE, scale = TRUE)

      group_factors <- sort(unique(theGroups))

      # Plotting the PCA using plotIndiv
      plotIndiv(cytokine.pca, group = theGroups, ind.names = FALSE, legend = TRUE,
                col = colors, title = paste("PCA of", theTrt), ellipse = ellipse,
                pch = pch.values, pch.levels = group_factors)

      if(!is.null(style) && comp.num == 3 && (style == "3D" || style == "3d")){
        # 3D scatter plot using plot3D package
        cytokine.scores <- cytokine.pca$variates$X
        plot3D::scatter3D(cytokine.scores[,1], cytokine.scores[,2], cytokine.scores[,3], pch=pch.values, col=colors,
                          xlab="PC1", ylab="PC2", zlab="PC3", main=paste("3D Plot", ":", theTrt),
                          theta=20, phi=30, bty="g", colkey=FALSE)
      } else if(!is.null(style)){
        stop("Please enter a valid style for 3D plot:'3d' or '3D' or enter valid number of components.")
      }

      # Scree Plot
      variances <- cytokine.pca$prop_expl_var$X
      cumulative_variances <- cytokine.pca$cum.var

      scree_data <- data.frame(
        Component = 1:comp.num,
        Variance = variances,
        Cumulative = cumulative_variances
      )

      scree_plot <- ggplot(scree_data, aes(x = Component)) +
        geom_line(aes(y = Variance, color = "Individual"), size = 1) +
        geom_point(aes(y = Variance, color = "Individual"), size = 2) +
        geom_line(aes(y = Cumulative, color = "Cumulative"), size = 1, linetype = "dashed") +
        geom_point(aes(y = Cumulative, color = "Cumulative"), size = 2) +
        scale_color_manual(values = c("Individual" = "blue", "Cumulative" = "green")) +
        labs(title = paste("Scree Plot", ":", theTrt), x = "Principal Components", y = "Explained Variance", color = "Variance Type") +
        theme_minimal() +
        scale_x_continuous(breaks = 1:comp.num) +
        geom_text(aes(y = Variance, label = paste0(round(Variance * 100, 1), "%")), vjust = -1.5, hjust = 0.5, size = 4) +
        geom_text(aes(y = Cumulative, label = paste0(round(Cumulative * 100, 1), "%")), vjust = 1.5, hjust = 0.5, size = 4)

      print(scree_plot)

      # Plotting the loadings
      for(comp in 1:comp.num) {
        plotLoadings(cytokine.pca, comp = 1, size.names = 1, size.legend = 1, col = "grey", legend.color = colors,
                     title = paste("Component", comp, ":", theTrt), legend = TRUE)
      }

      # Creating biplot using mixOmics biplot function
      b <- biplot(cytokine.pca, comp = c(1, 2), group = theGroups, col.per.group = colors, var.arrow.col = "blue",
                  var.arrow.size = 0.5, var.arrow.length = 0.2, var.names = TRUE, var.names.col = "blue",
                  var.names.size = 3, ind.names = FALSE, legend = TRUE, legend.title = "Group")
      print(b)
      # Variable plot: correlation circle plot
      c <- plotVar(cytokine.pca, comp = c(1,2), var.names = TRUE, cex = 4, col = "black",
                   overlap = T, title = paste("Correlation Circle Plot", ":", theTrt), style = "ggplot2")
    }
    dev.off()
  } else {
    pdf(file=title)
    for(i in 1:length(Stimulation.vec)) {
      theTrts <- Stimulation.vec[i]
      condt  <- data.df[, "stimulation"] == theTrt
      theData.df <- data.df[condt, -c(1:2)]
      theGroups  <- data.df[condt, "group"]

      cytokine.pca <- pca(theData.df, ncomp = comp.num, center = TRUE, scale = TRUE)

      group_factors <- sort(unique(theGroups))

      # Plotting the PCA using plotIndiv
      plotIndiv(cytokine.pca, group = theGroups, ind.names = FALSE, legend = TRUE,
                col = colors, title = paste("PCA of", theTrt), ellipse = ellipse,
                pch = pch.values, pch.levels = group_factors)

      if(!is.null(style) && comp.num == 3 && (style == "3D" || style == "3d")){
        # 3D scatter plot using plot3D package
        cytokine.scores <- cytokine.pca$variates$X
        plot3D::scatter3D(cytokine.scores[,1], cytokine.scores[,2], cytokine.scores[,3], pch=pch.values, col=colors,
                          xlab="PC1", ylab="PC2", zlab="PC3", main=paste("3D Plot", ":", theTrt),
                          theta=20, phi=30, bty="g", colkey=FALSE)
      } else if(!is.null(style)){
        stop("Please enter a valid style for 3D plot:'3d' or '3D' or enter valid number of components.")
      }

      # Scree Plot
      variances <- cytokine.pca$prop_expl_var$X
      cumulative_variances <- cytokine.pca$cum.var

      scree_data <- data.frame(
        Component = 1:comp.num,
        Variance = variances,
        Cumulative = cumulative_variances
      )

      scree_plot <- ggplot(scree_data, aes(x = Component)) +
        geom_line(aes(y = Variance, color = "Individual"), size = 1) +
        geom_point(aes(y = Variance, color = "Individual"), size = 2) +
        geom_line(aes(y = Cumulative, color = "Cumulative"), size = 1, linetype = "dashed") +
        geom_point(aes(y = Cumulative, color = "Cumulative"), size = 2) +
        scale_color_manual(values = c("Individual" = "blue", "Cumulative" = "green")) +
        labs(title = paste("Scree Plot", ":", theTrt), x = "Principal Components", y = "Explained Variance", color = "Variance Type") +
        theme_minimal() +
        scale_x_continuous(breaks = 1:comp.num) +
        geom_text(aes(y = Variance, label = paste0(round(Variance * 100, 1), "%")), vjust = -1.5, hjust = 0.5, size = 4) +
        geom_text(aes(y = Cumulative, label = paste0(round(Cumulative * 100, 1), "%")), vjust = 1.5, hjust = 0.5, size = 4)

      print(scree_plot)

      # Plotting the loadings
      for(comp in 1:comp.num) {
        plotLoadings(cytokine.pca, comp = 1, size.names = 1, size.legend = 1, col = "grey", legend.color = colors,
                     title = paste("Component", comp, ":", theTrt), legend = TRUE)
      }

      # Creating biplot using mixOmics biplot function
      b <- biplot(cytokine.pca, comp = c(1, 2), group = theGroups, col.per.group = colors, var.arrow.col = "blue",
                  var.arrow.size = 0.5, var.arrow.length = 0.2, var.names = TRUE, var.names.col = "blue",
                  var.names.size = 3, ind.names = FALSE, legend = TRUE, legend.title = "Group")
      print(b)
      # Variable plot: correlation circle plot
      c <- plotVar(cytokine.pca, comp = c(1,2), var.names = TRUE, cex = 4, col = "black",
                   overlap = T, title = paste("Correlation Circle Plot", ":", theTrt), style = "ggplot2")
    }
    dev.off()
  }
}
