#' Analyze data with Principal Component Analysis (PCA) for cytokines.
#'
#' @param data A data frame containing cytokine data. It should include at least one column representing grouping information and optionally a second column representing treatment or stimulation.
#' @param group.col Character. The name of the column containing the grouping information. If not specified and \code{trt.col} is provided, the treatment column will be used as the grouping variable.
#' @param trt.col Character. The name of the column containing the treatment (or stimulation) information. If not specified and \code{group.col} is provided, the grouping column will be used as the treatment variable.
#' @param colors A vector of colors corresponding to the groups. If set to NULL, a palette is generated using \code{rainbow()} based on the number of unique groups.
#' @param title A string specifying the file name of the PDF where the PCA plots will be saved.
#' @param ellipse Logical. If TRUE, a 95% confidence ellipse is drawn on the PCA individuals plot. Default is FALSE.
#' @param comp.num Numeric. The number of principal components to compute and display. Default is 2.
#' @param scale Character. If set to "log2", a log2 transformation is applied to the numeric cytokine measurements (excluding the grouping columns). Default is NULL.
#' @param pch.values A vector of plotting symbols (pch values) to be used in the PCA plots. Default is NULL.
#' @param style Character. If set to "3d" or "3D" and \code{comp.num} equals 3, a 3D scatter plot is generated using the plot3D package. Default is NULL.
#'
#' @description
#' This function performs Principal Component Analysis (PCA) on cytokine data and generates several types of plots,
#' including:
#' \itemize{
#'   \item 2D PCA plots using mixOmics's \code{plotIndiv} function,
#'   \item 3D scatter plots (if \code{style} is "3d" or "3D" and \code{comp.num} is 3) via the plot3D package,
#'   \item Scree plots showing both individual and cumulative explained variance,
#'   \item Loadings plots, and
#'   \item Biplots and correlation circle plots.
#' }
#' The function optionally applies a log2 transformation to the numeric data and handles analyses based on either treatment
#' or stimulation groups.
#'
#' @return A PDF file containing the PCA plots is generated and saved.
#'
#' @export
#' @import mixOmics
#' @import ggplot2
#' @import plot3D
#'
#' @examples
#' \dontrun{
#' data <- cytodata[, -c(1, 4, 24)]
#' data.df <- filter(data, Group != "ND" & Treatment != "Unstimulated")
#' cyt.pca(data.df,
#'   title = "Example PCA Analysis.pdf", colors = c("black", "red2"),
#'   scale = "log2", comp.num = 3, pch.values = c(16, 4), style = "3D",
#'   group.col = "Group", trt.col = "Treatment"
#' )
#' cyt.pca(data.df,
#'   title = "Example PCA Analysis 2.pdf", colors = c("black", "red2"),
#'   scale = "log2", comp.num = 2, pch.values = c(16, 4), group.col = "Group"
#' )
#' }
cyt.pca <- function(data, group.col = NULL, trt.col = NULL, colors = NULL, title,
                    ellipse = FALSE, comp.num = 2, scale = NULL, pch.values = NULL, style = NULL) {
  # If one factor is missing, use the provided column for both grouping and treatment.
  if (is.null(group.col) && !is.null(trt.col)) {
    message("No group column provided; using the treatment column as the grouping variable.")
    group.col <- trt.col
  }
  if (is.null(trt.col) && !is.null(group.col)) {
    message("No treatment column provided; using the group column as the treatment variable.")
    trt.col <- group.col
  }
  if (is.null(group.col) && is.null(trt.col)) {
    stop("At least one factor column must be provided.")
  }

  # Optionally apply log2 transformation only to numeric columns
  if (!is.null(scale) && scale == "log2") {
    # Identify numeric columns not corresponding to the factor columns
    numeric_idx <- sapply(data, is.numeric)
    # Exclude the group and treatment columns (if present, even if numeric)
    numeric_idx[names(data) %in% unique(c(group.col, trt.col))] <- FALSE
    if (sum(numeric_idx) == 0) {
      warning("No numeric columns available for log2 transformation.")
    }
    data <- data.frame(
      data[, unique(c(group.col, trt.col)), drop = FALSE],
      log2(data[, numeric_idx, drop = FALSE])
    )
    print("Results based on log2 transformation:")
  } else {
    print("Results based on no transformation:")
  }

  # Convert factor column names to lowercase for consistency
  names(data)[names(data) %in% unique(c(group.col, trt.col))] <-
    tolower(names(data)[names(data) %in% unique(c(group.col, trt.col))])
  group.col <- tolower(group.col)
  trt.col <- tolower(trt.col)

  # Generate a color palette if not provided (based on the grouping variable levels)
  if (is.null(colors)) {
    num_groups <- length(unique(data[[group.col]]))
    colors <- rainbow(num_groups)
  }

  pdf(file = title, width = 8.5, height = 8)

  # Case 1: Overall PCA when both factors are the same.
  if (group.col == trt.col) {
    Title <- "Overall Analysis"

    # Remove the factor column(s) from predictors and keep only numeric columns
    theData.df <- data[, !(names(data) %in% unique(c(group.col, trt.col)))]
    theData.df <- theData.df[, sapply(theData.df, is.numeric)]

    theGroups <- as.vector(data[[group.col]])
    if (length(unique(theGroups)) < 2) {
      stop("The grouping variable must have at least two levels for PCA. Please provide an appropriate grouping column.")
    }

    cytokine.pca <- mixOmics::pca(theData.df, ncomp = comp.num, center = TRUE, scale = TRUE)

    group_factors <- sort(unique(theGroups))

    # Plot PCA individuals plot
    plotIndiv(cytokine.pca,
      group = theGroups, ind.names = FALSE, legend = TRUE,
      col = colors, title = paste("PCA:", Title), ellipse = ellipse,
      pch = pch.values, pch.levels = group_factors
    )

    # 3D Plot if requested and exactly three components are used
    if (!is.null(style) && comp.num == 3 && (tolower(style) == "3d")) {
      cytokine.scores <- cytokine.pca$variates$X
      plot3D::scatter3D(cytokine.scores[, 1], cytokine.scores[, 2], cytokine.scores[, 3],
        pch = pch.values, col = colors,
        xlab = "PC1", ylab = "PC2", zlab = "PC3",
        main = paste("3D Plot:", Title),
        theta = 20, phi = 30, bty = "g", colkey = FALSE
      )
    } else if (!is.null(style)) {
      stop("Please enter a valid style for 3D plot: '3d' or '3D' or enter a valid number of components.")
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
      labs(title = paste("Scree Plot:", Title), x = "Principal Components", y = "Explained Variance", color = "Variance Type") +
      theme_minimal() +
      scale_x_continuous(breaks = 1:comp.num) +
      geom_text(aes(y = Variance, label = paste0(round(Variance * 100, 1), "%")),
        vjust = -1.5, hjust = 0.5, size = 4
      ) +
      geom_text(aes(y = Cumulative, label = paste0(round(Cumulative * 100, 1), "%")),
        vjust = 1.5, hjust = 0.5, size = 4
      )

    print(scree_plot)

    # Plot loadings for each component
    for (comp in 1:comp.num) {
      plotLoadings(cytokine.pca,
        comp = comp, size.names = 1, size.legend = 1, col = "grey",
        legend.color = colors, title = paste("Loadings Component", comp, ":", Title), legend = TRUE
      )
    }

    # Biplot for components 1 and 2
    biplot_obj <- biplot(cytokine.pca,
      comp = c(1, 2), group = theGroups, col.per.group = colors,
      var.arrow.col = "blue", var.arrow.size = 0.5, var.arrow.length = 0.2,
      var.names = TRUE, var.names.col = "blue", var.names.size = 3,
      ind.names = FALSE, legend = TRUE, legend.title = "Group"
    )
    print(biplot_obj)

    # Correlation circle plot
    plotVar(cytokine.pca,
      comp = c(1, 2), var.names = TRUE, cex = 4, col = "black",
      overlap = TRUE, title = paste("Correlation Circle Plot:", Title), style = "ggplot2"
    )
  } else {
    # Case 2: When grouping and treatment columns differ
    levels.vec <- unique(data[[trt.col]])
    for (i in seq_along(levels.vec)) {
      current.level <- levels.vec[i]
      Title_sub <- current.level
      condt <- data[[trt.col]] == current.level

      theData.df <- data[condt, !(names(data) %in% unique(c(group.col, trt.col)))]
      theData.df <- theData.df[, sapply(theData.df, is.numeric)]
      theGroups <- as.vector(data[condt, group.col])

      if (length(unique(theGroups)) < 2) {
        stop("The grouping variable must have at least two levels for PCA. Please provide an appropriate grouping column.")
      }

      cytokine.pca <- mixOmics::pca(theData.df, ncomp = comp.num, center = TRUE, scale = TRUE)

      group_factors <- sort(unique(theGroups))

      # PCA individuals plot for current treatment level
      plotIndiv(cytokine.pca,
        group = theGroups, ind.names = FALSE, legend = TRUE,
        col = colors, title = paste("PCA of", current.level), ellipse = ellipse,
        pch = pch.values, pch.levels = group_factors
      )

      # 3D Plot if applicable
      if (!is.null(style) && comp.num == 3 && (tolower(style) == "3d")) {
        cytokine.scores <- cytokine.pca$variates$X
        plot3D::scatter3D(cytokine.scores[, 1], cytokine.scores[, 2], cytokine.scores[, 3],
          pch = pch.values, col = colors, xlab = "PC1", ylab = "PC2", zlab = "PC3",
          main = paste("3D Plot:", current.level),
          theta = 20, phi = 30, bty = "g", colkey = FALSE
        )
      } else if (!is.null(style)) {
        stop("Please enter a valid style for 3D plot: '3d' or '3D' or enter a valid number of components.")
      }

      # Scree Plot for the current treatment level
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
        labs(title = paste("Scree Plot:", current.level), x = "Principal Components", y = "Explained Variance", color = "Variance Type") +
        theme_minimal() +
        scale_x_continuous(breaks = 1:comp.num) +
        geom_text(aes(y = Variance, label = paste0(round(Variance * 100, 1), "%")),
          vjust = -1.5, hjust = 0.5, size = 4
        ) +
        geom_text(aes(y = Cumulative, label = paste0(round(Cumulative * 100, 1), "%")),
          vjust = 1.5, hjust = 0.5, size = 4
        )

      print(scree_plot)

      # Plot loadings for each component
      for (comp in 1:comp.num) {
        plotLoadings(cytokine.pca,
          comp = comp, size.names = 1, size.legend = 1, col = "grey",
          legend.color = colors, title = paste("Loadings Component", comp, ":", current.level), legend = TRUE
        )
      }

      # Biplot for components 1 and 2
      biplot_obj <- biplot(cytokine.pca,
        comp = c(1, 2), group = theGroups, col.per.group = colors,
        var.arrow.col = "blue", var.arrow.size = 0.5, var.arrow.length = 0.2,
        var.names = TRUE, var.names.col = "blue", var.names.size = 3,
        ind.names = FALSE, legend = TRUE, legend.title = "Group"
      )
      print(biplot_obj)

      # Correlation circle plot
      plotVar(cytokine.pca,
        comp = c(1, 2), var.names = TRUE, cex = 4, col = "black",
        overlap = TRUE, title = paste("Correlation Circle Plot:", current.level), style = "ggplot2"
      )
    }
  }

  dev.off()
}
