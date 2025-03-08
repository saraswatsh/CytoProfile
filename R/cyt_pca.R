#' Analyze Data with Principal Component Analysis (PCA) for Cytokines.
#'
#' @param data A data frame containing cytokine data. It should include at
#' least one column representing grouping information and optionally a second
#' column representing treatment or stimulation.
#' @param group_col Character. The name of the column containing the
#'  grouping information. If not specified and \code{trt_col} is provided, the
#'   treatment column will be used as the grouping variable.
#' @param trt_col Character. The name of the column containing the
#' treatment information. If not specified and \code{group_col} is provided,
#'  the grouping column will be used as the treatment variable.
#' @param colors A vector of colors corresponding to the groups.
#'  If set to NULL, a palette is generated using \code{rainbow()} based on the
#'   number of unique groups.
#' @param pdf_title A string specifying the file name of the PDF where the
#'  PCA plots will be saved.
#' @param ellipse Logical. If TRUE, a 95% confidence ellipse is drawn on the
#'  PCA individuals plot. Default is FALSE.
#' @param comp_num Numeric. The number of principal components to compute and
#'  display. Default is 2.
#' @param scale Character. If set to "log2", a log2 transformation is applied
#'  to the numeric cytokine measurements (excluding the grouping columns).
#'  Default is NULL.
#' @param pch_values A vector of plotting symbols (pch values) to be used in
#' the PCA plots. Default is NULL.
#' @param style Character. If set to "3d" or "3D" and \code{comp_num} equals 3,
#' a 3D scatter plot is generated using the plot3D package. Default is NULL.
#'
#' @description
#' This function performs Principal Component Analysis (PCA) on cytokine data
#' and generates several types of plots,
#' including:
#' \itemize{
#'   \item 2D PCA plots using mixOmics's \code{plotIndiv} function,
#'   \item 3D scatter plots (if \code{style} is "3d" or "3D" and
#'   \code{comp_num} is 3) via the plot3D package,
#'   \item Scree plots showing both individual and cumulative
#'   explained variance,
#'   \item Loadings plots, and
#'   \item Biplots and correlation circle plots.
#' }
#' The function optionally applies a log2 transformation to the numeric
#' data and handles analyses based treatment groups.
#'
#' @return A PDF file containing the PCA plots is generated and saved.
#'
#' @export
#' @import mixOmics
#' @import ggplot2
#' @import plot3D
#'
#' @examples
#' # Load sample data
#' data("cytodata")
#' # Subset data to exclude columns 1, 4, and 24, then filter out rows
#' # where Group is "ND" and Treatment is "Unstimulated"
#' data_subset <- cytodata[, -c(1, 4, 24)]
#' data_df <- dplyr::filter(data_subset,
#' Group != "ND" & Treatment != "Unstimulated")
#' # Run PCA analysis and save plots to a PDF file
#' cyt_pca(
#'   data = data_df,
#'   pdf_title = "Example_PCA_Analysis.pdf",
#'   colors = c("black", "red2"),
#'   scale = "log2",
#'   comp_num = 3,
#'   pch_values = c(16, 4),
#'   style = "3D",
#'   group_col = "Group",
#'   trt_col = "Treatment",
#'   ellipse = FALSE
#' )
#'
cyt_pca <- function(data, group_col = NULL, trt_col = NULL,
                    colors = NULL, pdf_title, ellipse = FALSE,
                    comp_num = 2, scale = NULL, pch_values = NULL,
                    style = NULL) {
  # If one factor is missing, use the provided column for both grouping
  # and treatment.
  if (is.null(group_col) && !is.null(trt_col)) {
    message("No group column provided; using the treatment column as
            the grouping variable.")
    group_col <- trt_col
  }
  if (is.null(trt_col) && !is.null(group_col)) {
    message("No treatment column provided; using the group column as
            the treatment variable.")
    trt_col <- group_col
  }
  if (is.null(group_col) && is.null(trt_col)) {
    stop("At least one factor column must be provided.")
  }

  # Optionally apply log2 transformation only to numeric columns
  if (!is.null(scale) && scale == "log2") {
    # Identify numeric columns not corresponding to the factor columns
    numeric_idx <- sapply(data, is.numeric)
    # Exclude the group and treatment columns (if present, even if numeric)
    numeric_idx[names(data) %in% unique(c(group_col, trt_col))] <- FALSE
    if (sum(numeric_idx) == 0) {
      warning("No numeric columns available for log2 transformation.")
    }
    data <- data.frame(
      data[, unique(c(group_col, trt_col)), drop = FALSE],
      log2(data[, numeric_idx, drop = FALSE])
    )
    print("Results based on log2 transformation:")
  } else {
    print("Results based on no transformation:")
  }

  # Convert factor column names to lowercase for consistency
  names(data)[names(data) %in% unique(c(group_col, trt_col))] <-
    tolower(names(data)[names(data) %in% unique(c(group_col, trt_col))])
  group_col <- tolower(group_col)
  trt_col <- tolower(trt_col)

  # Generate a color palette if not provided (based on the
  # grouping variable levels)
  if (is.null(colors)) {
    num_groups <- length(unique(data[[group_col]]))
    colors <- rainbow(num_groups)
  }

  pdf(file = pdf_title, width = 8.5, height = 8)

  # Case 1: Overall PCA when both factors are the same.
  if (group_col == trt_col) {
    overall_analysis <- "Overall Analysis"

    # Remove the factor column(s) from predictors and keep only numeric columns
    the_data_df <- data[, !(names(data) %in% unique(c(group_col, trt_col)))]
    the_data_df <- the_data_df[, sapply(the_data_df, is.numeric)]

    the_groups <- as.vector(data[[group_col]])
    if (length(unique(the_groups)) < 2) {
      stop("The grouping variable must have at least two levels for PCA.
           Please provide an appropriate grouping column.")
    }

    cytokine_pca <- mixOmics::pca(the_data_df,
      ncomp = comp_num,
      center = TRUE, scale = TRUE
    )

    group_factors <- sort(unique(the_groups))

    # Plot PCA individuals plot
    plotIndiv(cytokine_pca,
      group = the_groups, ind.names = FALSE, legend = TRUE,
      col = colors, title = paste("PCA:", overall_analysis), ellipse = ellipse,
      pch = pch_values, pch.levels = group_factors
    )

    # 3D Plot if requested and exactly three components are used
    if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
      cytokine_scores <- cytokine_pca$variates$X
      plot3D::scatter3D(cytokine_scores[, 1], cytokine_scores[, 2],
        cytokine_scores[, 3],
        pch = pch_values, col = colors,
        xlab = "PC1", ylab = "PC2", zlab = "PC3",
        main = paste("3D Plot:", overall_analysis),
        theta = 20, phi = 30, bty = "g", colkey = FALSE
      )
    } else if (!is.null(style)) {
      stop("Please enter a valid style for 3D plot: '3d' or '3D' or
           enter a valid number of components.")
    }

    # Scree Plot
    variances <- cytokine_pca$prop_expl_var$X
    cumulative_variances <- cytokine_pca$cum.var
    scree_data <- data.frame(
      Component = 1:comp_num,
      Variance = variances,
      Cumulative = cumulative_variances
    )

    scree_plot <- ggplot(scree_data, aes(x = Component)) +
      geom_line(aes(y = Variance, color = "Individual"), size = 1) +
      geom_point(aes(y = Variance, color = "Individual"), size = 2) +
      geom_line(aes(y = Cumulative, color = "Cumulative"),
        size = 1, linetype = "dashed"
      ) +
      geom_point(aes(y = Cumulative, color = "Cumulative"), size = 2) +
      scale_color_manual(values = c(
        "Individual" = "blue",
        "Cumulative" = "green"
      )) +
      labs(
        title = paste("Scree Plot:", overall_analysis),
        x = "Principal Components",
        y = "Explained Variance", color = "Variance Type"
      ) +
      theme_minimal() +
      scale_x_continuous(breaks = 1:comp_num) +
      geom_text(
        aes(y = Variance, label = paste0(
          round(Variance * 100, 1),
          "%"
        )),
        vjust = -1.5, hjust = 0.5, size = 4
      ) +
      geom_text(
        aes(y = Cumulative, label = paste0(
          round(Cumulative * 100, 1),
          "%"
        )),
        vjust = 1.5, hjust = 0.5, size = 4
      )

    print(scree_plot)

    # Plot loadings for each component
    for (comp in 1:comp_num) {
      plotLoadings(cytokine_pca,
        comp = comp, size.names = 1, size.legend = 1, col = "grey",
        legend.color = colors, title = paste(
          "Loadings Component", comp, ":",
          overall_analysis
        ),
        legend = TRUE
      )
    }

    # Biplot for components 1 and 2
    biplot_obj <- biplot(cytokine_pca,
      comp = c(1, 2), group = the_groups, col.per.group = colors,
      var.arrow.col = "blue", var.arrow.size = 0.5, var.arrow.length = 0.2,
      var.names = TRUE, var.names.col = "blue", var.names.size = 3,
      ind.names = FALSE, legend = TRUE, legend.title = "Group"
    )
    print(biplot_obj)

    # Correlation circle plot
    mixOmics::plotVar(cytokine_pca,
      comp = c(1, 2), var.names = TRUE, cex = 4, col = "black",
      overlap = TRUE, title = paste(
        "Correlation Circle Plot:",
        overall_analysis
      ),
      style = "ggplot2"
    )
  } else {
    # Case 2: When grouping and treatment columns differ
    levels_vec <- unique(data[[trt_col]])
    for (i in seq_along(levels_vec)) {
      current_level <- levels_vec[i]
      title_sub <- current_level
      condt <- data[[trt_col]] == current_level

      the_data_df <- data[condt, !(names(data) %in% unique(c(
        group_col,
        trt_col
      )))]
      the_data_df <- the_data_df[, sapply(the_data_df, is.numeric)]
      the_groups <- as.vector(data[condt, group_col])

      if (length(unique(the_groups)) < 2) {
        stop("The grouping variable must have at least two levels for PCA.
             Please provide an appropriate grouping column.")
      }

      cytokine_pca <- mixOmics::pca(the_data_df,
        ncomp = comp_num,
        center = TRUE, scale = TRUE
      )

      group_factors <- sort(unique(the_groups))

      # PCA individuals plot for current treatment level
      plotIndiv(cytokine_pca,
        group = the_groups, ind.names = FALSE, legend = TRUE,
        col = colors, title = paste("PCA:", title_sub),
        ellipse = ellipse, pch = pch_values, pch.levels = group_factors
      )

      # 3D Plot if applicable
      if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
        cytokine_scores <- cytokine_pca$variates$X
        plot3D::scatter3D(cytokine_scores[, 1], cytokine_scores[, 2],
          cytokine_scores[, 3],
          pch = pch_values, col = colors,
          xlab = "PC1", ylab = "PC2", zlab = "PC3",
          main = paste("3D Plot:", current_level),
          theta = 20, phi = 30, bty = "g", colkey = FALSE
        )
      } else if (!is.null(style)) {
        stop("Please enter a valid style for 3D plot: '3d' or '3D' or
             enter a valid number of components.")
      }

      # Scree Plot for the current treatment level
      variances <- cytokine_pca$prop_expl_var$X
      cumulative_variances <- cytokine_pca$cum.var
      scree_data <- data.frame(
        Component = 1:comp_num,
        Variance = variances,
        Cumulative = cumulative_variances
      )

      scree_plot <- ggplot(scree_data, aes(x = Component)) +
        geom_line(aes(y = Variance, color = "Individual"), size = 1) +
        geom_point(aes(y = Variance, color = "Individual"), size = 2) +
        geom_line(aes(y = Cumulative, color = "Cumulative"),
          size = 1, linetype = "dashed"
        ) +
        geom_point(aes(y = Cumulative, color = "Cumulative"), size = 2) +
        scale_color_manual(values = c(
          "Individual" = "blue",
          "Cumulative" = "green"
        )) +
        labs(
          title = paste("Scree Plot:", current_level),
          x = "Principal Components", y = "Explained Variance",
          color = "Variance Type"
        ) +
        theme_minimal() +
        scale_x_continuous(breaks = 1:comp_num) +
        geom_text(
          aes(y = Variance, label = paste0(
            round(Variance * 100, 1),
            "%"
          )),
          vjust = -1.5, hjust = 0.5, size = 4
        ) +
        geom_text(
          aes(y = Cumulative, label = paste0(
            round(Cumulative * 100, 1),
            "%"
          )),
          vjust = 1.5, hjust = 0.5, size = 4
        )

      print(scree_plot)

      # Plot loadings for each component
      for (comp in 1:comp_num) {
        plotLoadings(cytokine_pca,
          comp = comp, size.names = 1, size.legend = 1, col = "grey",
          legend.color = colors, title = paste(
            "Loadings Component", comp, ":",
            current_level
          ),
          legend = TRUE
        )
      }

      # Biplot for components 1 and 2
      biplot_obj <- biplot(cytokine_pca,
        comp = c(1, 2), group = the_groups, col.per.group = colors,
        var.arrow.col = "blue", var.arrow.size = 0.5, var.arrow.length = 0.2,
        var.names = TRUE, var.names.col = "blue", var.names.size = 3,
        ind.names = FALSE, legend = TRUE, legend.title = "Group"
      )
      print(biplot_obj)

      # Correlation circle plot
      mixOmics::plotVar(cytokine_pca,
        comp = c(1, 2), var.names = TRUE, cex = 4, col = "black",
        overlap = TRUE, title = paste(
          "Correlation Circle Plot:",
          current_level
        ),
        style = "ggplot2"
      )
    }
  }

  if (dev.cur() > 1) dev.off()
}
