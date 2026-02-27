#' Analyze Data with Principal Component Analysis (PCA) for Cytokines.
#'
#' @param data A data frame containing cytokine data. It should include at
#' least one column representing grouping information and optionally a second
#' column representing treatment or stimulation.
#' @param group_col A string specifying the column name that contains the first group
#'   information. If \code{group_col2} is not provided, an overall analysis will
#'   be performed.
#' @param group_col2 A string specifying the second grouping column. Default is
#'   \code{NULL}.
#' @param colors A vector of colors corresponding to the groups.
#'  If set to NULL, a palette is generated using \code{rainbow()} based on the
#'   number of unique groups.
#' @param output_file Optional string specifying the name of the file
#'   to be created.  When `NULL` (default), plots are drawn on
#'   the current graphics device. Ensure that the file
#'   extension matches the desired format (e.g., ".pdf" for PDF output
#'   or ".png" for PNG output or .tiff for TIFF output).
#' @param ellipse Logical. If TRUE, a 95% confidence ellipse is drawn on the
#'  PCA individuals plot. Default is FALSE.
#' @param comp_num Numeric. The number of principal components to compute and
#'  display. Default is 2.
#' @param scale Character string specifying a transformation to apply to
#'   numeric variables before PCA.  Options are "none" (no
#'   transformation), "log2", "log10", "zscore", or "custom".  When
#'   "custom" is selected, a user supplied function must be given via
#'   \code{custom_fn}.  Defaults to "none".
#' @param custom_fn A custom function used when \code{scale = "custom"}.
#'   Should take a numeric vector and return a numeric vector.  Ignored
#'   otherwise.
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
#'   \item 2D PCA plots using mixOmics' \code{plotIndiv} function,
#'   \item 3D scatter plots (if \code{style} is "3d" or "3D" and
#'   \code{comp_num} is 3) via the plot3D package,
#'   \item Scree plots showing both individual and cumulative
#'   explained variance,
#'   \item Loadings plots, and
#'   \item Biplots and correlation circle plots.
#' }
#'
#' @return A PDF file containing the PCA plots is generated and saved when
#' \code{output_file} is provided. Otherwise, plots are displayed on the current
#' graphics device.
#' @author Shubh Saraswat
#'
#' @export
#' @importFrom mixOmics pca plotIndiv plotLoadings plotVar
#' @import ggplot2
#' @importFrom plot3D scatter3D
#'
#' @examples
#' # Load sample data
#' data <- ExampleData1[, -c(3,23)]
#' data_df <- dplyr::filter(data, Group != "ND" & Treatment != "Unstimulated")
#' # Run PCA analysis and save plots to a PDF file
#' cyt_pca(
#'   data = data_df,
#'   output_file = NULL,
#'   colors = c("black", "red2"),
#'   scale = "log2",
#'   comp_num = 3,
#'   pch_values = c(16, 4),
#'   style = "3D",
#'   group_col = "Group",
#'   group_col2 = "Treatment",
#'   ellipse = FALSE
#' )
#'

cyt_pca <- function(
  data,
  group_col = NULL,
  group_col2 = NULL,
  colors = NULL,
  output_file,
  ellipse = FALSE,
  comp_num = 2,
  scale = c("none", "log2", "log10", "zscore", "custom"),
  custom_fn = NULL,
  pch_values = NULL,
  style = NULL
) {
  names(data) <- make.names(names(data), unique = TRUE)
  # If one factor is missing, use the provided column for both grouping
  # and treatment.
  if (!is.null(group_col) && is.null(group_col2)) {
    message("No second grouping column provided; performing overall analysis.")
    group_col2 <- group_col
  }
  if (is.null(group_col) && !is.null(group_col2)) {
    stop(
      "No first grouping column provided; must provide the first grouping column."
    )
  }
  if (is.null(group_col) && is.null(group_col2)) {
    stop("At least one grouping column must be provided.")
  }
  # Identify numeric columns to transform (exclude grouping columns)
  scale <- match.arg(scale)
  id_cols <- if (is.null(group_col2)) group_col else c(group_col, group_col2)
  numeric_cols <- setdiff(names(data)[sapply(data, is.numeric)], id_cols)
  if (length(numeric_cols) > 0) {
    data <- apply_scale(
      data,
      columns = numeric_cols,
      scale = scale,
      custom_fn = custom_fn
    )
  }

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
  } else if (length(pch_values) < length(unique(group_vec))) {
    stop(
      "Please ensure the number of pch values provided (",
      length(pch_values),
      ") is at least equal to the number of unique groups (",
      length(unique(group_vec)),
      ") from the grouping column."
    )
  }
  # Generate a color palette if not provided (based on the
  # grouping variable levels)
  if (is.null(colors)) {
    num_groups <- length(unique(data[[group_col]]))
    colors <- rainbow(num_groups)
  }

  # Case 1: Overall PCA when both factors are the same.
  if (group_col == group_col2) {
    overall_analysis <- "Overall Analysis"

    # Remove the factor column(s) from predictors and keep only numeric columns
    the_data_df <- data[, !(names(data) %in% unique(c(group_col, group_col2)))]
    the_data_df <- the_data_df[, sapply(the_data_df, is.numeric)]

    the_groups <- as.vector(data[[group_col]])
    if (length(unique(the_groups)) < 2) {
      stop(
        "The grouping variable must have at least two levels for PCA.
           Please provide an appropriate grouping column."
      )
    }

    cytokine_pca <- mixOmics::pca(
      the_data_df,
      ncomp = comp_num,
      center = TRUE,
      scale = TRUE
    )

    group_factors <- seq_len(length(levels(factor(the_groups))))

    # Plot PCA individuals plot

    plot_args <- list(
      cytokine_pca,
      ind.names = NA,
      legend = TRUE,
      group = the_groups,
      col = colors,
      pch = pch_values[group_factors],
      title = paste("PCA:", overall_analysis),
      legend.title = group_col
    )
    if (ellipse) {
      plot_args$ellipse <- TRUE
    }
    overall_indiv_plot <- do.call(mixOmics::plotIndiv, plot_args)
    overall_3D <- NULL
    # 3D Plot if requested and exactly three components are used
    if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
      cytokine_scores <- cytokine_pca$variates$X
      overall_3D <- function() {
        plot3D::scatter3D(
          cytokine_scores[, 1],
          cytokine_scores[, 2],
          cytokine_scores[, 3],
          pch = pch_values,
          col = colors,
          xlab = "Component 1",
          ylab = "Component 2",
          zlab = "Component 3",
          main = paste("3D Plot:", overall_analysis),
          theta = 20,
          phi = 30,
          bty = "g",
          colkey = FALSE
        )
      }
      overall_3D()
    } else if (!is.null(style)) {
      stop(
        "Please enter a valid style for 3D plot: '3d' or '3D' or
           enter a valid number of components."
      )
    }

    # Scree Plot
    variances <- cytokine_pca$prop_expl_var$X
    cumulative_variances <- cytokine_pca$cum.var
    scree_data <- data.frame(
      Component = 1:comp_num,
      Variance = variances,
      Cumulative = cumulative_variances
    )

    scree_plot <- ggplot2::ggplot(scree_data, aes(x = Component)) +
      ggplot2::geom_line(
        aes(y = Variance, color = "Individual"),
        linewidth = 1
      ) +
      ggplot2::geom_point(aes(y = Variance, color = "Individual"), size = 2) +
      ggplot2::geom_line(
        aes(y = Cumulative, color = "Cumulative"),
        linewidth = 1,
        linetype = "dashed"
      ) +
      ggplot2::geom_point(aes(y = Cumulative, color = "Cumulative"), size = 2) +
      ggplot2::scale_color_manual(
        values = c(
          "Individual" = "blue",
          "Cumulative" = "green"
        )
      ) +
      ggplot2::labs(
        title = paste("Scree Plot:", overall_analysis),
        x = "Principal Components",
        y = "Explained Variance",
        color = "Variance Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::scale_x_continuous(breaks = 1:comp_num) +
      ggplot2::geom_text(
        aes(
          y = Variance,
          label = paste0(
            round(Variance * 100, 1),
            "%"
          )
        ),
        vjust = -1.5,
        hjust = 0.5,
        size = 4
      ) +
      ggplot2::geom_text(
        aes(
          y = Cumulative,
          label = paste0(
            round(Cumulative * 100, 1),
            "%"
          )
        ),
        vjust = 1.5,
        hjust = 0.5,
        size = 4
      )

    invisible(print(scree_plot))

    # Plot loadings for each component
    loadings_list <- setNames(
      lapply(seq_len(comp_num), function(comp) {
        force(comp) # capture comp in the closure
        function() {
          mixOmics::plotLoadings(
            cytokine_pca,
            comp = comp,
            size.names = 1,
            size.legend = 1,
            col = "grey",
            legend.color = colors,
            title = paste(
              "Loadings Component",
              comp,
              ":",
              overall_analysis
            ),
            legend = TRUE
          )
        }
      }),
      nm = paste0("Comp", seq_len(comp_num))
    )

    invisible(lapply(loadings_list, function(plot_fn) plot_fn()))

    # Biplot for components 1 and 2
    biplot_obj <- biplot(
      cytokine_pca,
      comp = c(1, 2),
      group = the_groups,
      col = colors,
      var.arrow.col = "blue",
      var.arrow.size = 0.5,
      var.arrow.length = 0.2,
      var.names = TRUE,
      var.names.col = "blue",
      var.names.size = 3,
      ind.names = FALSE,
      legend = TRUE,
      legend.title = "Group"
    )

    invisible(print(biplot_obj))

    # Correlation circle plot
    corr_plot <- function() {
      mixOmics::plotVar(
        cytokine_pca,
        comp = c(1, 2),
        var.names = TRUE,
        cex = 4,
        col = "black",
        overlap = TRUE,
        title = paste("Correlation Circle Plot:", overall_analysis),
        style = "ggplot2"
      )
    }
    corr_plot()

    result_list <- list(
      overall_indiv_plot = overall_indiv_plot$graph,
      loadings = loadings_list,
      plot_3d = if (!is.null(overall_3D)) overall_3D else NULL,
      biplot = biplot_obj,
      correlation_circle = corr_plot,
      scree_plot = scree_plot
    )
    if (!is.null(output_file)) {
      plot_list <- Filter(
        Negate(is.null),
        c(
          list(overall_indiv_plot$graph),
          list(scree_plot),
          unname(loadings_list), # list of closures
          list(function() print(biplot_obj)),
          list(corr_plot), # closure
          if (!is.null(overall_3D)) list(overall_3D)
        )
      )
      cyt_export(
        plot_list,
        filename = tools::file_path_sans_ext(output_file),
        format = tolower(tools::file_ext(output_file)),
        width = 8.5,
        height = 8
      )
    }
    invisible(result_list)
  } else {
    # Case 2: When grouping and treatment columns differ
    levels_vec <- unique(data[[group_col2]])
    indiv_plots <- list()
    all_plots <- list()
    overall_3D <- NULL
    for (i in seq_along(levels_vec)) {
      current_level <- levels_vec[i]
      title_sub <- current_level
      condt <- data[[group_col2]] == current_level

      the_data_df <- data[
        condt,
        !(names(data) %in%
          unique(c(
            group_col,
            group_col2
          )))
      ]
      the_data_df <- the_data_df[, sapply(the_data_df, is.numeric)]
      the_groups <- as.vector(data[condt, group_col])

      if (length(unique(the_groups)) < 2) {
        stop(
          "The grouping variable must have at least two levels for PCA.
             Please provide an appropriate grouping column."
        )
      }

      cytokine_pca <- mixOmics::pca(
        the_data_df,
        ncomp = comp_num,
        center = TRUE,
        scale = TRUE
      )

      group_factors <- seq_len(length(levels(factor(the_groups))))

      # PCA individuals plot for current treatment level

      plot_args <- list(
        cytokine_pca,
        ind.names = NA,
        legend = TRUE,
        group = the_groups,
        col = colors,
        pch = pch_values[group_factors],
        title = paste("PCA:", title_sub),
        legend.title = group_col
      )
      if (ellipse) {
        plot_args$ellipse <- TRUE
      }
      overall_indiv_plot <- do.call(mixOmics::plotIndiv, plot_args)
      indiv_plots[[as.character(
        current_level
      )]] <- overall_indiv_plot$graph
      # 3D Plot if applicable
      if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
        cytokine_scores <- cytokine_pca$variates$X
        overall_3D <- function() {
          plot3D::scatter3D(
            cytokine_scores[, 1],
            cytokine_scores[, 2],
            cytokine_scores[, 3],
            pch = pch_values,
            col = colors,
            xlab = "Component 1",
            ylab = "Component 2",
            zlab = "Component 3",
            main = paste("3D Plot:", current_level),
            theta = 20,
            phi = 30,
            bty = "g",
            colkey = FALSE
          )
        }
        overall_3D()
      } else if (!is.null(style)) {
        stop(
          "Please enter a valid style for 3D plot: '3d' or '3D' or
             enter a valid number of components."
        )
      }

      # Scree Plot for the current treatment level
      variances <- cytokine_pca$prop_expl_var$X
      cumulative_variances <- cytokine_pca$cum.var
      scree_data <- data.frame(
        Component = 1:comp_num,
        Variance = variances,
        Cumulative = cumulative_variances
      )

      scree_plot <- ggplot2::ggplot(scree_data, aes(x = Component)) +
        ggplot2::geom_line(
          aes(y = Variance, color = "Individual"),
          linewidth = 1
        ) +
        ggplot2::geom_point(aes(y = Variance, color = "Individual"), size = 2) +
        ggplot2::geom_line(
          aes(y = Cumulative, color = "Cumulative"),
          linewidth = 1,
          linetype = "dashed"
        ) +
        ggplot2::geom_point(
          aes(y = Cumulative, color = "Cumulative"),
          size = 2
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "Individual" = "blue",
            "Cumulative" = "green"
          )
        ) +
        ggplot2::labs(
          title = paste("Scree Plot:", current_level),
          x = "Principal Components",
          y = "Explained Variance",
          color = "Variance Type"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::scale_x_continuous(breaks = 1:comp_num) +
        ggplot2::geom_text(
          aes(
            y = Variance,
            label = paste0(
              round(Variance * 100, 1),
              "%"
            )
          ),
          vjust = -1.5,
          hjust = 0.5,
          size = 4
        ) +
        ggplot2::geom_text(
          aes(
            y = Cumulative,
            label = paste0(
              round(Cumulative * 100, 1),
              "%"
            )
          ),
          vjust = 1.5,
          hjust = 0.5,
          size = 4
        )

      invisible(print(scree_plot))

      # Plot loadings for each component
      loadings_list <- setNames(
        lapply(seq_len(comp_num), function(comp) {
          force(comp) # capture comp in the closure
          function() {
            mixOmics::plotLoadings(
              cytokine_pca,
              comp = comp,
              size.names = 1,
              size.legend = 1,
              col = "grey",
              legend.color = colors,
              title = paste(
                "Loadings Component",
                comp,
                ":",
                current_level
              ),
              legend = TRUE
            )
          }
        }),
        nm = paste0("Comp", seq_len(comp_num))
      )

      invisible(lapply(loadings_list, function(plot_fn) plot_fn()))

      # Biplot for components 1 and 2
      biplot_obj <- biplot(
        cytokine_pca,
        comp = c(1, 2),
        group = the_groups,
        col = colors,
        var.arrow.col = "blue",
        var.arrow.size = 0.5,
        var.arrow.length = 0.2,
        var.names = TRUE,
        var.names.col = "blue",
        var.names.size = 3,
        ind.names = FALSE,
        legend = TRUE,
        legend.title = "Group"
      )
      invisible(print(biplot_obj))

      # Correlation circle plot
      corr_plot <- function() {
        mixOmics::plotVar(
          cytokine_pca,
          comp = c(1, 2),
          var.names = TRUE,
          cex = 4,
          col = "black",
          overlap = TRUE,
          title = paste("Correlation Circle Plot:", current_level),
          style = "ggplot2"
        )
      }
      corr_plot()
      # After overall_indiv_plot is created:
      all_plots <- c(all_plots, list(overall_indiv_plot$graph))
      # After scree_plot:
      all_plots <- c(all_plots, list(scree_plot))
      # After loadings_list is created:
      all_plots <- c(all_plots, unname(loadings_list))
      # After biplot_obj:
      all_plots <- c(all_plots, list(function() print(biplot_obj)))
      # After corr_plot:
      all_plots <- c(all_plots, list(corr_plot))
      # After overall_3D (inside the 3D if block):
      if (!is.null(overall_3D)) {
        all_plots <- c(all_plots, list(overall_3D))
      }
    }
    result_list <- list(
      overall_indiv_plot = overall_indiv_plot,
      all_indiv_plots = indiv_plots,
      loadings = loadings_list,
      plot_3d = if (!is.null(overall_3D)) overall_3D else NULL,
      biplot = biplot_obj,
      correlation_circle = corr_plot,
      scree_plot = scree_plot
    )
    # AFTER the for loop, before return():
    if (!is.null(output_file)) {
      cyt_export(
        Filter(Negate(is.null), all_plots),
        filename = tools::file_path_sans_ext(output_file),
        format = tolower(tools::file_ext(output_file)),
        width = 8.5,
        height = 8
      )
    }
    invisible(result_list)
  }
}
