#' Analyze data with Multivariate INTegration (MINT) Sparse Partial Least Squares Discriminant Analysis (sPLS-DA).
#'
#' This function performs a MINT (Multivariate INTegrative) sPLS-DA to handle
#' batch effects by modeling a global biological signal across different studies or batches.
#' If a second grouping column (`group_col2`) is provided, the analysis is stratified
#' and performed for each level of that column.
#'
#' @param data A matrix or data frame containing the variables. Columns not
#'   specified by \code{group_col}, \code{group_col2}, or \code{multilevel_col} are assumed to be continuous
#'   variables for analysis.
#' @param group_col A string specifying the first grouping column name that contains grouping
#'   information. If \code{group_col2} is not provided, it will be used for both
#'   grouping and treatment.
#' @param group_col2 A string specifying the second grouping column name. Default is
#'   \code{NULL}.
#' @param batch_col A string specifying the column that identifies the batch or study for each sample.
#' @param colors A vector of splsda_colors for the groups or treatments. If
#'   \code{NULL}, a random palette (using \code{rainbow}) is generated based on
#'   the number of groups.
#' @param pdf_title A string specifying the file name for saving the PDF output.
#'                If set to NULL, the function runs in IDE plots pane.
#' @param ellipse Logical. Whether to draw a 95\% confidence ellipse on the figures.
#'   Default is \code{FALSE}.
#' @param bg Logical. Whether to draw the prediction background in the figures.
#'   Default is \code{FALSE}.
#' @param var_num Numeric. The number of variables to be used in the PLS-DA model.
#' @param scale Character. Option for data transformation; if set to \code{"log2"}, a log2
#'   transformation is applied to the continuous variables. Default is \code{NULL}.
#' @param comp_num Numeric. The number of components to calculate in the sPLS-DA model. Default is 2.
#' @param cim Logical. Whether to compute and plot the Clustered Image Map (CIM) heatmap. Default is \code{FALSE}.
#' @param roc Logical. Whether to compute and plot the ROC curve for the model.
#'   Default is \code{FALSE}.
#' @param verbose A logical value indicating whether to print additional
#'   informational output to the console. When \code{TRUE}, the function will
#'   display progress messages, and intermediate results when
#'   \code{FALSE} (the default), it runs quietly.
#' @return Plots consisting of the classification figures, ROC curves, correlation circle plots, and heatmaps.
#' @details
#' When \code{verbose} is set to \code{TRUE}, additional information about the analysis and confusion matrices
#' are printed to the console. These can be suppressed by keeping \code{verbose = FALSE}.
#'
#' @author Shubh Saraswat
#' @references Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017).
#' MINT: A multivariate integrative approach to identify a reproducible
#' biomarker signature across multiple experiments and platforms. BMC
#' Bioinformatics 18:128.
#' 
#' @examples
#' # Loading ExampleData5 dataset with batch column
#' data_df <- ExampleData5[,-c(2,4)]
#' data_df <- dplyr::filter(data_df, Group != "ND")
#'
#' cyt_mint_splsda(data_df, group_col = "Group",
#'  batch_col = "Batch", colors = c("black", "purple"),
#'  ellipse = TRUE, var_num = 25, comp_num = 2,
#'  scale = "log2", verbose = FALSE)
#'
#' @export
#' @importFrom mixOmics mint.splsda auroc plotIndiv plotLoadings plotVar cim background.predict
#' @import ggplot2
#' @importFrom grDevices rainbow pdf dev.off
cyt_mint_splsda <- function(
  data,
  group_col,
  batch_col,
  group_col2 = NULL,
  colors = NULL,
  pdf_title = NULL,
  ellipse = TRUE,
  bg = FALSE,
  var_num = 20,
  comp_num = 2,
  scale = NULL,
  cim = FALSE,
  roc = FALSE,
  verbose = FALSE
) {
  # If one factor is missing, use the provided column for
  # both grouping and treatment.
  if (!is.null(group_col) && is.null(group_col2)) {
    if (verbose) {
      cat("No second grouping column provided; performing overall analysis.\n")
    }
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
  # Identify columns given by user
  id_cols <- c(group_col, group_col2, batch_col)
  id_cols <- id_cols[!is.na(id_cols) & id_cols %in% names(data)]

  # Optionally apply log2 transformation
  if (!is.null(scale) && scale == "log2") {
    data <- data.frame(
      data[, id_cols, drop = FALSE],
      log2(data[, !(names(data) %in% id_cols), drop = FALSE])
    )
    if (verbose) message("Applied log2 transformation.\n")
  } else if (is.null(scale) && verbose) {
    if (verbose) message("No data transformation applied.\n")
  }

  if (!is.null(batch_col)) {
    if (!(batch_col %in% names(data))) {
      stop(sprintf("Batch column '%s' not found in your data.", batch_col))
    }
  } else if (is.null(batch_col)) {
    stop("Batch column must be provided.")
  }

  # Generate a color palette if not provided (based on the
  # grouping variable levels in the entire dataset)
  if (is.null(colors)) {
    num_groups <- length(unique(data[[group_col]]))
    colors <- rainbow(num_groups)
  }

  if (!is.null(pdf_title)) {
    grDevices::pdf(file = pdf_title, width = 8.5, height = 8)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  # Perform MINT sPLS-DA
  # Case 1: If group_col == group_col2, perform overall analysis

  if (group_col == group_col2) {
    overall_analysis <- "Overall Analysis"
    # --- 1. Data Preparation ---
    Y <- data[[group_col]]
    study <- data[[batch_col]]
    X <- data[,
      !(names(data) %in% c(group_col, group_col2, batch_col)),
      drop = FALSE
    ]
    X <- X[, sapply(X, is.numeric)]

    # --- 2. Run MINT sPLS-DA Model ---
    final_model <- mixOmics::mint.splsda(
      X = X,
      Y = Y,
      study = study,
      ncomp = comp_num,
      keepX = rep(var_num, comp_num),
      scale = TRUE
    )

    # --- 3. Calculate Prediction Accuracy ---
    mint_predict <- stats::predict(
      final_model,
      X,
      study.test = study,
      dist = 'max.dist'
    )
    # Get predictions from the final component for the best accuracy
    final_predictions <- mint_predict$class$max.dist[, comp_num]
    accuracy <- sum(Y == final_predictions) / length(Y)
    acc_percent <- signif(accuracy * 100, 2)

    # --- 4. Prepare plot titles and background object ---
    title_label <- if (nzchar(overall_analysis)) {
      paste(overall_analysis)
    } else {
      "MINT sPLS-DA Global Plot"
    }
    bg_obj <- if (bg) {
      try(
        mixOmics::background.predict(
          final_model,
          comp.predicted = 2,
          dist = "max.dist"
        ),
        silent = TRUE
      )
    } else {
      NULL
    }
    # -- Global individual plot --
    plot_args <- list(
      final_model,
      study = "global",
      group = Y,
      col = colors,
      legend = TRUE,
      legend.title = group_col,
      subtitle = paste0(title_label, " - Accuracy: ", acc_percent, "%")
    )
    if (ellipse) {
      plot_args$ellipse <- ellipse
    }
    if (!is.null(bg_obj)) {
      plot_args$background <- bg_obj
    }
    global_indiv_plot <- do.call(mixOmics::plotIndiv, plot_args)

    # -- Global loadings plots --
    global_loadings_plots <- setNames(
      lapply(seq_len(comp_num), function(comp) {
        force(comp)
        function() {
          mixOmics::plotLoadings(
            final_model,
            comp = comp,
            study = "global",
            contrib = "max",
            method = "mean",
            size.name = 1,
            size.legend = 1,
            size.title = 1,
            legend.color = colors,
            legend.title = group_col,
            title = paste("Global Loadings - Comp", comp)
          )
        }
      }),
      nm = paste0("Comp", seq_len(comp_num))
    )
    invisible(lapply(global_loadings_plots, function(plot_fn) plot_fn()))

    plot_args2 <- list(
      final_model,
      study = "all.partial",
      group = Y,
      col = colors,
      legend = TRUE,
      legend.title = group_col,
      title = paste("Partial Plots:", overall_analysis)
    )
    if (!is.null(bg_obj)) {
      plot_args2$background <- bg_obj
    }
    partial_indiv_plot <- do.call(
      mixOmics::plotIndiv,
      plot_args2
    )

    partial_loadings_plots <- setNames(
      lapply(seq_len(comp_num), function(comp) {
        force(comp)
        function() {
          mixOmics::plotLoadings(
            final_model,
            comp = comp,
            study = "all.partial",
            contrib = "max",
            method = "mean",
            legend.color = colors,
            legend.title = group_col,
            title = paste("Partial Loadings - Comp", comp)
          )
        }
      }),
      nm = paste0("Comp", seq_len(comp_num))
    )
    invisible(lapply(partial_loadings_plots, function(plot_fn) plot_fn()))

    # -- CIM heatmap --
    cim_plot <- NULL
    if (cim) {
      cim_plot <- function() {
        for (comp in seq_len(comp_num)) {
          mixOmics::cim(
            final_model,
            comp = comp,
            row.sideColors = colors[as.numeric(Y)],
            row.names = FALSE,
            title = paste0(
              "Clustered Image Map (CIM) for Comp: ",
              comp,
              " ",
              overall_analysis
            )
          )
        }
      }
      cim_plot()
    }

    # -- Correlation circle --
    correlation_circle_plot <- NULL
    correlation_circle_plot <- function() {
      mixOmics::plotVar(
        final_model,
        comp = c(1, 2),
        cex = 4,
        col = "black",
        overlap = TRUE,
        style = "ggplot2",
        var.names = TRUE,
        legend = TRUE,
        title = paste("Correlation Circle:", overall_analysis)
      )
    }
    correlation_circle_plot()

    # -- ROC curve, if requested --
    roc_plot <- NULL
    if (roc) {
      roc_plot <- mixOmics::auroc(
        object = final_model,
        newdata = X,
        outcome.test = Y,
        study.test = study,
        roc.study = "global",
        plot = TRUE,
        print = FALSE
      )
    }
    results_list <- list(
      global_indiv_plot = global_indiv_plot$graph,
      global_loadings_plots = global_loadings_plots,
      partial_indiv_plot = partial_indiv_plot$graph,
      partial_loadings_plots = partial_loadings_plots,
      cim_plot = cim_plot,
      correlation_circle_plot = correlation_circle_plot,
      roc_plot = roc_plot
    )
    return(invisible(results_list))
  } else {
    # Case 2: If group_col != group_col2, perform stratified analysis
    levels_vec <- unique(data[[group_col2]])
    indiv_plots <- list()
    results_list <- list()
    for (lvl in levels_vec) {
      # 1) subset to only those rows where group_col2 == lvl
      df_sub <- data[data[[group_col2]] == lvl, , drop = FALSE]

      # 2) build Y, study, and X *from the subset*
      Y <- df_sub[[group_col]]
      study <- df_sub[[batch_col]]
      X <- df_sub[,
        !(names(df_sub) %in% c(group_col, group_col2, batch_col)),
        drop = FALSE
      ]
      # only the numeric predictors:
      X <- X[, sapply(X, is.numeric), drop = FALSE]

      # --- 2. Run MINT sPLS-DA Model ---
      final_model <- mixOmics::mint.splsda(
        X = X,
        Y = Y,
        study = study,
        ncomp = comp_num,
        keepX = rep(var_num, comp_num),
        scale = TRUE
      )

      # --- 3. Calculate Prediction Accuracy ---
      mint_predict <- stats::predict(
        final_model,
        X,
        study.test = study,
        dist = 'max.dist'
      )
      # Get predictions from the final component for the best accuracy
      final_predictions <- mint_predict$class$max.dist[, comp_num]
      accuracy <- sum(Y == final_predictions) / length(Y)
      acc_percent <- signif(accuracy * 100, 2)

      # --- 4. Prepare plot titles and background object ---
      title_label <- if (nzchar(lvl)) {
        paste(lvl)
      } else {
        "MINT sPLS-DA Global Plot"
      }
      bg_obj <- if (bg) {
        try(
          mixOmics::background.predict(
            final_model,
            comp.predicted = 2,
            dist = "max.dist"
          ),
          silent = TRUE
        )
      } else {
        NULL
      }
      # -- Global individual plot --
      plot_args <- list(
        final_model,
        study = "global",
        group = Y,
        col = colors,
        legend = TRUE,
        legend.title = group_col,
        subtitle = paste0(title_label, " - Accuracy: ", acc_percent, "%")
      )
      if (ellipse) {
        plot_args$ellipse <- ellipse
      }
      if (!is.null(bg_obj)) {
        plot_args$background <- bg_obj
      }
      global_indiv_plot <- do.call(mixOmics::plotIndiv, plot_args)
      indiv_plots[[as.character(lvl)]] <- global_indiv_plot$graph

      # -- Global loadings plots --
      global_loadings_plots <- setNames(
        lapply(seq_len(comp_num), function(comp) {
          force(comp)
          function() {
            mixOmics::plotLoadings(
              final_model,
              comp = comp,
              study = "global",
              contrib = "max",
              method = "mean",
              size.name = 1,
              size.legend = 1,
              size.title = 1,
              legend.color = colors,
              legend.title = group_col,
              title = paste("Global Loadings - Comp", comp)
            )
          }
        }),
        nm = paste0("Comp", seq_len(comp_num))
      )
      invisible(lapply(global_loadings_plots, function(plot_fn) plot_fn()))

      plot_args2 <- list(
        final_model,
        study = "all.partial",
        group = Y,
        col = colors,
        legend = TRUE,
        legend.title = group_col,
        title = paste("Partial Plots:", lvl)
      )
      if (!is.null(bg_obj)) {
        plot_args2$background <- bg_obj
      }
      partial_indiv_plot <- do.call(
        mixOmics::plotIndiv,
        plot_args2
      )

      partial_loadings_plots <- setNames(
        lapply(seq_len(comp_num), function(comp) {
          force(comp)
          function() {
            mixOmics::plotLoadings(
              final_model,
              comp = comp,
              study = "all.partial",
              contrib = "max",
              method = "mean",
              legend.color = colors,
              legend.title = group_col,
              title = paste("Partial Loadings - Comp", comp)
            )
          }
        }),
        nm = paste0("Comp", seq_len(comp_num))
      )
      invisible(lapply(partial_loadings_plots, function(plot_fn) plot_fn()))

      # -- CIM heatmap --
      cim_plot <- NULL
      if (cim) {
        cim_plot <- function() {
          for (comp in seq_len(comp_num)) {
            mixOmics::cim(
              final_model,
              comp = comp,
              row.sideColors = colors[as.numeric(Y)],
              row.names = FALSE,
              title = paste0(
                "Clustered Image Map (CIM) for Comp: ",
                comp,
                " ",
                overall_analysis
              )
            )
          }
        }
        cim_plot()
      }

      # -- Correlation circle --
      correlation_circle_plot <- NULL
      correlation_circle_plot <- function() {
        mixOmics::plotVar(
          final_model,
          comp = c(1, 2),
          cex = 4,
          col = "black",
          overlap = TRUE,
          style = "ggplot2",
          var.names = TRUE,
          legend = TRUE,
          title = paste("Correlation Circle:", lvl)
        )
      }
      correlation_circle_plot()

      # -- ROC curve, if requested --
      roc_plot <- NULL
      if (roc) {
        roc_plot <- mixOmics::auroc(
          object = final_model,
          newdata = X,
          outcome.test = Y,
          study.test = study,
          roc.study = "global",
          plot = TRUE,
          print = FALSE
        )
      }
      results_list[[lvl]] <- list(
        global_indiv_plot = global_indiv_plot$graph,
        all_indiv_plots = indiv_plots,
        global_loadings_plots = global_loadings_plots,
        partial_indiv_plot = partial_indiv_plot$graph,
        partial_loadings_plots = partial_loadings_plots,
        cim_plot = cim_plot,
        correlation_circle_plot = correlation_circle_plot,
        roc_plot = roc_plot
      )
    }
    return(invisible(results_list))
  }
}
