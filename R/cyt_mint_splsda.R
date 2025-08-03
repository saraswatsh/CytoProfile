#' Analyze data with MINT Sparse Partial Least Squares Discriminant Analysis (sPLS-DA).
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
#' @param comp_num Numeric. The number of components to calculate in the sPLS-DA model.
#'   Default is 2.
#'   to be used in the plots.
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
#' @examples
#' pdf(file = tempfile(fileext = ".pdf"))
#' # Loading ExampleData5 dataset with batch column
#' data_df <- ExampleData5[,-c(2,4)]
#' data_df <- dplyr::filter(data_df, Group != "ND")
#'
#' results <- cyt_mint_splsda(data_df, group_col = "Group",
#'  batch_col = "Batch", colors = c("black", "purple"),
#'  ellipse = TRUE, var_num = 25, comp_num = 2,
#'  scale = "log2", verbose = FALSE)
#'
#' dev.off()
#'
#' @export
#' @importFrom mixOmics mint.splsda auroc plotIndiv plotLoadings plotVar cim background.predict
#' @import ggplot2
#' @importFrom grDevices rainbow pdf dev.off recordPlot replayPlot
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
  roc = FALSE,
  verbose = FALSE
) {
  # --- Helper function to run the core analysis ---
  run_mint_analysis <- function(
    data_subset,
    analysis_label = "",
    is_pdf_mode = !is.null(pdf_title)
  ) {
    # --- 1. Data Preparation ---
    data_subset[[group_col]] <- as.factor(data_subset[[group_col]])
    data_subset[[batch_col]] <- as.factor(data_subset[[batch_col]])
    if (nlevels(data_subset[[group_col]]) < 2) {
      if (verbose) {
        message(paste(
          "Skipping '",
          analysis_label,
          "': requires at least two group levels.",
          sep = ""
        ))
      }
      return(NULL)
    }
    if (nlevels(data_subset[[batch_col]]) < 2) {
      if (verbose) {
        message(paste(
          "Skipping '",
          analysis_label,
          "': MINT requires at least two batches.",
          sep = ""
        ))
      }
      return(NULL)
    }
    Y <- data_subset[[group_col]]
    study <- data_subset[[batch_col]]
    X <- data_subset[,
      !(names(data_subset) %in% c(group_col, group_col2, batch_col)),
      drop = FALSE
    ]
    X <- X[, sapply(X, is.numeric)]

    # --- 2. Run MINT sPLS-DA Model ---
    final_model <- mixOmics::mint.splsda(
      X = X,
      Y = Y,
      study = study,
      ncomp = comp_num,
      keepX = rep(var_num, comp_num)
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
    title_label <- if (nzchar(analysis_label)) {
      paste("MINT sPLS-DA:", analysis_label)
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

    record_base_plot <- function(expr) {
      # enable recording on the current device
      grDevices::dev.control(displaylist = "enable")
      expr # draw the plot
      grDevices::recordPlot() # capture and return it
    }
    do_partial <- !is.null(group_col2) && group_col2 != group_col

    # --- 5. Handle Output: PDF vs. Direct Call ---
    if (is_pdf_mode) {
      # ---- strictly global plots ----
      if (!do_partial) {
        mixOmics::plotIndiv(
          final_model,
          study = "global",
          group = Y,
          col = colors,
          legend = TRUE,
          legend.title = group_col,
          subtitle = paste0(title_label, " - Accuracy: ", acc_percent, "%"),
          ellipse = ellipse,
          background = bg_obj
        )
        # global loadings:
        for (i in seq_len(comp_num)) {
          mixOmics::plotLoadings(
            final_model,
            comp = i,
            legend.color = colors,
            legend.title = group_col,
            study = "global",
            contrib = "max",
            method = "mean",
            title = paste("Global Loadings - Comp", i)
          )
        }
        # CIM heatmap and correlation circle:
        mixOmics::cim(
          final_model,
          comp = 1,
          row.sideColors = colors[as.numeric(Y)],
          row.names = FALSE,
          title = paste("CIM (Comp 1):", analysis_label)
        )
        mixOmics::plotVar(
          final_model,
          var.names = TRUE,
          legend = TRUE,
          title = paste("Correlation Circle:", analysis_label)
        )
        if (roc) {
          try(
            mixOmics::auroc(
              object = final_model,
              newdata = X,
              outcome.test = Y,
              study.test = study,
              roc.study = "global",
              plot = TRUE,
              print = FALSE
            ),
            silent = TRUE
          )
        }
        return(NULL)
      }

      # ---- stratified (partial) branch ----
      mixOmics::plotIndiv(
        final_model,
        study = "global",
        group = Y,
        col = colors,
        legend = TRUE,
        legend.title = group_col,
        subtitle = paste0(title_label, " - Accuracy: ", acc_percent, "%"),
        ellipse = ellipse,
        background = bg_obj
      )
      # global loadings
      for (i in seq_len(comp_num)) {
        mixOmics::plotLoadings(
          final_model,
          comp = i,
          legend.color = colors,
          legend.title = group_col,
          study = "global",
          contrib = "max",
          method = "mean",
          title = paste("Global Loadings - Comp", i)
        )
      }
      # partial loadings
      for (i in seq_len(comp_num)) {
        mixOmics::plotLoadings(
          final_model,
          comp = i,
          legend.color = colors,
          legend.title = group_col,
          study = "all.partial",
          contrib = "max",
          method = "mean",
          title = paste("Partial Loadings - Comp", i)
        )
      }
      # partial indiv
      mixOmics::plotIndiv(
        final_model,
        study = "all.partial",
        group = Y,
        col = colors,
        legend = TRUE,
        title = paste("Partial Plots:", analysis_label)
      )
      # CIM heatmap and correlation circle
      mixOmics::cim(
        final_model,
        comp = 1,
        row.sideColors = colors[as.numeric(Y)],
        row.names = FALSE,
        title = paste("CIM (Comp 1):", analysis_label)
      )
      mixOmics::plotVar(
        final_model,
        var.names = TRUE,
        legend = TRUE,
        title = paste("Correlation Circle:", analysis_label)
      )
      if (roc) {
        try(
          mixOmics::auroc(
            object = final_model,
            newdata = X,
            outcome.test = Y,
            study.test = study,
            roc.study = "global",
            plot = TRUE,
            print = FALSE
          ),
          silent = TRUE
        )
      }
      return(NULL)
    }

    ## ----------------------------------------------------------------------------
    ## 6. list-of-recorded-plots Mode
    ## ----------------------------------------------------------------------------
    results_list <- list()

    # --- global indiv ---
    mixOmics::plotIndiv(
      final_model,
      study = "global",
      group = Y,
      col = colors,
      legend = TRUE,
      legend.title = group_col,
      subtitle = paste0(title_label, " - Accuracy: ", acc_percent, "%"),
      ellipse = ellipse,
      background = bg_obj
    )
    results_list$global_indiv_plot <- grDevices::recordPlot()

    # --- global loadings ---
    results_list$global_loadings_plots <- lapply(
      seq_len(comp_num),
      function(i) {
        mixOmics::plotLoadings(
          final_model,
          comp = i,
          legend.color = colors,
          legend.title = group_col,
          study = "global",
          contrib = "max",
          method = "mean",
          title = paste("Global Loadings - Comp", i)
        )
        grDevices::recordPlot()
      }
    )

    # --- partial branch? ---
    if (do_partial) {
      # partial loadings
      results_list$partial_loadings_plots <- lapply(
        seq_len(comp_num),
        function(i) {
          mixOmics::plotLoadings(
            final_model,
            comp = i,
            legend.color = colors,
            legend.title = group_col,
            study = "all.partial",
            contrib = "max",
            method = "mean",
            title = paste("Partial Loadings - Comp", i)
          )
          grDevices::recordPlot()
        }
      )
      # partial indiv
      mixOmics::plotIndiv(
        final_model,
        study = "all.partial",
        group = Y,
        col = colors,
        legend = TRUE,
        title = paste("Partial Plots:", analysis_label)
      )
      results_list$partial_indiv_plot <- grDevices::recordPlot()
    }

    # --- CIM heatmap & correlation circle ---
    mixOmics::cim(
      final_model,
      comp = 1,
      row.sideColors = colors[as.numeric(Y)],
      row.names = FALSE,
      title = paste("CIM (Comp 1):", analysis_label)
    )
    results_list$cim_plot <- grDevices::recordPlot()

    mixOmics::plotVar(
      final_model,
      var.names = TRUE,
      legend = TRUE,
      title = paste("Correlation Circle:", analysis_label)
    )
    results_list$correlation_circle_plot <- grDevices::recordPlot()

    # --- ROC ---
    if (roc) {
      results_list$roc_plot <- mixOmics::auroc(
        object = final_model,
        newdata = X,
        outcome.test = Y,
        study.test = study,
        roc.study = "global",
        plot = TRUE,
        print = FALSE
      )
    }

    return(results_list)
  }

  # --- Main Execution ---
  if (!is.null(scale) && scale == "log2") {
    id_cols <- unique(c(group_col, group_col2, batch_col))
    id_cols <- id_cols[!sapply(id_cols, is.null)]
    numeric_cols <- data[, !(names(data) %in% id_cols), drop = FALSE]
    data <- data.frame(data[, id_cols, drop = FALSE], log2(numeric_cols))
  }
  num_groups <- nlevels(as.factor(data[[group_col]]))
  if (is.null(colors) || length(colors) < num_groups) {
    colors <- grDevices::rainbow(num_groups)
  }
  is_pdf <- !is.null(pdf_title)
  if (is_pdf) {
    grDevices::pdf(file = pdf_title, width = 11, height = 8.5)
    on.exit(grDevices::dev.off())
  }
  if (is.null(group_col2) || group_col == group_col2) {
    results <- run_mint_analysis(
      data,
      analysis_label = "Overall",
      is_pdf_mode = is_pdf
    )
  } else {
    treatments <- levels(as.factor(data[[group_col2]]))
    if (is_pdf) {
      for (trt in treatments) {
        run_mint_analysis(
          data[data[[group_col2]] == trt, , drop = FALSE],
          trt,
          TRUE
        )
      }
      results <- paste("Output file generated:", normalizePath(pdf_title))
    } else {
      results <- lapply(treatments, function(trt) {
        run_mint_analysis(
          data[data[[group_col2]] == trt, , drop = FALSE],
          trt,
          FALSE
        )
      })
      names(results) <- treatments
    }
  }
  invisible(results)
}
