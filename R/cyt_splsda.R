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
#' @param multilevel_col A string specifying the column name that identifies
#'   repeated measurements (e.g., patient or sample IDs). If provided, a
#'   multilevel analysis will be performed. Default is \code{NULL}.
#' @param batch_col A string specifying the column that identifies the batch or study for each sample.
#' @param ind_names  If \code{TRUE}, the row names of the first (or second) data matrix is used as names.
#'   Default is \code{FALSE}. If a character vector is provided, these values will be used as names.
#'   If 'pch' is set this will overwrite the names as shapes. See ?mixOmics::plotIndiv for details.
#' @param colors A vector of colors for the groups or treatments. If
#'   \code{NULL}, a random palette (using \code{rainbow}) is generated based on
#'   the number of groups.
#' @param output_file Optional string specifying the name of the file
#'   to be created.  When `NULL` (default), plots are drawn on
#'   the current graphics device. Ensure that the file
#'   extension matches the desired format (e.g., ".pdf" for PDF output
#'   or ".png" for PNG output or .tiff for TIFF output).
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
#' @param comp_num Numeric. The number of components to calculate in the sPLS-DA
#'   model. Default is 2.
#' @param pch_values A vector of integers specifying the plotting characters
#'   (pch values) to be used in the plots.
#' @param style Character. If set to \code{"3D"} or \code{"3d"} and
#'   \code{comp_num} equals 3, a 3D plot is generated using the
#'   \code{plot3D} package. Default is \code{NULL}.
#' @param roc Logical. Whether to compute and plot the ROC curve for the model.
#'   Default is \code{FALSE}.
#' @param verbose A logical value indicating whether to print additional
#'   informational output to the console. When \code{TRUE}, the function will
#'   display progress messages, and intermediate results when
#'   \code{FALSE} (the default), it runs quietly.
#' @param seed An integer specifying the seed for reproducibility (default is 123).
#' @param tune Logical.  If `TRUE`, performs tuning of `ncomp` and
#'   `keepX` via cross‑validation.  Default is `FALSE`.
#' @param tune_folds Integer.  Number of folds in cross‑validation when
#'   tuning.  Default is 5.
#' @param scale Character string specifying a transformation to apply to the
#'   numeric predictor columns prior to model fitting.  Options are
#'   "none", "log2", "log10", "zscore", or "custom".  When
#'   "custom" is selected a user defined function must be supplied via
#'   `custom_fn`.  Defaults to "none".
#' @param custom_fn A custom transformation function used when
#'   `scale = "custom"`.  Ignored otherwise.  It should take a numeric
#'   vector and return a numeric vector of the same length.
#' @description
#' This function conducts Sparse Partial Least Squares Discriminant Analysis
#' (sPLS-DA) on the provided data. It uses the specified \code{group_col} (and
#' optionally \code{group_col2}) to define class labels while assuming the remaining
#' columns contain continuous variables. The function supports transformations
#' via the \code{scale} parameter and generates a series of plots,
#' including classification plots, scree plots, loadings plots, and VIP score
#' plots. Optionally, ROC curves are produced when \code{roc} is \code{TRUE}.
#' Additionally, cross-validation is supported via LOOCV or Mfold methods. When
#' both \code{group_col} and \code{group_col2} are provided and differ, the function
#' analyzes each treatment level separately.
#'
#' @return Plots consisting of the classification figures, component figures
#'   with Variable of Importance in Projection (VIP) scores, and classifications
#'   based on VIP scores greater than 1. ROC curves and confusion matrices are also
#'   produced if requested.
#' @details
#' When \code{verbose} is set to \code{TRUE}, additional information about the analysis and confusion matrices
#' are printed to the console. These can be suppressed by keeping \code{verbose = FALSE}.
#' @author Xiaohua Douglas Zhang and Shubh Saraswat
#'
#' @references Lê Cao, K.-A., Boitard, S. and Besse, P. (2011).
#' Sparse PLS Discriminant Analysis: biologically relevant feature selection
#' and graphical displays for multiclass problems. \emph{BMC Bioinformatics}
#' \bold{12}:253.
#' @examples
#' # Loading Sample Data
#' data_df <- ExampleData1[,-c(3)]
#' data_df <- dplyr::filter(data_df, Group != "ND", Treatment != "Unstimulated")
#'
#' cyt_splsda(data_df, output_file = NULL,
#' colors = c("black", "purple"), bg = FALSE, scale = "log2",
#' conf_mat = FALSE, var_num = 25, cv_opt = NULL, comp_num = 2,
#' pch_values = c(16, 4), style = NULL, ellipse = TRUE,
#' group_col = "Group", group_col2 = "Treatment", roc = FALSE, verbose = FALSE)
#'
#' @export
#' @importFrom mixOmics splsda background.predict perf vip auroc plotIndiv plotLoadings
#' @import ggplot2
#' @import dplyr
#' @importFrom plot3D scatter3D
#' @importFrom reshape2 melt
#' @importFrom caret confusionMatrix

cyt_splsda <- function(
  data,
  group_col = NULL,
  group_col2 = NULL,
  multilevel_col = NULL,
  batch_col = NULL,
  ind_names = FALSE,
  colors = NULL,
  output_file = NULL,
  ellipse = FALSE,
  bg = FALSE,
  conf_mat = FALSE,
  var_num,
  cv_opt = NULL,
  fold_num = 5,
  scale = c("none", "log2", "log10", "zscore", "custom"),
  custom_fn = NULL,
  tune = FALSE,
  tune_folds = 5,
  comp_num = 2,
  pch_values,
  style = NULL,
  roc = FALSE,
  verbose = FALSE,
  seed = 123
) {
  names(data) <- make.names(names(data), unique = TRUE)
  data <- as.data.frame(data)
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
  id_cols <- c(group_col, group_col2, multilevel_col, batch_col)
  id_cols <- id_cols[!is.na(id_cols) & id_cols %in% names(data)]

  if (!is.null(batch_col)) {
    if (!(batch_col %in% names(data))) {
      stop(sprintf("Batch column '%s' not found in your data.", batch_col))
    }
    if (verbose) {
      message(
        "Applying per-batch z-score normalization on numeric predictors.\n"
      )
    }
    id_cols2 <- unique(c(id_cols, batch_col))
    num_cols <- setdiff(
      names(data)[sapply(data, is.numeric)],
      id_cols2
    )

    data <- data |>
      dplyr::group_by(!!dplyr::sym(batch_col)) |>
      dplyr::mutate(
        dplyr::across(
          .cols = dplyr::all_of(num_cols),
          .fns = ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)
        )
      ) |>
      dplyr::ungroup()
  }
  scale <- match.arg(scale)
  num_cols <- setdiff(names(data)[sapply(data, is.numeric)], id_cols)
  if (length(num_cols) > 0) {
    data <- apply_scale(
      data,
      columns = num_cols,
      scale = scale,
      custom_fn = custom_fn
    )
  }
  # Now perform the check for pch_values:
  if (is.null(pch_values)) {
    stop("Please enter a vector of pch values, e.g. c(16, 4).")
  }
  if (group_col == group_col2) {
    if (length(pch_values) < length(unique(data[[group_col]]))) {
      stop(
        "Please ensure the number of pch values provided (",
        length(pch_values),
        ") is at least equal to the number of unique groups (",
        length(unique(data[[group_col]])),
        ") from the grouping column."
      )
    }
  } else {
    # When group_col and group_col2 differ, use the levels of group_col for pch
    if (length(pch_values) < length(unique(data[[group_col]]))) {
      stop(
        "Please ensure the number of pch values provided (",
        length(pch_values),
        ") is at least equal to the number of unique groups (",
        length(unique(data[[group_col]])),
        ") from the first grouping column."
      )
    }
  }

  # Generate a color palette if not provided (based on the
  # grouping variable levels in the entire dataset)
  if (is.null(colors)) {
    num_groups <- length(unique(data[[group_col]]))
    colors <- rainbow(num_groups)
  }

  # Case 1: Only one factor provided (both columns are the same)
  if (group_col == group_col2) {
    overall_analysis <- "Overall Analysis"

    # Remove the factor column from predictors and keep only numeric columns
    the_data_df <- data[, !(names(data) %in% c(group_col))]
    the_data_df <- the_data_df[, sapply(the_data_df, is.numeric)]

    the_groups <- as.vector(data[[group_col]])
    if (length(unique(the_groups)) < 2) {
      stop(
        "The grouping variable must have at least two levels for PLS-DA.
           Please provide an appropriate grouping column."
      )
    }
    multilevel_vec <- NULL
    if (!is.null(multilevel_col)) {
      if (!(multilevel_col %in% names(data))) {
        stop(sprintf("Multilevel column '%s' not found.", multilevel_col))
      }
      multilevel_vec <- data[[multilevel_col]]
      if (verbose) {
        message(sprintf("Using multilevel design: '%s'.\n"), multilevel_col)
      }
    }
    # Tuning sPLS-DA
    outcome <- factor(data[[group_col]])
    if (nlevels(outcome) < 2) {
      stop("Outcome must have at least two levels.")
    }
    predictors <- as.matrix(the_data_df)

    tune_res <- NULL
    if (tune) {
      tune_res <- mixOmics::tune.splsda(
        X = predictors,
        Y = outcome,
        ncomp = comp_num,
        test.keepX = c(5, 10, 15, 20, 25),
        validation = "Mfold",
        dist = "max.dist",
        nrepeat = 100,
        folds = tune_folds,
        multilevel = multilevel_vec,
        progressBar = FALSE
      )
      # Plot tuning results with proper title and labels
      p <- plot(tune_res) +
        ggplot2::ggtitle(paste("sPLS-DA Tuning Results:", overall_analysis)) +
        ggplot2::xlab("Number of Variables") +
        ggplot2::ylab("Balanced Error Rate (BER)")
      print(p)
      ncomp_final <- tune_res$choice.ncomp$ncomp
      keepX_final <- tune_res$choice.keepX[1:ncomp_final]
    }
    if (tune) {
      cytokine_splsda <- mixOmics::splsda(
        the_data_df,
        the_groups,
        scale = TRUE,
        ncomp = ncomp_final,
        keepX = keepX_final,
        multilevel = multilevel_vec
      )
    } else {
      cytokine_splsda <- mixOmics::splsda(
        the_data_df,
        the_groups,
        scale = TRUE,
        ncomp = comp_num,
        keepX = rep(var_num, comp_num),
        multilevel = multilevel_vec
      )
    }

    splsda_predict <- predict(cytokine_splsda, the_data_df, dist = "max.dist")
    prediction1 <- cbind(original = the_groups, splsda_predict$class$max.dist)
    accuracy1 <- (sum(prediction1[, 1] == prediction1[, 3]) /
      length(prediction1[, 1]))
    acc1 <- 100 * signif(accuracy1, digits = 2)

    if (bg) {
      bg_maxdist <- mixOmics::background.predict(
        cytokine_splsda,
        comp.predicted = 2,
        dist = "max.dist",
        xlim = c(-15, 15),
        ylim = c(-15, 15),
        resolution = 200
      )
    }
    group_factors <- seq_len(length(levels(factor(the_groups))))

    .build_indiv_args <- function(
      obj,
      colors,
      ind_names,
      pch_vec,
      title,
      legend_title,
      ellipse = FALSE,
      bg_obj = NULL,
      extra = list()
    ) {
      args <- c(
        list(
          obj,
          legend = TRUE,
          col = colors,
          title = title,
          legend.title = legend_title
        ),
        extra
      )

      if (ellipse) {
        args$ellipse <- TRUE
      }
      if (!is.null(bg_obj)) {
        args$background <- bg_obj
      }

      # labels OR shapes (never both)
      if (isTRUE(ind_names) || is.character(ind_names)) {
        args$ind.names <- ind_names
      } else {
        args$ind.names <- FALSE
        args$pch <- pch_vec
      }
      args
    }

    plot_args <- .build_indiv_args(
      obj = cytokine_splsda,
      colors = colors,
      ind_names = ind_names,
      pch_vec = pch_values, # overall case uses full vector
      title = paste(overall_analysis, "With Accuracy:", acc1, "%"),
      legend_title = group_col,
      ellipse = ellipse,
      bg_obj = if (bg) bg_maxdist else NULL
    )
    overall_indiv_plot <- do.call(mixOmics::plotIndiv, plot_args)
    overall_indiv_plots <- list(Overall = overall_indiv_plot$graph)

    overall_3D <- NULL
    if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
      cytokine_scores <- cytokine_splsda$variates$X
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
    }

    # If roc = TRUE, compute and plot ROC curve for the overall model
    overall_ROC <- NULL
    if (roc) {
      overall_ROC <- mixOmics::auroc(
        object = cytokine_splsda,
        newdata = the_data_df,
        outcome.test = the_groups,
        plot = TRUE,
        roc.comp = comp_num,
        title = paste0("ROC Curve:", overall_analysis),
        print = FALSE
      )
    }

    overall_CV <- NULL

    # Cross-validation methods
    if (!is.null(cv_opt)) {
      if (cv_opt == "loocv") {
        set.seed(seed)
        loocv_results <- mixOmics::perf(cytokine_splsda, validation = "loo")
        loocv_error_rate <- loocv_results$error.rate$overall[
          "comp2",
          "max.dist"
        ]
        loocv_acc <- 100 * signif(1 - loocv_error_rate, digits = 2)
        if (verbose) {
          cat("LOOCV Accuracy: ", paste0(loocv_acc, "%"), "\n")
        }

        error_rates <- loocv_results$error.rate$overall[, "max.dist"]
        error_df <- as.data.frame(error_rates)
        error_df$Component <- rownames(error_df)
        error_df <- reshape2::melt(
          error_df,
          id.vars = "Component",
          variable.name = "Distance",
          value.name = "ErrorRate"
        )
        overall_CV <- ggplot2::ggplot(
          error_df,
          ggplot2::aes(
            x = Component,
            y = ErrorRate,
            color = Distance,
            group = 1
          )
        ) +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 3) +
          ggplot2::labs(
            title = paste("LOOCV Error Rate:", overall_analysis),
            x = "Number of Components",
            y = "Error Rate"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
          ) +
          ggplot2::scale_color_manual(values = "red", labels = "max.dist")
        print(overall_CV)
      } else if (cv_opt == "Mfold") {
        set.seed(seed)
        fold_results <- mixOmics::perf(
          cytokine_splsda,
          validation = "Mfold",
          folds = fold_num,
          nrepeat = 100
        )
        fold_error_rate <- fold_results$error.rate$overall[
          "comp2",
          "max.dist"
        ]
        fold_acc <- 100 * signif(1 - fold_error_rate, digits = 2)
        if (verbose) {
          cat("Mfold Accuracy: ", paste0(fold_acc, "%"), "\n")
        }

        error_rates <- fold_results$error.rate$overall[, "max.dist"]
        error_df <- as.data.frame(error_rates)
        error_df$Component <- rownames(error_df)
        error_df <- reshape2::melt(
          error_df,
          id.vars = "Component",
          variable.name = "Distance",
          value.name = "ErrorRate"
        )
        overall_CV <- ggplot2::ggplot(
          error_df,
          ggplot2::aes(
            x = Component,
            y = ErrorRate,
            color = Distance,
            group = 1
          )
        ) +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 3) +
          ggplot2::labs(
            title = paste("Mfold Error Rate:", overall_analysis),
            x = "Number of Components",
            y = "Error Rate"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
          ) +
          ggplot2::scale_color_manual(values = "red", labels = "max.dist")
        print(overall_CV)
      }
    }

    # Loadings plot for each component
    loadings_list <- setNames(
      lapply(seq_len(comp_num), function(comp) {
        force(comp) # capture comp in the closure
        function() {
          mixOmics::plotLoadings(
            cytokine_splsda,
            comp = comp,
            contrib = "max",
            method = "mean",
            size.name = 1,
            size.legend = 1,
            legend.color = colors,
            title = paste(
              "Loadings for Component",
              comp,
              ":",
              overall_analysis
            ),
            size.title = 1,
            legend.title = group_col
          )
        }
      }),
      nm = paste0("Comp", seq_len(comp_num))
    )
    invisible(lapply(loadings_list, function(plot_fn) plot_fn()))

    # VIP scores and plot for PLS-DA with VIP > 1
    all_vip_scores <- mixOmics::vip(cytokine_splsda)
    vip_scores <- setNames(
      lapply(seq_len(comp_num), function(comp) {
        force(comp) # capture `comp` in the closure
        function() {
          # recreate  data frame
          vscore <- as.data.frame(all_vip_scores[, comp, drop = FALSE])
          vscore$metabo <- rownames(vscore)
          vscore$comp <- vscore[, 1]
          bar <- vscore[
            order(vscore$comp, decreasing = TRUE),
            c("metabo", "comp")
          ]
          bar <- bar[is.finite(bar$comp), , drop = FALSE]
          bar <- bar[bar$comp > 0, , drop = FALSE]

          # Match number with keepX
          if (tune) {
            bar <- head(bar, min(keepX_final, nrow(bar)))
          } else {
            bar <- head(bar, min(var_num, nrow(bar)))
          }
          # Lock the ordering to exactly what's being plotted (prevents gaps)
          bar$metabo <- factor(bar$metabo, levels = rev(bar$metabo))

          # build and print the plot
          p <- ggplot2::ggplot(bar, ggplot2::aes(x = metabo, y = comp)) +
            ggplot2::geom_bar(stat = "identity", position = "dodge") +
            ggplot2::scale_y_continuous(limits = c(0, max(bar$comp))) +
            ggplot2::geom_hline(yintercept = 1, color = "grey") +
            ggplot2::scale_x_discrete(limits = factor(bar$metabo)) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(
                angle = 45,
                hjust = 1,
                size = 15
              ),
              panel.grid = ggplot2::element_blank(),
              panel.background = ggplot2::element_rect(
                color = "black",
                fill = "transparent"
              )
            ) +
            ggplot2::coord_flip() +
            ggplot2::labs(x = "Variables", y = "VIP score") +
            ggplot2::ggtitle(paste("Component", comp))

          print(p)
        }
      }),
      nm = paste0("Comp", seq_len(comp_num))
    )
    invisible(lapply(vip_scores, function(draw_fn) draw_fn()))

    vip_indiv_plot <- NULL
    vip_loadings <- NULL
    vip_3D <- NULL
    vip_ROC <- NULL
    vip_CV <- NULL

    # PLS-DA on VIP > 1: Subset predictors with VIP > 1
    condt_variable <- all_vip_scores[, 1] > 1
    keep_x <- sum(condt_variable)
    the_data_mat <- the_data_df[, condt_variable, drop = FALSE]
    if (keep_x < 2) {
      if (verbose) {
        warning(
          "Only ",
          keep_x,
          " variable has VIP > 1. Skipping VIP > 1 PLS-DA Model."
        )
      }
      if (conf_mat == TRUE) {
        if (verbose) {
          cat(paste0(
            "Confusion Matrix for PLS-DA Comparison: ",
            overall_analysis,
            "\n"
          ))
        }

        # Confusion Matrix for main model
        cm_overall <- caret::confusionMatrix(
          data = as.factor(prediction1[, 3]), # predicted
          reference = as.factor(prediction1[, 1]) # actual
        )
        if (verbose) {
          print(cm_overall$table)
          cat("Accuracy:", signif(cm_overall$overall["Accuracy"], 2), "\n")

          # Check if binary or multi-class
          if (nlevels(as.factor(prediction1[, 1])) == 2) {
            cat(
              "Sensitivity:",
              signif(cm_overall$byClass["Sensitivity"], 2),
              "\n"
            )
            cat(
              "Specificity:",
              signif(cm_overall$byClass["Specificity"], 2),
              "\n"
            )
          } else {
            cat("\nPer-Class Sensitivity:\n")
            print(signif(cm_overall$byClass[, "Sensitivity"], 2))
            cat("\nPer-Class Specificity:\n")
            print(signif(cm_overall$byClass[, "Specificity"], 2))
            macro_sens <- mean(
              cm_overall$byClass[, "Sensitivity"],
              na.rm = TRUE
            )
            macro_spec <- mean(
              cm_overall$byClass[, "Specificity"],
              na.rm = TRUE
            )
            cat("\nMacro-Averaged Sensitivity:", signif(macro_sens, 2), "\n")
            cat("Macro-Averaged Specificity:", signif(macro_spec, 2), "\n")
          }
        }
      }
    } else {
      cytokine_splsda2 <- mixOmics::splsda(
        the_data_mat,
        the_groups,
        scale = TRUE,
        ncomp = comp_num,
        keepX = rep(keep_x, comp_num),
        multilevel = multilevel_vec
      )

      splsda_predict2 <- predict(
        cytokine_splsda2,
        the_data_mat,
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
      if (bg) {
        bg_maxdist2 <- mixOmics::background.predict(
          cytokine_splsda2,
          comp.predicted = 2,
          dist = "max.dist",
          xlim = c(-15, 15),
          ylim = c(-15, 15),
          resolution = 200
        )
      }
      plot_args2 <- .build_indiv_args(
        obj = cytokine_splsda2,
        colors = colors,
        ind_names = ind_names,
        pch_vec = pch_values, # overall case
        title = paste(overall_analysis, "(VIP>1)", "With Accuracy:", acc2, "%"),
        legend_title = group_col,
        ellipse = ellipse,
        bg_obj = if (bg) bg_maxdist2 else NULL
      )
      vip_indiv_plot <- do.call(mixOmics::plotIndiv, plot_args2)
      vip_indiv_plots <- list(VIP = vip_indiv_plot$graph)

      if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
        cytokine_scores2 <- cytokine_splsda2$variates$X
        vip_3D <- function() {
          plot3D::scatter3D(
            cytokine_scores2[, 1],
            cytokine_scores2[, 2],
            cytokine_scores2[, 3],
            pch = pch_values,
            col = colors,
            xlab = "Component 1",
            ylab = "Component 2",
            zlab = "Component 3",
            main = paste("3D Plot:", overall_analysis, "(VIP>1)"),
            theta = 20,
            phi = 30,
            bty = "g",
            colkey = FALSE
          )
        }
        vip_3D()
      }
      # Loadings plot for each component with VIP > 1
      vip_loadings <- setNames(
        lapply(seq_len(comp_num), function(comp) {
          force(comp) # capture comp in the closure
          function() {
            mixOmics::plotLoadings(
              cytokine_splsda2,
              comp = comp,
              contrib = "max",
              method = "mean",
              size.name = 1,
              size.legend = 1,
              legend.color = colors,
              title = paste(
                "Loadings for Component",
                comp,
                "(VIP > 1):",
                overall_analysis
              ),
              size.title = 1,
              legend.title = group_col
            )
          }
        }),
        nm = paste0("Comp", seq_len(comp_num))
      )
      invisible(lapply(vip_loadings, function(plot_fn) plot_fn()))

      if (!is.null(cv_opt)) {
        if (cv_opt == "loocv") {
          set.seed(seed)
          loocv_results2 <- mixOmics::perf(cytokine_splsda2, validation = "loo")
          loocv_error_rate2 <- loocv_results2$error.rate$overall[
            "comp2",
            "max.dist"
          ]
          loocv_acc2 <- 100 * signif(1 - loocv_error_rate2, digits = 2)
          if (verbose) {
            cat(paste0("LOOCV Accuracy (VIP>1): ", loocv_acc2, "%\n"))
          }

          error_rates2 <- loocv_results2$error.rate$overall[, "max.dist"]
          error_df2 <- as.data.frame(error_rates2)
          error_df2$Component <- rownames(error_df2)
          error_df2 <- reshape2::melt(
            error_df2,
            id.vars = "Component",
            variable.name = "Distance",
            value.name = "ErrorRate"
          )

          vip_CV <- ggplot2::ggplot(
            error_df2,
            ggplot2::aes(
              x = Component,
              y = ErrorRate,
              color = Distance,
              group = 1
            )
          ) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 3) +
            ggplot2::labs(
              title = paste("LOOCV Error Rate (VIP>1):", overall_analysis),
              x = "Number of Components",
              y = "Error Rate"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
            ) +
            ggplot2::scale_color_manual(values = "red", labels = "max.dist")
          print(vip_CV)
        } else if (cv_opt == "Mfold") {
          set.seed(seed)
          fold_results2 <- mixOmics::perf(
            cytokine_splsda2,
            validation = "Mfold",
            folds = fold_num,
            nrepeat = 100
          )
          fold_error_rate2 <- fold_results2$error.rate$overall[
            "comp2",
            "max.dist"
          ]
          fold_acc2 <- 100 * signif(1 - fold_error_rate2, digits = 2)
          if (verbose) {
            cat(paste0("Mfold Accuracy (VIP>1): ", fold_acc2, "%\n"))
          }

          error_rates2 <- fold_results2$error.rate$overall[, "max.dist"]
          error_df2 <- as.data.frame(error_rates2)
          error_df2$Component <- rownames(error_df2)
          error_df2 <- reshape2::melt(
            error_df2,
            id.vars = "Component",
            variable.name = "Distance",
            value.name = "ErrorRate"
          )

          vip_CV <- ggplot2::ggplot(
            error_df2,
            ggplot2::aes(
              x = Component,
              y = ErrorRate,
              color = Distance,
              group = 1
            )
          ) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 3) +
            ggplot2::labs(
              title = paste("Mfold Error Rate (VIP>1):", overall_analysis),
              x = "Number of Components",
              y = "Error Rate"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
            ) +
            ggplot2::scale_color_manual(values = "red", labels = "max.dist")
          print(vip_CV)
        }
      }

      if (conf_mat == TRUE) {
        if (verbose) {
          cat(paste0(
            "Confusion Matrix for PLS-DA Comparison: ",
            overall_analysis,
            "\n"
          ))
        }

        # Confusion Matrix for main model
        cm_overall <- caret::confusionMatrix(
          data = as.factor(prediction1[, 3]), # predicted
          reference = as.factor(prediction1[, 1]) # actual
        )
        if (verbose) {
          print(cm_overall$table)
          cat("Accuracy:", signif(cm_overall$overall["Accuracy"], 2), "\n")

          # Check if binary or multi-class
          if (nlevels(as.factor(prediction1[, 1])) == 2) {
            cat(
              "Sensitivity:",
              signif(cm_overall$byClass["Sensitivity"], 2),
              "\n"
            )
            cat(
              "Specificity:",
              signif(cm_overall$byClass["Specificity"], 2),
              "\n"
            )
          } else {
            cat("\nPer-Class Sensitivity:\n")
            print(signif(cm_overall$byClass[, "Sensitivity"], 2))
            cat("\nPer-Class Specificity:\n")
            print(signif(cm_overall$byClass[, "Specificity"], 2))
            macro_sens <- mean(
              cm_overall$byClass[, "Sensitivity"],
              na.rm = TRUE
            )
            macro_spec <- mean(
              cm_overall$byClass[, "Specificity"],
              na.rm = TRUE
            )
            cat("\nMacro-Averaged Sensitivity:", signif(macro_sens, 2), "\n")
            cat("Macro-Averaged Specificity:", signif(macro_spec, 2), "\n")
          }
        }

        if (verbose) {
          cat(paste0(
            "Confusion Matrix for PLS-DA Comparison with VIP > 1: ",
            overall_analysis,
            "\n"
          ))
        }

        # Confusion Matrix for VIP>1 model
        cm_vip <- caret::confusionMatrix(
          data = as.factor(prediction2[, 3]), # predicted
          reference = as.factor(prediction2[, 1]) # actual
        )
        if (verbose) {
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
            macro_sens_vip <- mean(
              cm_vip$byClass[, "Sensitivity"],
              na.rm = TRUE
            )
            macro_spec_vip <- mean(
              cm_vip$byClass[, "Specificity"],
              na.rm = TRUE
            )
            cat(
              "\nMacro-Averaged Sensitivity:",
              signif(macro_sens_vip, 2),
              "\n"
            )
            cat("Macro-Averaged Specificity:", signif(macro_spec_vip, 2), "\n")
          }
        }
      }

      # If roc = TRUE, compute and plot ROC curve for the overall model of VIP > 1
      if (roc) {
        vip_ROC <- mixOmics::auroc(
          object = cytokine_splsda2,
          newdata = the_data_mat,
          outcome.test = the_groups,
          plot = TRUE,
          roc.comp = comp_num,
          title = paste0("ROC Curve (VIP>1):", overall_analysis),
          print = FALSE
        )
      }
    }
    # Return a list of results
    result_list <- list(
      tune_res = tune_res,
      overall_indiv_plot = overall_indiv_plot$graph,
      overall_3D = overall_3D,
      overall_ROC = overall_ROC,
      overall_CV = overall_CV,
      loadings = loadings_list,
      vip_scores = vip_scores,
      vip_indiv_plot = if (!is.null(vip_indiv_plot)) {
        vip_indiv_plot$graph
      } else {
        NULL
      },
      vip_loadings = vip_loadings,
      vip_3D = vip_3D,
      vip_ROC = vip_ROC,
      vip_CV = vip_CV
    )
    if (!is.null(output_file)) {
      plot_list <- Filter(
        Negate(is.null),
        c(
          list(overall_indiv_plot$graph),
          if (!is.null(overall_CV)) list(overall_CV),
          unname(loadings_list), # flatten named list of closures
          unname(vip_scores), # flatten named list of closures
          list(vip_indiv_plot$graph),
          unname(vip_loadings),
          if (!is.null(overall_3D)) list(overall_3D),
          if (!is.null(vip_3D)) list(vip_3D),
          if (!is.null(overall_ROC)) {
            list(function() {
              mixOmics::auroc(
                object = cytokine_splsda,
                newdata = the_data_df,
                outcome.test = the_groups,
                plot = TRUE,
                roc.comp = comp_num,
                print = FALSE
              )
            })
          },
          if (!is.null(vip_ROC)) {
            list(function() {
              mixOmics::auroc(
                object = cytokine_splsda2,
                newdata = the_data_mat,
                outcome.test = the_groups,
                plot = TRUE,
                roc.comp = comp_num,
                print = FALSE
              )
            })
          },
          if (!is.null(vip_CV)) list(vip_CV)
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
    # Case 2: Both group and treatment columns are provided and they differ.
    levels_vec <- unique(data[[group_col2]])
    indiv_plots <- list()
    vip_indiv_plots <- list()
    for (i in seq_along(levels_vec)) {
      current_level <- levels_vec[i]
      overall_analysis <- current_level
      condt <- data[[group_col2]] == current_level

      the_data_df <- data[
        condt,
        -which(
          names(data) %in%
            c(
              group_col,
              group_col2,
              multilevel_col,
              batch_col
            )
        )
      ]
      the_data_df <- the_data_df[, sapply(the_data_df, is.numeric)]
      the_groups <- as.vector(data[[group_col]][condt])

      if (length(unique(the_groups)) < 2) {
        stop(
          "The grouping variable must have at least two levels for PLS-DA.
             Please provide an appropriate grouping column."
        )
      }
      multilevel_vec <- NULL
      if (!is.null(multilevel_col)) {
        if (!(multilevel_col %in% names(data))) {
          stop(sprintf("Multilevel column '%s' not found.", multilevel_col))
        }
        multilevel_vec <- data[[multilevel_col]]
        if (verbose) {
          message(sprintf("Using multilevel design: '%s'.\n"), multilevel_col)
        }
      }
      # Tuning sPLS-DA
      outcome <- factor(data[[group_col]])
      if (nlevels(outcome) < 2) {
        stop("Outcome must have at least two levels.")
      }
      predictors <- data[,
        setdiff(names(data), id_cols),
        drop = FALSE
      ]
      tune_res <- NULL
      if (tune) {
        tune_res <- mixOmics::tune.splsda(
          X = predictors,
          Y = outcome,
          ncomp = comp_num,
          test.keepX = c(5, 10, 15, 20, 25),
          validation = "Mfold",
          dist = "max.dist",
          nrepeat = 100,
          folds = tune_folds,
          multilevel = multilevel_vec,
          progressBar = FALSE
        )
        # Plot tuning results with proper title and labels
        p <- plot(tune_res) +
          ggplot2::ggtitle(paste("sPLS-DA Tuning Results:", overall_analysis)) +
          ggplot2::xlab("Number of Variables") +
          ggplot2::ylab("Balanced Error Rate (BER)")
        print(p)
        ncomp_final <- tune_res$choice.ncomp$ncomp
        keepX_final <- tune_res$choice.keepX[1:ncomp_final]
      }
      if (tune) {
        cytokine_splsda <- mixOmics::splsda(
          the_data_df,
          the_groups,
          scale = TRUE,
          ncomp = ncomp_final,
          keepX = keepX_final,
          multilevel = multilevel_vec
        )
      } else {
        cytokine_splsda <- mixOmics::splsda(
          the_data_df,
          the_groups,
          scale = TRUE,
          ncomp = comp_num,
          keepX = rep(var_num, comp_num),
          multilevel = multilevel_vec
        )
      }

      splsda_predict <- predict(cytokine_splsda, the_data_df, dist = "max.dist")
      prediction1 <- cbind(
        original = the_groups,
        splsda_predict$class$max.dist
      )
      accuracy1 <- (sum(prediction1[, 1] == prediction1[, 3]) /
        length(prediction1[, 1]))
      acc1 <- 100 * signif(accuracy1, digits = 2)

      # Create a grid of values
      if (bg) {
        bg_maxdist <- mixOmics::background.predict(
          cytokine_splsda,
          comp.predicted = 2,
          dist = "max.dist",
          xlim = c(-15, 15),
          ylim = c(-15, 15),
          resolution = 200
        )
      }

      group_factors <- seq_len(length(levels(factor(the_groups))))

      .build_indiv_args <- function(
        obj,
        colors,
        ind_names,
        pch_vec,
        title,
        legend_title,
        ellipse = FALSE,
        bg_obj = NULL,
        extra = list()
      ) {
        args <- c(
          list(
            obj,
            legend = TRUE,
            col = colors,
            title = title,
            legend.title = legend_title
          ),
          extra
        )

        if (ellipse) {
          args$ellipse <- TRUE
        }
        if (!is.null(bg_obj)) {
          args$background <- bg_obj
        }

        # labels OR shapes (never both)
        if (isTRUE(ind_names) || is.character(ind_names)) {
          args$ind.names <- ind_names
        } else {
          args$ind.names <- FALSE
          args$pch <- pch_vec
        }
        args
      }

      plot_args <- .build_indiv_args(
        obj = cytokine_splsda,
        colors = colors,
        ind_names = ind_names,
        pch_vec = pch_values[group_factors], # match groups present in this level
        title = paste(overall_analysis, "With Accuracy:", acc1, "%"),
        legend_title = group_col,
        ellipse = ellipse,
        bg_obj = if (bg) bg_maxdist else NULL
      )
      overall_indiv_plot <- do.call(mixOmics::plotIndiv, plot_args)
      indiv_plots[[as.character(current_level)]] <- overall_indiv_plot$graph

      overall_3D <- NULL

      if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
        cytokine_scores <- cytokine_splsda$variates$X
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
      }

      # If roc = TRUE, compute and plot ROC curve for the overall model
      overall_ROC <- NULL

      if (roc) {
        overall_ROC <- mixOmics::auroc(
          object = cytokine_splsda,
          newdata = the_data_df,
          outcome.test = the_groups,
          plot = TRUE,
          roc.comp = comp_num,
          title = paste0("ROC Curve:", overall_analysis),
          print = FALSE
        )
      }

      overall_CV <- NULL

      if (!is.null(cv_opt)) {
        if (cv_opt == "loocv") {
          set.seed(seed)
          loocv_results <- mixOmics::perf(cytokine_splsda, validation = "loo")
          loocv_error_rate <- loocv_results$error.rate$overall[
            "comp2",
            "max.dist"
          ]
          loocv_acc <- 100 * signif(1 - loocv_error_rate, digits = 2)
          if (verbose) {
            cat(paste0(current_level, " LOOCV Accuracy: ", loocv_acc, "%\n"))
          }

          error_rates <- loocv_results$error.rate$overall[, "max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(
            error_df,
            id.vars = "Component",
            variable.name = "Distance",
            value.name = "ErrorRate"
          )
          overall_CV <- ggplot2::ggplot(
            error_df,
            ggplot2::aes(
              x = Component,
              y = ErrorRate,
              color = Distance,
              group = 1
            )
          ) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 3) +
            ggplot2::labs(
              title = paste("LOOCV Error Rate:", overall_analysis),
              x = "Number of Components",
              y = "Error Rate"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
            ) +
            ggplot2::scale_color_manual(values = "red", labels = "max.dist")
          print(overall_CV)
        } else if (cv_opt == "Mfold") {
          set.seed(seed)
          fold_results <- mixOmics::perf(
            cytokine_splsda,
            validation = "Mfold",
            folds = fold_num,
            nrepeat = 100
          )
          fold_error_rate <- fold_results$error.rate$overall[
            "comp2",
            "max.dist"
          ]
          fold_acc <- 100 * signif(1 - fold_error_rate, digits = 2)
          if (verbose) {
            cat(paste0(current_level, " Mfold Accuracy: ", fold_acc, "%\n"))
          }

          error_rates <- fold_results$error.rate$overall[, "max.dist"]
          error_df <- as.data.frame(error_rates)
          error_df$Component <- rownames(error_df)
          error_df <- reshape2::melt(
            error_df,
            id.vars = "Component",
            variable.name = "Distance",
            value.name = "ErrorRate"
          )

          overall_CV <- ggplot2::ggplot(
            error_df,
            ggplot2::aes(
              x = Component,
              y = ErrorRate,
              color = Distance,
              group = 1
            )
          ) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 3) +
            ggplot2::labs(
              title = paste("Mfold Error Rate:", overall_analysis),
              x = "Number of Components",
              y = "Error Rate"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
            ) +
            ggplot2::scale_color_manual(values = "red", labels = "max.dist")
          print(overall_CV)
        }
      }

      # Loadings plot for each component

      loadings_list <- setNames(
        lapply(seq_len(comp_num), function(comp) {
          force(comp) # capture comp in the closure
          function() {
            mixOmics::plotLoadings(
              cytokine_splsda,
              comp = comp,
              contrib = "max",
              method = "mean",
              size.name = 1,
              size.legend = 1,
              legend.color = colors,
              title = paste(
                "Loadings for Component",
                comp,
                ":",
                overall_analysis
              ),
              size.title = 1,
              legend.title = group_col
            )
          }
        }),
        nm = paste0("Comp", seq_len(comp_num))
      )
      invisible(lapply(loadings_list, function(plot_fn) plot_fn()))

      all_vip_scores <- mixOmics::vip(cytokine_splsda)
      vip_scores <- setNames(
        lapply(seq_len(comp_num), function(comp) {
          force(comp) # capture `comp` in the closure
          function() {
            # recreate data frame
            vscore <- as.data.frame(all_vip_scores[, comp, drop = FALSE])
            vscore$metabo <- rownames(vscore)
            vscore$comp <- vscore[, 1]
            bar <- vscore[
              order(vscore$comp, decreasing = TRUE),
              c("metabo", "comp")
            ]
            bar <- bar[is.finite(bar$comp), , drop = FALSE]
            bar <- bar[bar$comp > 0, , drop = FALSE]

            # Match number with keepX
            if (tune) {
              bar <- head(bar, min(keepX_final, nrow(bar)))
            } else {
              bar <- head(bar, min(var_num, nrow(bar)))
            }
            # Lock the ordering to exactly what's being plotted (prevents gaps)
            bar$metabo <- factor(bar$metabo, levels = rev(bar$metabo))

            # build and print the plot
            p <- ggplot2::ggplot(bar, ggplot2::aes(x = metabo, y = comp)) +
              ggplot2::geom_bar(stat = "identity", position = "dodge") +
              ggplot2::scale_y_continuous(limits = c(0, max(bar$comp))) +
              ggplot2::geom_hline(yintercept = 1, color = "grey") +
              ggplot2::scale_x_discrete(limits = factor(bar$metabo)) +
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                  angle = 45,
                  hjust = 1,
                  size = 15
                ),
                panel.grid = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(
                  color = "black",
                  fill = "transparent"
                )
              ) +
              ggplot2::coord_flip() +
              ggplot2::labs(x = "Variables", y = "VIP score") +
              ggplot2::ggtitle(paste("Component", comp))

            print(p)
          }
        }),
        nm = paste0("Comp", seq_len(comp_num))
      )
      invisible(lapply(vip_scores, function(draw_fn) draw_fn()))

      vip_indiv_plot <- NULL
      vip_loadings <- NULL
      vip_3D <- NULL
      vip_ROC <- NULL
      vip_CV <- NULL

      condt_variable <- all_vip_scores[, 1] > 1
      keep_x <- sum(condt_variable)
      the_data_mat <- the_data_df[, condt_variable, drop = FALSE]
      cytokine_splsda2 <- mixOmics::splsda(
        the_data_mat,
        the_groups,
        scale = TRUE,
        ncomp = comp_num,
        keepX = rep(keep_x, comp_num),
        multilevel = multilevel_vec
      )

      splsda_predict2 <- predict(
        cytokine_splsda2,
        the_data_mat,
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
      if (bg) {
        bg_maxdist2 <- mixOmics::background.predict(
          cytokine_splsda2,
          comp.predicted = 2,
          dist = "max.dist",
          xlim = c(-15, 15),
          ylim = c(-15, 15),
          resolution = 200
        )
      }

      plot_args2 <- .build_indiv_args(
        obj = cytokine_splsda2,
        colors = colors,
        ind_names = ind_names,
        pch_vec = pch_values[group_factors], # keep consistent with main loop
        title = paste(overall_analysis, "(VIP>1)", "With Accuracy:", acc2, "%"),
        legend_title = group_col,
        ellipse = ellipse,
        bg_obj = if (bg) bg_maxdist2 else NULL
      )
      vip_indiv_plot <- do.call(mixOmics::plotIndiv, plot_args2)
      vip_indiv_plots[[as.character(current_level)]] <- vip_indiv_plot$graph

      if (!is.null(style) && comp_num == 3 && (tolower(style) == "3d")) {
        cytokine_scores2 <- cytokine_splsda2$variates$X
        vip_3D <- function() {
          plot3D::scatter3D(
            cytokine_scores2[, 1],
            cytokine_scores2[, 2],
            cytokine_scores2[, 3],
            pch = pch_values,
            col = colors,
            xlab = "Component 1",
            ylab = "Component 2",
            zlab = "Component 3",
            main = paste("3D Plot:", overall_analysis, "(VIP>1)"),
            theta = 20,
            phi = 30,
            bty = "g",
            colkey = FALSE
          )
        }
        vip_3D()
      }
      # Loadings plot for each component with VIP > 1

      vip_loadings <- setNames(
        lapply(seq_len(comp_num), function(comp) {
          force(comp) # capture comp in the closure
          function() {
            mixOmics::plotLoadings(
              cytokine_splsda2,
              comp = comp,
              contrib = "max",
              method = "mean",
              size.name = 1,
              size.legend = 1,
              legend.color = colors,
              title = paste(
                "Loadings for Component",
                comp,
                "(VIP > 1):",
                overall_analysis
              ),
              size.title = 1,
              legend.title = group_col
            )
          }
        }),
        nm = paste0("Comp", seq_len(comp_num))
      )
      invisible(lapply(vip_loadings, function(plot_fn) plot_fn()))

      if (roc) {
        vip_ROC <- mixOmics::auroc(
          object = cytokine_splsda2,
          newdata = the_data_mat,
          outcome.test = the_groups,
          plot = TRUE,
          roc.comp = comp_num,
          title = paste0("ROC Curve (VIP>1):", overall_analysis),
          print = FALSE
        )
      }

      if (!is.null(cv_opt)) {
        if (cv_opt == "loocv") {
          set.seed(seed)
          loocv_results2 <- mixOmics::perf(cytokine_splsda2, validation = "loo")
          loocv_error_rate2 <- loocv_results2$error.rate$overall[
            "comp2",
            "max.dist"
          ]
          loocv_acc2 <- 100 * signif(1 - loocv_error_rate2, digits = 2)
          if (verbose) {
            cat(paste0(
              current_level,
              " LOOCV Accuracy (VIP>1): ",
              loocv_acc2,
              "%\n"
            ))
          }

          error_rates2 <- loocv_results2$error.rate$overall[, "max.dist"]
          error_df2 <- as.data.frame(error_rates2)
          error_df2$Component <- rownames(error_df2)
          error_df2 <- reshape2::melt(
            error_df2,
            id.vars = "Component",
            variable.name = "Distance",
            value.name = "ErrorRate"
          )
          vip_CV <- ggplot2::ggplot(
            error_df2,
            ggplot2::aes(
              x = Component,
              y = ErrorRate,
              color = Distance,
              group = 1
            )
          ) +
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
          print(vip_CV)
        } else if (cv_opt == "Mfold") {
          set.seed(seed)
          fold_results2 <- mixOmics::perf(
            cytokine_splsda2,
            validation = "Mfold",
            folds = fold_num,
            nrepeat = 100
          )
          fold_error_rate2 <- fold_results2$error.rate$overall[
            "comp2",
            "max.dist"
          ]
          fold_acc2 <- 100 * signif(1 - fold_error_rate2, digits = 2)
          if (verbose) {
            cat(paste0(
              current_level,
              " Mfold Accuracy (VIP>1): ",
              fold_acc2,
              "%\n"
            ))
          }

          error_rates2 <- fold_results2$error.rate$overall[, "max.dist"]
          error_df2 <- as.data.frame(error_rates2)
          error_df2$Component <- rownames(error_df2)
          error_df2 <- reshape2::melt(
            error_df2,
            id.vars = "Component",
            variable.name = "Distance",
            value.name = "ErrorRate"
          )
          vip_CV <- ggplot2::ggplot(
            error_df2,
            ggplot2::aes(
              x = Component,
              y = ErrorRate,
              color = Distance,
              group = 1
            )
          ) +
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
          print(vip_CV)
        }
      }
      if (conf_mat == TRUE) {
        if (verbose) {
          cat(paste0(
            "Confusion Matrix for PLS-DA Comparison: ",
            current_level,
            "\n"
          ))
        }

        # Confusion Matrix for main model
        cm_overall <- caret::confusionMatrix(
          data = as.factor(prediction1[, 3]), # predicted
          reference = as.factor(prediction1[, 1]) # actual
        )
        if (verbose) {
          print(cm_overall$table)
          cat("Accuracy:", signif(cm_overall$overall["Accuracy"], 2), "\n")
          if (nlevels(as.factor(prediction1[, 1])) == 2) {
            cat(
              "Sensitivity:",
              signif(cm_overall$byClass["Sensitivity"], 2),
              "\n"
            )
            cat(
              "Specificity:",
              signif(cm_overall$byClass["Specificity"], 2),
              "\n"
            )
          } else {
            cat("\nPer-Class Sensitivity:\n")
            print(signif(cm_overall$byClass[, "Sensitivity"], 2))
            cat("\nPer-Class Specificity:\n")
            print(signif(cm_overall$byClass[, "Specificity"], 2))
            macro_sens <- mean(
              cm_overall$byClass[, "Sensitivity"],
              na.rm = TRUE
            )
            macro_spec <- mean(
              cm_overall$byClass[, "Specificity"],
              na.rm = TRUE
            )
            cat("\nMacro-Averaged Sensitivity:", signif(macro_sens, 2), "\n")
            cat("Macro-Averaged Specificity:", signif(macro_spec, 2), "\n")
          }
        }

        if (verbose) {
          cat(paste0(
            "Confusion Matrix for PLS-DA Comparison with VIP > 1: ",
            current_level,
            "\n"
          ))
        }

        # Confusion Matrix for VIP>1 model
        cm_vip <- caret::confusionMatrix(
          data = as.factor(prediction2[, 3]), # predicted
          reference = as.factor(prediction2[, 1]) # actual
        )
        if (verbose) {
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
            macro_sens_vip <- mean(
              cm_vip$byClass[, "Sensitivity"],
              na.rm = TRUE
            )
            macro_spec_vip <- mean(
              cm_vip$byClass[, "Specificity"],
              na.rm = TRUE
            )
            cat(
              "\nMacro-Averaged Sensitivity:",
              signif(macro_sens_vip, 2),
              "\n"
            )
            cat("Macro-Averaged Specificity:", signif(macro_spec_vip, 2), "\n")
          }
        }
      }
    }
    # Return a list of results
    result_list <- list(
      tune_res = tune_res,
      overall_indiv_plot = overall_indiv_plot$graph,
      all_indiv_plots = indiv_plots,
      overall_3D = overall_3D,
      overall_ROC = overall_ROC,
      overall_CV = overall_CV,
      loadings = loadings_list,
      vip_scores = vip_scores,
      vip_indiv_plot = if (!is.null(vip_indiv_plot)) {
        vip_indiv_plot$graph
      } else {
        NULL
      },
      vip_indiv_plots = vip_indiv_plots,
      vip_loadings = vip_loadings,
      vip_3D = vip_3D,
      vip_ROC = vip_ROC,
      vip_CV = vip_CV
    )
    if (!is.null(output_file)) {
      plot_list <- Filter(
        Negate(is.null),
        c(
          list(overall_indiv_plot$graph),
          if (!is.null(overall_CV)) list(overall_CV),
          unname(loadings_list), # flatten named list of closures
          unname(vip_scores), # flatten named list of closures
          list(vip_indiv_plot$graph),
          unname(vip_loadings),
          if (!is.null(overall_3D)) list(overall_3D),
          if (!is.null(vip_3D)) list(vip_3D),
          if (!is.null(overall_ROC)) {
            list(function() {
              mixOmics::auroc(
                object = cytokine_splsda,
                newdata = the_data_df,
                outcome.test = the_groups,
                plot = TRUE,
                roc.comp = comp_num,
                print = FALSE
              )
            })
          },
          if (!is.null(vip_ROC)) {
            list(function() {
              mixOmics::auroc(
                object = cytokine_splsda2,
                newdata = the_data_mat,
                outcome.test = the_groups,
                plot = TRUE,
                roc.comp = comp_num,
                print = FALSE
              )
            })
          },
          if (!is.null(vip_CV)) list(vip_CV)
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
  }
}
