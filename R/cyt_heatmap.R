#' Heat Map
#'
#' @param data A data frame.  Only numeric columns are used to
#'   construct the heat map.
#' @param scale Character specifying an optional scaling.  Accepts
#'   `NULL` (no scaling), "log2", "log10", "row_zscore",
#'   "col_zscore" or "zscore" (apply both row and column
#'   z‑scoring).  Default is `NULL`.
#' @param annotation_col Optional.  Either the name of a column in
#'   `data` or a vector of length equal to the number of rows or
#'   columns of the numeric matrix.  If a column name is supplied
#'   the function determines whether it annotates rows or columns based
#'   on its length or the value of `annotation_side`.
#' @param annotation_side Character.  One of "auto", "row" or
#'   "col".  When "auto" (default) the side is determined by
#'   matching the length of `annotation_col` to rows or columns.
#' @param show_row_names Logical.  If `TRUE` row names are shown
#'   Default is `FALSE`.
#' @param show_col_names Logical.  If `FALSE` column names are
#'   hidden.  Default is `TRUE`.
#' @param fontsize_row, fontsize_col Numeric.  Font sizes for row and
#'   column names respectively.  Default is 10.
#' @param cluster_rows Logical.  If `TRUE` (default), rows are
#'   clustered.
#' @param cluster_cols Logical.  If `TRUE` (default), columns are
#'   clustered.
#' @param title Character.  The heat map title or file name.  If
#'   `title` ends with ".pdf" or ".png" (case insensitive), the
#'   heat map is saved to that file and no title is printed on
#'   screen.  If `NULL` (default), the heat map is drawn on the
#'   active device without saving and without a main title.
#'
#' @description
#' This function creates a heatmap using the numeric columns from the
#' provided data frame. It provides the ability to hide row and
#' column names, adjust font sizes and clustering, and apply
#' additional transformations such as log₁₀ or combined z‑scoring.  A
#' file name with extension may be provided via `title` to save the
#' heat map to disk; otherwise the plot is drawn on the active
#' graphics device.
#'
#' @return Invisibly returns the pheatmap object created by
#'   `pheatmap::pheatmap()`.
#' @author Shubh Saraswat
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
#'
#' @examples
#' # Load sample data
#' data("ExampleData1")
#' data_df <- ExampleData1
#' # Generate a heatmap with log2 scaling and annotation based on
#' # the "Group" column
#' cyt_heatmap(
#'   data = data_df[, -c(2:3)],
#'   scale = "log2",  # Optional scaling
#'   annotation_col = "Group",
#'   annotation_side = "auto",
#'   title = NULL
#' )
#'
cyt_heatmap <- function(
  data,
  scale = c(NULL, "log2", "log10", "row_zscore", "col_zscore", "zscore"),
  annotation_col = NULL,
  annotation_side = c("auto", "row", "col"),
  show_row_names = FALSE,
  show_col_names = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  title = NULL
) {
  # Match arguments
  scale <- match.arg(scale)
  annotation_side <- match.arg(annotation_side)
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }
  # Extract numeric columns
  num <- data[, vapply(data, is.numeric, logical(1)), drop = FALSE]
  if (!ncol(num)) {
    stop("No numeric columns found in `data`.")
  }
  mat <- as.matrix(num)
  if (is.null(rownames(mat))) {
    rownames(mat) <- seq_len(nrow(mat))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  }
  # Apply scaling
  pm_scale <- "none"
  if (!is.null(scale)) {
    if (scale == "log2") {
      mat <- log2(mat)
    } else if (scale == "log10") {
      mat <- log10(mat)
    } else if (scale == "row_zscore") {
      pm_scale <- "row"
    } else if (scale == "col_zscore") {
      pm_scale <- "column"
    } else if (scale == "zscore") {
      # Standardise rows then columns
      mat <- t(apply(mat, 1, function(x) {
        mu <- mean(x, na.rm = TRUE)
        sdv <- sd(x, na.rm = TRUE)
        if (is.na(sdv) || sdv == 0) {
          return(x - mu)
        }
        (x - mu) / sdv
      }))
      mat <- apply(mat, 2, function(x) {
        mu <- mean(x, na.rm = TRUE)
        sdv <- sd(x, na.rm = TRUE)
        if (is.na(sdv) || sdv == 0) {
          return(x - mu)
        }
        (x - mu) / sdv
      })
    }
  }
  # Annotation handling
  ann_row <- ann_col <- NULL
  ann_colors <- NULL
  ann_title <- NULL
  if (!is.null(annotation_col)) {
    if (
      is.character(annotation_col) &&
        length(annotation_col) == 1 &&
        annotation_col %in% names(data)
    ) {
      ann <- factor(data[[annotation_col]])
      ann_title <- annotation_col
    } else if (length(annotation_col) %in% c(nrow(mat), ncol(mat))) {
      ann <- factor(annotation_col)
      ann_title <- "Annotation"
    } else {
      warning(
        "`annotation_col` must be a column in `data` or a vector matching rows or columns; skipping annotation."
      )
      ann <- NULL
    }
    if (!is.null(ann)) {
      # Determine side if auto
      side <- if (annotation_side == "auto") {
        if (length(ann) == nrow(mat)) {
          "row"
        } else if (length(ann) == ncol(mat)) {
          "col"
        } else {
          "row"
        }
      } else {
        annotation_side
      }
      levs <- levels(ann)
      cols <- grDevices::rainbow(length(levs))
      cmap <- stats::setNames(cols, levs)
      if (side == "row" && length(ann) == nrow(mat)) {
        ann_row <- stats::setNames(
          data.frame(ann, row.names = rownames(mat)),
          ann_title
        )
        ann_colors <- list()
        ann_colors[[ann_title]] <- cmap
      } else if (side == "col" && length(ann) == ncol(mat)) {
        ann_col <- stats::setNames(
          data.frame(ann, row.names = colnames(mat)),
          ann_title
        )
        ann_colors <- list()
        ann_colors[[ann_title]] <- cmap
      } else {
        warning(
          "`annotation_col` length does not match the chosen side; skipping annotation."
        )
      }
    }
  }
  # Determine filename and main title
  filename <- if (
    !is.null(title) && grepl("\\.(pdf|png)$", title, ignore.case = TRUE)
  ) {
    title
  } else {
    NA
  }
  main <- if (
    !is.null(title) && !grepl("\\.(pdf|png)$", title, ignore.case = TRUE)
  ) {
    title
  } else {
    NA
  }
  # Draw heat map
  ph <- pheatmap::pheatmap(
    mat,
    scale = pm_scale,
    color = grDevices::colorRampPalette(c("#253494", "#f7f7f7", "#b30000"))(
      255
    ),
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    border_color = NA,
    annotation_row = ann_row,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    legend = TRUE,
    annotation_legend = TRUE,
    show_rownames = show_row_names,
    show_colnames = show_col_names,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    filename = filename,
    main = main
  )
  invisible(ph)
}
