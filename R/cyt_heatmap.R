#' Heat Map.
#'
#' @param data A data frame containing the input data. Only numeric columns
#' will be used to generate the heatmap.
#' @param scale Character. An optional scaling option. Options are NULL
#' (no scaling), "log2" (log2 transformation), "row_zscore" (z-score scaling by row),
#' or "col_zscore" (z-score scaling by column). Default is NULL.
#' @param annotation_col Character. An optional column name from
#' \code{data} to be used for generating annotation colors. Default is NULL.
#' @param annotation_side Character. Specifies whether the annotation should
#' be applied to rows or columns. Options are "auto", "row", or "col".
#' @param title Character. The title of the heatmap and the file name for
#' saving the plot. The file extension (".pdf" or ".png") determines the
#' output format. If \code{NULL}, the plot is generated on the current
#' graphics device. Default is \code{NULL}.
#'
#' @description
#' This function creates a heatmap using the numeric columns from the
#' provided data frame. It supports various scaling options and allows for row or
#' column annotations. The heatmap is saved as a file,
#' with the format determined by the file extension in \code{title}.
#'
#' @return The function does not return a value. It saves the heatmap to a file.
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
  scale = c(NULL, "log2", "row_zscore", "col_zscore"),
  annotation_col = NULL,
  annotation_side = c("auto", "row", "col"),
  title = NULL
) {
  scale <- match.arg(scale)
  annotation_side <- match.arg(annotation_side)

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }
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

  # ---- transform / scale (mapped to pheatmap's 'scale') ----
  pm_scale <- "none"
  if (identical(scale, "log2")) {
    mat <- log2(mat + 0.005) # small offset to avoid log(0)
  } else if (identical(scale, "row_zscore")) {
    pm_scale <- "row"
  } else if (identical(scale, "col_zscore")) {
    pm_scale <- "column"
  }

  # ---- annotation handling ----
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
      ann <- NULL
      warning(
        "`annotation_col` must be a column in `data` or a vector matching rows or columns; skipping."
      )
    }

    if (!is.null(ann)) {
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

  # ---- draw with pheatmap ----
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

  pheatmap::pheatmap(
    mat,
    scale = pm_scale,
    color = grDevices::colorRampPalette(c("#253494", "#f7f7f7", "#b30000"))(
      255
    ),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    annotation_row = ann_row,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    legend = TRUE,
    annotation_legend = TRUE,
    filename = filename,
    main = main
  )

  invisible(list(annotation_map = ann_colors))
}
