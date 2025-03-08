#' Heat Map.
#'
#' @param data A data frame containing the input data. Only numeric columns
#' will be used to generate the heatmap.
#' @param scale Character. An optional scaling option. If set to "log2",
#' the numeric data will be log2-transformed (with non-positive values
#' set to NA). Default is NULL.
#' @param annotation_col_name Character. An optional column name from
#' \code{data} to be used for generating annotation colors. Default is NULL.
#' @param title Character. The title of the heatmap and the file name for
#' saving the plot. The file extension (".pdf" or ".png") determines the
#' output format.
#'
#' @description
#' This function creates a heatmap using the numeric columns from the
#' provided data frame. If requested via the \code{scale} parameter,
#' the function applies a log2 transformation to the data (with non-positive
#' values replaced by NA). The heatmap is saved as a file,
#' with the format determined by the file extension in \code{title}.
#'
#' @return The function does not return a value. It saves the heatmap to a file.
#'
#' @export
#' @import gplots
#'
#' @examples
#' # Load sample data
#' data("cytodata")
#' data_df <- cytodata
#' # Generate a heatmap with log2 scaling and annotation based on
#' # the "Group" column
#' cyt_heatmap(
#'   data = data_df[, -c(1,3,4)],
#'   scale = "log2",  # Optional scaling
#'   annotation_col_name = "Group",
#'   title = "Heatmap.png"
#' )
#'
cyt_heatmap <- function(data, scale = NULL, annotation_col_name = NULL, title) {
  if (grepl("\\.pdf$", title, ignore.case = TRUE)) {
    pdf(file = title)
  } else if (grepl("\\.png$", title, ignore.case = TRUE)) {
    png(filename = title, res = 300, width = 2100, height = 2100, units = "px")
  } else {
    stop("Title must end with .pdf or .png")
  }

  # Ensure data is a data frame and extract only numeric data
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }
  numeric_data <- data[, sapply(data, is.numeric), drop = FALSE]
  if (ncol(data) != ncol(numeric_data)) {
    warning("Non-numeric columns detected. Subsetting to
            numeric columns only.")
  }

  # Apply log2 transformation if requested
  if (!is.null(scale) && scale == "log2") {
    numeric_data[numeric_data <= 0] <- NA  # Set non-positive values to NA
    numeric_data <- log2(numeric_data)
  }

  # Generate annotation colors if annotation_col_name is provided
  # and exists in data
  side_colors <- NULL
  if (!is.null(annotation_col_name) && annotation_col_name %in% names(data)) {
    ann_data <- as.factor(data[[annotation_col_name]])
    # Check if we have enough annotations for each numeric column
    if (length(ann_data) >= ncol(numeric_data)) {
      num_levels <- length(levels(ann_data))
      # Generate a color for each level and subset to match the number of
      # numeric cols
      side_colors <- rainbow(num_levels)[as.integer(
        ann_data[seq_len(ncol(numeric_data))])]
    } else {
      warning("Length of annotation column is less than the number of
               numeric columns. ",
              "Skipping annotation colors.")
    }
  }

  # Build argument list for heatmap.2: only include ColSideColors if
  # side_colors is not NULL
  heatmap_args <- list(
    as.matrix(numeric_data),
    distfun = function(x) dist(x, method = "euclidean"),
    hclustfun = function(x) hclust(x, method = "complete"),
    dendrogram = "both",
    trace = "column",
    key = TRUE,
    cexCol = 1,
    margins = c(10, 10)
  )
  if (!is.null(side_colors)) {
    heatmap_args$ColSideColors <- side_colors
  }

  do.call(gplots::heatmap.2, heatmap_args)
  if (dev.cur() > 1) dev.off()
}
