#' Heat Map
#'
#' @param data A data frame containing the input data. Only numeric columns will be used to generate the heatmap.
#' @param scale Character. An optional scaling option. If set to "log2", the numeric data will be log2-transformed (with non-positive values set to NA). Default is NULL.
#' @param annotation_col_name Character. An optional column name from \code{data} to be used for generating annotation colors. Default is NULL.
#' @param title Character. The title of the heatmap and the file name for saving the plot. The file extension (".pdf" or ".png") determines the output format.
#'
#' @description
#' This function creates a heatmap using the numeric columns from the provided data frame. If requested via
#' the \code{scale} parameter, the function applies a log2 transformation to the data (with non-positive values replaced by NA).
#' Optionally, if an annotation column is specified and exists in \code{data}, the function attempts to generate a color annotation
#' (although the annotation is not passed to \code{heatmap.2} in the current implementation). The heatmap is saved as a file,
#' with the format determined by the file extension in \code{title}.
#'
#' @return The function does not return a value. It saves the heatmap to a file.
#'
#' @export
#' @import gplots
#' @examples
#' cyt.heatmap(cytodata[,-4], scale = "log2", annotation_col_name = "Group")

cyt.heatmap <- function(data, scale = NULL, annotation_col_name = NULL, title) {

  if(grepl(".pdf", title)){
    pdf(file = title)
  }
  else if(grepl(".png", title)){
    png(filename = title, res = 300, width = 2100, height = 2100, units = "px")
  }
  # Ensure data is a data frame and extract only numeric data
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }

  numeric_data <- data[, sapply(data, is.numeric)]
  if (ncol(data) != ncol(numeric_data)) {
    warning("Non-numeric columns detected. Subsetting to numeric columns only.")
  }

  # Apply log2 transformation if requested
  if (!is.null(scale) && scale == "log2") {
    numeric_data[numeric_data <= 0] <- NA  # Set non-positive values to NA
    numeric_data <- log2(numeric_data)
  }
  # Generate the annotation color bar v2
  side_colors <- NULL
  palette_colors <- NULL  # Initialize a variable for palette colors
  if (!is.null(annotation_col_name) && annotation_col_name %in% names(data)) {
    ann_data <- as.factor(data[[annotation_col_name]])
    num_levels <- length(levels(ann_data))
    palette_colors <- rainbow(num_levels)  # Generate colors
    side_colors <- rainbow(num_levels)  # Use rainbow to generate a color for each level
    # Ensure the correct assignment of colors to the side bar
    side_colors <- side_colors[as.integer(ann_data[1:ncol(numeric_data)])]
  }


  # Define the color palette for the heatmap
  max_colors <- min(9, length(unique(side_colors)))

  # Plot the heatmap using heatmap.2
  gplots::heatmap.2(as.matrix(numeric_data),
                    distfun = function(x) dist(x, method = "euclidean"),
                    hclustfun = function(x) hclust(x, method = "complete"),
                    dendrogram = "both",
                    #ColSideColors = as.vector(side_colors),
                    #col = heatmap_colors,
                    trace = "column",
                    key = TRUE,
                    cexCol = 1,
                    margins = c(10,10))
  dev.off()
}
