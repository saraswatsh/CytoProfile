#' Heat Map
#'
#' @param data A data frame to be entered.
#' @param scale Log2 transformation for scaling. Default set to NULL.
#' @param annotation_col_name A string of a column name to be chosen for annotation. Default set to NULL.
#' @param title Title of the heatmap to be saved. Can save it as png or pdf.
#' @description
#' This function takes a data frame and subsets the numeric columns or variables to generate a heatmap.
#' @return Prints a heatmap.
#' @export
#'
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
