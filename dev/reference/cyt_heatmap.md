# Heat Map.

This function creates a heatmap using the numeric columns from the
provided data frame. It supports various scaling options and allows for
row or column annotations. The heatmap is saved as a file, with the
format determined by the file extension in `title`.

## Usage

``` r
cyt_heatmap(
  data,
  scale = c(NULL, "log2", "row_zscore", "col_zscore"),
  annotation_col = NULL,
  annotation_side = c("auto", "row", "col"),
  title = NULL
)
```

## Arguments

- data:

  A data frame containing the input data. Only numeric columns will be
  used to generate the heatmap.

- scale:

  Character. An optional scaling option. Options are NULL (no scaling),
  "log2" (log2 transformation), "row_zscore" (z-score scaling by row),
  or "col_zscore" (z-score scaling by column). Default is NULL.

- annotation_col:

  Character. An optional column name from `data` to be used for
  generating annotation colors. Default is NULL.

- annotation_side:

  Character. Specifies whether the annotation should be applied to rows
  or columns. Options are "auto", "row", or "col".

- title:

  Character. The title of the heatmap and the file name for saving the
  plot. The file extension (".pdf" or ".png") determines the output
  format. If `NULL`, the plot is generated on the current graphics
  device. Default is `NULL`.

## Value

The function does not return a value. It saves the heatmap to a file.

## Author

Shubh Saraswat

## Examples

``` r
# Load sample data
data("ExampleData1")
data_df <- ExampleData1
# Generate a heatmap with log2 scaling and annotation based on
# the "Group" column
cyt_heatmap(
  data = data_df[, -c(2:3)],
  scale = "log2",  # Optional scaling
  annotation_col = "Group",
  annotation_side = "auto",
  title = NULL
)

```
