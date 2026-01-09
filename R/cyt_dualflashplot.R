#' Dual-flashlight Plot
#'
#' @description
#' This function reshapes the input data and computes summary statistics
#' (mean and variance) for each variable grouped by a specified factor column.
#' It then calculates the SSMD (Strictly Standardized Mean Difference)
#' and log2 fold change between two groups (group1 and group2) and categorizes
#' the effect strength as "Strong Effect", "Moderate Effect", or "Weak Effect".
#' A dual flash plot is generated using ggplot2 where the x-axis represents
#' the average log2 fold change and the y-axis represents the SSMD.
#' Additionally, the function prints the computed statistics to the console.
#' A ggplot object is returned for further modification.
#' @param data A data frame containing at least one numeric column
#'   and a grouping column.
#' @param group_var Character.  Name of the grouping column.
#' @param group1 Character strings identifying the two levels
#'   of `group_var` to compare.
#' @param group2 Character strings identifying the two levels
#'   of `group_var` to compare.
#' @param ssmd_thresh Numeric.  Absolute SSMD threshold for
#'   highlighting observations as significant.  Default is 1.
#' @param log2fc_thresh Numeric.  Absolute log2 fold change threshold
#'   for significance.  Default is 1.
#' @param top_labels Integer.  Number of variables with the largest
#'   absolute SSMD to label on the plot.  Default is 15.
#' @param category_labels Optional named character vector of length
#'   three providing alternative labels for the SSMD effect
#'   categories.  Names must include "Strong", "Moderate" and
#'   "Weak".  Defaults are `c(Strong = "Strong Effect", Moderate =
#'   "Moderate Effect", Weak = "Weak Effect")`.
#' @param colors Optional named character vector mapping the effect
#'   categories to colors.  Names must match those in
#'   `category_labels`.  Defaults are red for strong, orange for
#'   moderate and blue for weak.
#' @param shapes Optional named numeric vector of length two giving
#'   the plotting characters for non‑significant and significant
#'   points.  Defaults are `c(`FALSE` = 16, `TRUE` = 17)` (solid and
#'   triangular symbols).
#' @param verbose Logical.  If `TRUE`, prints the computed summary
#'   statistics.  Default is `FALSE`.
#'
#' @return A ggplot object representing the dual‑flashlight plot.
#'   When `verbose` is `TRUE`, the summary statistics data frame is
#'   printed to the console.
#'
#' @author Xiaohua Douglas Zhang and Shubh Saraswat
#'
#' @export
#' @import dplyr
#' @importFrom tidyr pivot_longer pivot_wider
#' @import ggplot2
#'
#' @examples
#' # Loading data
#' data_df <- ExampleData1[, -c(2:3)]
#' cyt_dualflashplot(data_df, group_var = "Group", group1 = "T2D", group2 = "ND",
#'                        ssmd_thresh = 0.2, log2fc_thresh = 1,
#'                        top_labels = 10)
#'
cyt_dualflashplot <- function(
  data,
  group_var,
  group1,
  group2,
  ssmd_thresh = 1,
  log2fc_thresh = 1,
  top_labels = 15,
  category_labels = c(
    Strong = "Strong Effect",
    Moderate = "Moderate Effect",
    Weak = "Weak Effect"
  ),
  colors = c(
    "Strong Effect" = "red",
    "Moderate Effect" = "orange",
    "Weak Effect" = "blue"
  ),
  shapes = c(`FALSE` = 16, `TRUE` = 17),
  verbose = FALSE
) {
  names(data) <- make.names(names(data), unique = TRUE)
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }
  if (!group_var %in% names(data)) {
    stop(paste0("group_var '", group_var, "' not found in data."))
  }
  # Check groups
  if (!all(c(group1, group2) %in% unique(data[[group_var]]))) {
    stop("Both group1 and group2 must be levels of group_var.")
  }
  # Select numeric variables
  num_vars <- names(data)[sapply(data, is.numeric) & names(data) != group_var]
  if (length(num_vars) == 0) {
    stop("No numeric variables found.")
  }
  # Pivot longer and compute means and variances
  data_long <- data |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of(group_var),
      names_to = "cytokine",
      values_to = "level"
    )
  stats_df <- data_long |>
    dplyr::group_by(cytokine, .data[[group_var]]) |>
    dplyr::summarise(
      mean = mean(level, na.rm = TRUE),
      variance = var(level, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = .data[[group_var]],
      values_from = c(mean, variance)
    )
  # Compute SSMD and log2FC
  stats_df <- stats_df |>
    dplyr::mutate(
      ssmd = (get(paste0("mean_", group1)) - get(paste0("mean_", group2))) /
        sqrt(
          (get(paste0("variance_", group1)) +
            get(paste0("variance_", group2))) /
            2
        ),
      log2FC = log2(
        get(paste0("mean_", group1)) / get(paste0("mean_", group2))
      ),
      Effect = dplyr::case_when(
        abs(ssmd) >= 1 ~ category_labels["Strong"],
        abs(ssmd) >= 0.5 ~ category_labels["Moderate"],
        TRUE ~ category_labels["Weak"]
      ),
      Significant = (abs(ssmd) >= ssmd_thresh) & (abs(log2FC) >= log2fc_thresh)
    )
  # Ensure effect labels factor levels match provided colors
  stats_df$Effect <- factor(stats_df$Effect, levels = category_labels)
  # Determine colors (extend if necessary)
  if (!is.null(colors)) {
    # If names missing, use defaults
    for (lab in category_labels) {
      if (!lab %in% names(colors)) {
        colors[lab] <- colors[[1]]
      }
    }
  }
  # Determine shapes vector
  shapes <- shapes[c("FALSE", "TRUE")]
  # Build the plot
  p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = log2FC, y = ssmd)) +
    ggplot2::geom_point(
      ggplot2::aes(color = Effect, shape = as.factor(Significant)),
      size = 2
    ) +
    ggplot2::geom_vline(
      xintercept = c(log2fc_thresh, -log2fc_thresh),
      linetype = "dashed",
      color = "blue"
    ) +
    ggplot2::geom_hline(
      yintercept = c(ssmd_thresh, -ssmd_thresh),
      linetype = "dashed",
      color = "blue"
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::labs(
      x = "Average log2 Fold Change",
      y = "SSMD",
      title = paste0("SSMD vs log2FC for ", group1, " vs ", group2),
      shape = "Significant"
    ) +
    ggplot2::theme_minimal()
  # Label top variables
  if (top_labels > 0) {
    top_df <- stats_df |>
      dplyr::arrange(dplyr::desc(abs(ssmd))) |>
      dplyr::slice_head(n = top_labels)
    p <- p +
      ggrepel::geom_text_repel(
        data = top_df,
        ggplot2::aes(label = cytokine),
        size = 3,
        vjust = 1.5,
        hjust = 1.1
      )
  }
  # Print statistics if verbose
  if (verbose) {
    print(as.data.frame(stats_df))
  }
  return(p)
}
