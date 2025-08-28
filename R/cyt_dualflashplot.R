#' Dual-flashlight Plot.
#'
#' This function reshapes the input data and computes summary statistics
#' (mean and variance) for each variable grouped by a specified factor column.
#' It then calculates the SSMD (Strictly Standardized Mean Difference)
#' and log2 fold change between two groups (group1 and group2) and categorizes
#' the effect strength as "Strong Effect", "Moderate Effect", or "Weak Effect".
#' A dual flash plot is generated using ggplot2 where the x-axis represents
#' the average log2 fold change and the y-axis represents the SSMD.
#' Additionally, the function prints the computed statistics to the console.
#'
#' @param data A data frame containing the input data.
#' @param group_var A string specifying the name of the grouping
#' column in the data.
#' @param group1 A string representing the name of the first group
#' for comparison.
#' @param group2 A string representing the name of the second group
#' for comparison.
#' @param ssmd_thresh A numeric threshold for the SSMD value used to
#' determine significance.
#'   Default is 1.
#' @param log2fc_thresh A numeric threshold for the log2 fold change
#' used to determine significance.
#'   Default is 1.
#' @param top_labels An integer specifying the number of top variables
#' (based on absolute SSMD) to label in the plot. Default is 15.
#' @param verbose A logical indicating whether to print the computed
#'  statistics to the console. Default is \code{FALSE}.
#'
#' @return A ggplot object representing the dual flash plot for the
#' comparisons between group1 and group2.
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
#'
#' cyt_dualflashplot(
#'   data_df,
#'   group_var = "Group",
#'   group1 = "T2D",
#'   group2 = "ND",
#'   ssmd_thresh = -0.2,
#'   log2fc_thresh = 1,
#'   top_labels = 10,
#'   verbose = FALSE
#' )
#'
cyt_dualflashplot <- function(
  data,
  group_var,
  group1,
  group2,
  ssmd_thresh = 1,
  log2fc_thresh = 1,
  top_labels = 15,
  verbose = FALSE
) {
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }

  data_long <- data %>%
    tidyr::pivot_longer(
      cols = -all_of(group_var),
      names_to = "cytokine",
      values_to = "level"
    )

  stats <- data_long %>%
    dplyr::group_by(cytokine, .data[[group_var]]) %>%
    dplyr::summarise(
      .groups = "drop",
      mean = mean(level, na.rm = TRUE),
      variance = var(level, na.rm = TRUE)
    ) %>%
    tidyr::pivot_wider(
      names_from = .data[[group_var]],
      values_from = c(mean, variance)
    ) %>%
    dplyr::mutate(
      ssmd = (get(paste0("mean_", group1)) - get(paste0("mean_", group2))) /
        sqrt(
          (get(paste0("variance_", group1)) +
            get(paste0("variance_", group2))) /
            2
        ),
      log2FC = log2(
        get(paste0("mean_", group1)) /
          get(paste0("mean_", group2))
      ),
      SSMD_Category = case_when(
        abs(ssmd) >= 1 ~ "Strong Effect",
        abs(ssmd) >= 0.5 ~ "Moderate Effect",
        TRUE ~ "Weak Effect"
      ),
      Significant = (abs(ssmd) >= ssmd_thresh) & (abs(log2FC) >= log2fc_thresh)
    )

  p <- ggplot2::ggplot(stats, aes(x = log2FC, y = ssmd, label = cytokine)) +
    ggplot2::geom_point(aes(color = SSMD_Category, shape = Significant)) +
    ggplot2::geom_text(
      data = top_n(stats, top_labels, abs(ssmd)),
      vjust = 1.5,
      hjust = 1.1,
      check_overlap = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Strong Effect" = "red",
        "Moderate Effect" = "orange",
        "Weak Effect" = "blue"
      )
    ) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(
      xintercept = c(log2fc_thresh, -log2fc_thresh),
      linetype = "dashed",
      color = "blue"
    ) +
    ggplot2::labs(
      x = "Average log2 Fold Change",
      y = "SSMD",
      title = paste0("SSMD vs log2FC for ", group1, " vs ", group2)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      legend.background = element_rect(fill = "white", colour = "white"),
      axis.title = element_text(color = "black", size = 12, face = "bold"),
      legend.title = element_text(color = "black", size = 10, face = "bold"),
      legend.text = element_text(color = "black")
    )

  if (verbose) {
    print(as.data.frame(stats), n = nrow(stats), na.print = "", quote = FALSE)
  }
  return(p)
}
