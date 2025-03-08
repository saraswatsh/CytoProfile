#' Volcano Plot.
#'
#' @param data A matrix or data frame containing the data to be analyzed.
#' @param group_col A character string specifying the column name used for
#'   comparisons (e.g., group, treatment, or stimulation).
#' @param cond1 A character string specifying the name of the first condition
#'   for comparison. Default is \code{NULL}.
#' @param cond2 A character string specifying the name of the second condition
#'   for comparison. Default is \code{NULL}.
#' @param fold_change_thresh A numeric threshold for the fold change.
#'   Default is \code{2}.
#' @param p_value_thresh A numeric threshold for the p-value.
#'   Default is \code{0.05}.
#' @param top_labels An integer specifying the number of top variables to label
#'   on the plot. Default is \code{10}.
#'
#' @description This function subsets the numeric columns from the input data
#'   and compares them based on a selected grouping column. It computes the fold
#'   changes (as the ratio of means) and associated p-values (using two-sample
#'   t-tests) for each numeric variable between two groups. The results are
#'   log2-transformed (for fold change) and -log10-transformed (for p-values)
#'   to generate a volcano plot.
#'
#' @note If \code{cond1} and \code{cond2} are not provided, the function
#'   automatically generates all possible pairwise combinations of groups from
#'   the specified \code{group_col} for comparisons.
#'
#' @return A list of volcano plots (as \code{ggplot} objects) for each pairwise
#'   comparison. Additionally, the function prints the data frame used for
#'   plotting (excluding the significance column) from the final comparison.
#'
#' @export
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#'
#' @examples
#' # Loading data
#' data_df <- cytodata[,-4]
#'
#' volc_plot <- cyt_volc(data_df, "Group", cond1 = "T2D", cond2 = "ND",
#' fold_change_thresh = 2.0, top_labels= 15)
#' print(volc_plot$`T2D vs ND`)

cyt_volc <- function(data, group_col, cond1 = NULL, cond2 = NULL,
                     fold_change_thresh = 2,
                     p_value_thresh = 0.05, top_labels = 10) {
  # Determine the pairs of conditions to compare
  if (!is.null(cond1) && !is.null(cond2)) {
    condition_pairs <- list(c(cond1, cond2))
  } else {
    # Extract unique conditions and generate all pairs
    conditions <- unique(data[[group_col]])
    condition_pairs <- combn(conditions, 2, simplify = FALSE)
  }

  # Initialize a list to store plots
  plots <- list()

  for (pair in condition_pairs) {
    cond1 <- pair[1]
    cond2 <- pair[2]

    # Separate the data based on conditions
    data_cond1 <- data[data[, group_col] == cond1, ]
    data_cond2 <- data[data[, group_col] == cond2, ]

    # Ensure only numeric columns are used for calculations
    numeric_columns <- sapply(data, is.numeric)
    data_ncond1 <- data_cond1[, numeric_columns]
    data_ncond2 <- data_cond2[, numeric_columns]

    # Calculate means for each cytokine in each group
    means_cond1 <- colMeans(data_ncond1, na.rm = TRUE)
    means_cond2 <- colMeans(data_ncond2, na.rm = TRUE)

    # Calculate fold changes and p-values
    fold_changes <- means_cond2 / means_cond1
    p_values <- mapply(function(x, y) t.test(x, y)$p.value,
                       data_ncond1, data_ncond2)

    # Log2 transform of fold changes and -log10 of p-values
    fc_log <- log2(fold_changes)
    p_log <- -log10(p_values)

    # Create dataframe for plotting (convert column names to snake_case)
    plot_data <- data.frame(
      cytokine = names(fold_changes),
      fc_log = fc_log,
      p_log = p_log,
      significant = (abs(fc_log) >= log2(fold_change_thresh)) &
        (p_log >= -log10(p_value_thresh))
    )

    # Order data by significance and p_log, then add labels for top points
    plot_data <- plot_data %>%
      arrange(desc(significant), desc(p_log)) %>%
      mutate(label = ifelse(row_number() <= top_labels,
                            as.character(cytokine), ""))

    # Create the volcano plot with labels for top significant points
    volcano_plot <- ggplot(plot_data, aes(x = fc_log, y = p_log, label = label,
                                          color = significant)) +
      geom_point(alpha = 1, size = 2) +
      geom_vline(xintercept = c(log2(fold_change_thresh),
                                -log2(fold_change_thresh)),
                 linetype = "dashed", color = "blue") +
      geom_hline(yintercept = -log10(p_value_thresh),
                 linetype = "dashed", color = "blue") +
      ggrepel::geom_text_repel(aes(label = label), size = 3,
                               vjust = 1.5, hjust = 0.5,
                               show.legend = FALSE) +
      scale_color_manual(values = c("grey2", "red")) +
      labs(x = "Log2 Fold Change", y = "-Log10 P-Value",
           title = paste("Volcano Plot of Cytokine Levels:",
                         cond1, "vs", cond2)) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        legend.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(color = "black", size = 12, face = "bold"),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black")
      )
    # Store the plot in the list
    plots[[paste(cond1, "vs", cond2)]] <- volcano_plot
  }

  # Print the final plot data (excluding the label column)
  print(plot_data[, -which(names(plot_data) == "label")],
        n = nrow(plot_data), na.print = "", quote = FALSE)
  return(plots)
}
