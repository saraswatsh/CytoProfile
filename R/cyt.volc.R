#' Volcano Plot
#'
#' @param x.df A matrix or data frame
#' @param group_col Column used for comparison (i.e. group, treatment, or stimulation).
#' @param cond1 A string of a group name to be chosen. Default set to NULL.
#' @param cond2 A string of a group name to be chosen. Default set to NULL.
#' @param fold_change_thresh Threshold value for fold change. Default set to 2.
#' @param p_value_thresh Threshold value for p-value. Default set to 0.05.
#' @param top_labels Number of labels to be print for top variables. Default set to 10
#' @description
#' This function takes a matrix or data frame and subsets the numeric columns or variables
#' to compare against each other based on a selected column of choice that defines a category
#' to draw a volcano plot.
#' @note
#' If cond1 and cond2 are set to null, by default the function creates combinations of pairs
#' to be used for comparisons.
#' @return Prints volcano plot.
#' @export
#'
#' @examples
#' cyt.volc(cytodata,group_col = "Group")
#' cyt.volc(cytodata, group_col = "Group", fold_change_thresh = 2, top_labels = 15)

cyt.volc <- function(data, group_col, cond1 = NULL, cond2 = NULL, fold_change_thresh = 2, p_value_thresh = 0.05, top_labels = 10) {

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
    p_values <- mapply(function(x, y) t.test(x, y)$p.value, data_ncond1, data_ncond2)

    # Log2 transform of fold changes and -log10 of p-values
    fc_log <- log2(fold_changes)
    p_log <- -log10(p_values)

    # Create dataframe for plotting
    plot_data <- data.frame(Cytokine = names(fold_changes),
                            FC_Log = fc_log,
                            P_Log = p_log,
                            Significant = (abs(fc_log) >= log2(fold_change_thresh)) & (p_log >= -log10(p_value_thresh)))

    # Order data by significance and P_Log
    plot_data <- plot_data %>%
      arrange(desc(Significant), desc(P_Log)) %>%
      mutate(Label = ifelse(row_number() <= top_labels, as.character(Cytokine), ""))

    # Create the volcano plot with labels for top significant points
    volcano_plot <- ggplot(plot_data, aes(x = FC_Log, y = P_Log, label = Label, color = Significant)) +
      geom_point(alpha = 1, size = 2) +
      geom_vline(xintercept = c(log2(fold_change_thresh), -log2(fold_change_thresh)), linetype = "dashed", color = "blue") +
      geom_hline(yintercept = -log10(p_value_thresh), linetype = "dashed", color = "blue") +
      ggrepel::geom_text_repel(aes(label = Label), size = 3, vjust = 1.5, hjust = .5, show.legend = FALSE) +
      scale_color_manual(values = c("grey2", "red")) +
      labs(x = "Log2 Fold Change", y = "-Log10 P-Value", title = paste("Volcano Plot of Cytokine Levels:", cond1, "vs", cond2)) +
      theme_minimal()+
      theme(panel.background = element_rect(fill = "white", colour = "white"),
            plot.background = element_rect(fill = "white", colour = "white"),  # Ensure plot background is white
            legend.background = element_rect(fill = "white", colour = "white"),  # Ensure legend background is white
            axis.title = element_text(color = "black", size = 12, face = "bold"),  # Customize axis titles
            legend.title = element_text(color = "black", size = 10, face = "bold"),  # Customize legend title
            legend.text = element_text(color = "black"))
    # Store the plot in the list
    plots[[paste(cond1, "vs", cond2)]] <- volcano_plot
  }
  # Return the list of plots
  print(plot_data[,-5], n = nrow(plot_data),na.print = "", quote = FALSE)
  return(plots)
}
