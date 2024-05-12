#' Dual flash plot.
#'
#' @param data A data frame to be entered.
#' @param group_var Grouping column to be used from the data frame.
#' @param group1 Name of the first group to be compared.
#' @param group2 Name of the second group to be compared.
#' @param ssmd_thresh Threshold value for SSMD. Default set to 1.
#' @param log2fc_thresh Threshold value for log2 fold change. Default set to 1.
#' @param top_labels Number of labels to be print for top variables. Default set to 15.
#'
#' @return Prints dual flash plot.
#' @export
#'
#' @examples
#' cyt.dualflashplot(cytdata.df[,-c(1,3)], group_var = "Group", group1 = "ND", group2 = "T2D", ssmd_thresh = -0.2, log2fc_thresh = 1, top_labels = 10)

cyt.dualflashplot <- function(data, group_var, group1, group2, ssmd_thresh = 1, log2fc_thresh = 1, top_labels = 15) {
  # Check that the data is a dataframe
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }

  # Reshape data to long format
  data_long <- data %>%
    pivot_longer(cols = -all_of(group_var), names_to = "cytokine", values_to = "level")

  # Calculate means, variances, ssmd, and log2 fold change
  stats <- data_long %>%
    group_by(cytokine, .data[[group_var]]) %>%
    summarise(.groups = "drop",
              mean = mean(level, na.rm = TRUE),
              variance = var(level, na.rm = TRUE)) %>%
    pivot_wider(names_from = .data[[group_var]], values_from = c(mean, variance)) %>%
    mutate(
      ssmd = (get(paste0("mean_", group1)) - get(paste0("mean_", group2))) /
        sqrt((get(paste0("variance_", group1)) + get(paste0("variance_", group2))) / 2),
      log2FC = log2(get(paste0("mean_", group1)) / get(paste0("mean_", group2))),
      Significant = (abs(ssmd) >= ssmd_thresh) & (abs(log2FC) >= log2fc_thresh)
    ) %>%
    ungroup()

  top_n <- top_labels  # Number of cytokines to label based on highest SSMD

  p <- ggplot(stats, aes(x = log2FC, y = ssmd, label = cytokine)) +
    geom_point(aes(color = Significant)) +
    geom_text(data = top_n(stats %>% arrange(desc(abs(ssmd))), top_n), check_overlap = TRUE, nudge_y = 0.03) +
    theme_minimal() +
    labs(x = "Average log2 Fold Change", y = "SSMD", title = paste0("SSMD vs log2FC for ", group1, " vs ", group2)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = log2fc_thresh, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = -log2fc_thresh, linetype = "dashed", color = "blue") +
    scale_color_manual(values = c("grey2", "red"))+
    theme(panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"),  # Ensure plot background is white
          legend.background = element_rect(fill = "white", colour = "white"),  # Ensure legend background is white
          axis.title = element_text(color = "black", size = 12, face = "bold"),  # Customize axis titles
          legend.title = element_text(color = "black", size = 10, face = "bold"),  # Customize legend title
          legend.text = element_text(color = "black"))

  # Return the plot
  print(stats, n = nrow(stats), na.print = "", quote = FALSE)
  return(p)
}
