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
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }

  data_long <- data %>%
    pivot_longer(cols = -all_of(group_var), names_to = "cytokine", values_to = "level")

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
      SSMD_Category = case_when(
        abs(ssmd) >= 1 ~ "Strong Effect",
        abs(ssmd) >= 0.5 ~ "Moderate Effect",
        TRUE ~ "Weak Effect"
      ),
      Significant = (abs(ssmd) >= ssmd_thresh) & (abs(log2FC) >= log2fc_thresh)
    )

  p <- ggplot(stats, aes(x = log2FC, y = ssmd, label = cytokine)) +
    geom_point(aes(color = SSMD_Category, shape = Significant)) +
    geom_text(data = top_n(stats, top_labels, abs(ssmd)), vjust = 1.5, hjust = 1.1, check_overlap = TRUE) +
    scale_color_manual(values = c("Strong Effect" = "red", "Moderate Effect" = "orange", "Weak Effect" = "blue")) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = c(log2fc_thresh, -log2fc_thresh), linetype = "dashed", color = "blue") +
    labs(x = "Average log2 Fold Change", y = "SSMD",
         title = paste0("SSMD vs log2FC for ", group1, " vs ", group2)) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"),
          legend.background = element_rect(fill = "white", colour = "white"),
          axis.title = element_text(color = "black", size = 12, face = "bold"),
          legend.title = element_text(color = "black", size = 10, face = "bold"),
          legend.text = element_text(color = "black"))

  print(stats, n = nrow(stats), na.print = "", quote = FALSE)
  return(p)
}
