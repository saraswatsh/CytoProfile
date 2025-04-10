#' Error-bar Plot.
#'
#' This function generates an error-bar plot to visually compare different
#' groups against a designated baseline group. It displays the central
#' tendency (mean or median) as a bar and overlays error bars to represent
#' the data's spread (e.g., standard deviation, MAD, or standard error).
#' The plot can also include p-value and effect size labels (based on SSMD),
#' presented either as symbols or numeric values, to highlight significant
#' differences and the magnitude of effects.
#'
#' @param data A data frame containing the data for each group. It should
#'   include at least one numeric column for the measurements and a column
#'   specifying the group membership.
#' @param group_col Character. The name of the column in `data` that
#'   specifies the group membership.
#' @param p_lab Logical. If `TRUE`, p-values are displayed on the plot.
#'   Default is `FALSE`.
#' @param es_lab Logical. If `TRUE`, effect sizes (SSMD) are displayed on
#'   the plot. Default is `FALSE`.
#' @param class_symbol Logical. If `TRUE`, significance and effect size are
#'   represented using symbolic notation (e.g., *, **, >, <<). If `FALSE`,
#'   numeric values are used. Default is `TRUE`.
#' @param x_lab Character. Label for the x-axis. If not provided, defaults
#'   to the name of the `group_col` or "Group" if `group_col` is `NULL`.
#' @param y_lab Character. Label for the y-axis. If not provided, defaults
#'   to "Value".
#' @param title Character. Title of the plot. If not provided, a default
#'   title is generated based on the measured variables.
#' @param log2 Logical. If `TRUE`, a log2 transformation (with a +1 offset)
#'   is applied to all numeric columns before analysis. Default is `FALSE`.
#' @param output_file Character. The file path to save the plot as a PDF.
#'   If `NULL`, the plot is displayed but not saved. Default is `NULL`.
#'
#' @return An error-bar plot (a `ggplot` object) is produced and optionally
#'   saved as a PDF. If `output_file` is specified, the function returns
#'   returns the `ggplot` object.
#' 
#' @details
#'   The function performs the following steps:
#'   \enumerate{
#'     \item Optionally applies a log2 transformation to numeric data.
#'     \item Determines the baseline group (the first level of `group_col`).
#'     \item Calculates summary statistics (sample size, mean, standard
#'           deviation) for each group and each numeric variable.
#'     \item Performs t-tests to compare each group against the baseline
#'           for each numeric variable.
#'     \item Computes effect sizes (SSMD) for each group compared to the
#'           baseline.
#'     \item Generates a faceted error-bar plot, with one facet per
#'           numeric variable.
#'     \item Optionally adds p-value and effect size labels to the plot.
#'     \item Optionally saves the plot as a PDF.
#'   }
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export
#' @examples
#' data <- ExampleData1
#' 
#' cyt_errbp(data[,c("Group", "CCL.20.MIP.3A", "IL.10")], group_col = "Group", 
#' p_lab = TRUE, es_lab = TRUE, class_symbol = TRUE, x_lab = "Cytokines", 
#' y_lab = "Concentrations in log2 scale", log2 = TRUE)

cyt_errbp <- function(data, group_col = NULL,
  p_lab = FALSE, es_lab = FALSE,
  class_symbol = TRUE,
  x_lab = "", y_lab = "", title = "",
  log2 = FALSE,
  output_file = NULL) {

if (missing(group_col) || is.null(group_col)) {
    stop("The 'group_col' argument is required. Please provide the name of the column containing group information.", call. = FALSE)
}
if (!is.character(group_col) || length(group_col) != 1) {
    stop("'group_col' must be a single character string naming a column in 'data'.", call. = FALSE)
}
if (!group_col %in% names(data)) {
    stop(paste0("The specified group_col '", group_col, "' was not found in the data frame."), call. = FALSE)
  }
  
# If log2 transformation is requested, transform all numeric columns (or a subset if desired)
if (log2) {
numeric_cols <- sapply(data, is.numeric)
data[numeric_cols] <- lapply(data[numeric_cols], function(x) log2(x))
}
  
# Convert the specified group column to a factor
data[[group_col]] <- factor(data[[group_col]])
  
# Identify all numeric columns, excluding the grouping column.
num_vars <- names(data)[sapply(data, is.numeric) & names(data) != group_col]
if (length(num_vars) == 0)
stop("No numeric columns found in the data.")
if (y_lab == "") y_lab <- "Value"

# Reshape the data into long format (one row per numeric measurement)
long_df <- data %>%
dplyr::select(dplyr::all_of(c(group_col, num_vars))) %>%
tidyr::pivot_longer(cols = dplyr::all_of(num_vars),
names_to = "Measure",
values_to = "Value")

# Calculate summary statistics (sample size, mean, standard deviation) per group and per measure.
metrics <- long_df %>%
dplyr::group_by(.data[[group_col]], Measure) %>%
dplyr::summarize(
n = sum(!is.na(Value)),
center = mean(Value, na.rm = TRUE),
sd = sd(Value, na.rm = TRUE),
.groups = "drop"
) %>%
dplyr::mutate(spread = sd / sqrt(n))

# Determine the baseline group: using the first level of the grouping variable.
group_levels <- levels(data[[group_col]])
baseline <- group_levels[1]

# Initialize columns for p-value and effect size.
metrics <- metrics %>%
dplyr::mutate(p.value = NA_real_,
effect.size = NA_real_)


# For each measure, perform a t-test comparing each group against the baseline and compute an effect size.
unique_measures <- unique(metrics$Measure)
for (m in unique_measures) {
baseline_row <- metrics %>% dplyr::filter(Measure == m, .data[[group_col]] == baseline)
if (nrow(baseline_row) == 0) next
base_mean <- baseline_row$center
base_sd <- baseline_row$sd
base_n <- baseline_row$n

# Get indices for non-baseline groups in the measure m.
non_base_idx <- which(metrics$Measure == m & metrics[[group_col]] != baseline)
for (i in non_base_idx) {
current_group <- metrics[[group_col]][i]
baseline_data <- long_df %>% 
dplyr::filter(Measure == m, .data[[group_col]] == baseline) %>% 
dplyr::pull(Value)
grp_data <- long_df %>% 
dplyr::filter(Measure == m, .data[[group_col]] == current_group) %>% 
dplyr::pull(Value)

tt <- t.test(grp_data, baseline_data)
metrics$p.value[i] <- tt$p.value

grp_mean <- metrics$center[i]
grp_sd <- metrics$sd[i]
grp_n <- metrics$n[i]
pooled_sd <- sqrt(((base_n - 1) * base_sd^2 + (grp_n - 1) * grp_sd^2) / (base_n + grp_n - 2))
d <- (grp_mean - base_mean) / pooled_sd
metrics$effect.size[i] <- d
}
}

# Create labels for p-values and effect sizes according to class_symbol.
if (p_lab) {
if (class_symbol) {
significance_mark_fn <- function(p_value) {
if (is.na(p_value)) return(NA_character_)
if (p_value <= 0.00001) return("*****")
if (p_value <= 0.0001)  return("****")
if (p_value <= 0.001)   return("***")
if (p_value <= 0.01)    return("**")
if (p_value <= 0.05)    return("*")
return("")
}
metrics <- metrics %>%
dplyr::mutate(p_label = sapply(p.value, significance_mark_fn))
} else {
metrics <- metrics %>%
dplyr::mutate(p_label = paste0("p=", ifelse(p.value > 0.001, round(p.value, 3),
                   formatC(p.value, format = "e", digits = 1))))
}
}

if (es_lab) {
if (class_symbol) {
effect_size_mark_fn <- function(es) {
if (is.na(es)) return(NA_character_)
if (es >= 5) return(">>>>>")
if (es >= 3) return(">>>>")
if (es >= 1.645) return(">>>")
if (es >= 1) return(">>")
if (es > 0.25) return(">")
if (es >= -0.25) return(" ")
if (es > -1) return("<")
if (es > -1.645) return("<<")
if (es > -3) return("<<<")
if (es > -5) return("<<<<")
return("<<<<<")
}
metrics <- metrics %>%
dplyr::mutate(es_label = sapply(effect.size, effect_size_mark_fn))
} else {
metrics <- metrics %>%
dplyr::mutate(es_label = round(effect.size, 3))
}
}

# Compute a rough y-range to position annotations.
y_range <- diff(range(metrics$center + metrics$spread, na.rm = TRUE))
metrics <- metrics %>%
  dplyr::mutate(
    p_text_y = center + ifelse(center >= 0, spread + y_range/20, -spread - y_range/20),
    es_text_y = center + ifelse(center >= 0, spread + y_range/4, -spread - y_range/4)  # increased from y_range/8 to y_range/6
  )

# Set default axis labels and title if not provided.
if (x_lab == "") x_lab <- group_col
if (title == "") title <- paste("Error Bar Plots for", paste(unique_measures, collapse = ", "))

# Build the faceted ggplot (one facet per numeric measure).
p <- ggplot2::ggplot(metrics, aes_string(x = group_col, y = "center")) +
  ggplot2::geom_bar(stat = "identity", fill = "gray", width = 0.7) +
  ggplot2::geom_errorbar(aes(ymin = center - spread, ymax = center + spread), width = 0.2) +
  ggplot2::facet_wrap(~ Measure, scales = "free_y") +
  ggplot2::labs(x = x_lab, y = y_lab, title = title) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.background = element_rect(fill = "white", colour = "white"),
                plot.background = element_rect(fill = "white", colour = "white"),
                legend.background = element_rect(fill = "white", colour = "white"),
                axis.title = element_text(color = "black", size = 12, face = "bold"),
                legend.title = element_text(color = "black", size = 10, face = "bold"),
                legend.text = element_text(color = "black"))

# Add text annotations if requested.
if (p_lab) {
p <- p + ggplot2::geom_text(
  data = metrics %>% dplyr::filter(.data[[group_col]] != baseline),
  aes_string(x = group_col, y = "p_text_y", label = "p_label"),
  size = 4, vjust = 0
  )
}
if (es_lab) {
    p <- p + ggplot2::geom_text(
    data = metrics %>% dplyr::filter(.data[[group_col]] != baseline),
    ggplot2::aes_string(x = group_col, y = "es_text_y", label = "es_label"),
    size = 4, vjust = 0
    )
}

if (!is.null(output_file)) {
    pdf(file = output_file, width = 7, height = 5)
    print(p)
    dev.off()
  } else {
    return(p)
  }
}
  