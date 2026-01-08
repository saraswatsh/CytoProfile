#' Error-bar Plot
#'
#' @description
#' This function generates an error-bar plot to visually compare different
#' groups against a designated baseline group. It displays the central
#' tendency (mean or median) as a bar and overlays error bars to represent
#' the data's spread (e.g., standard deviation, MAD, or standard error).
#' The plot can also include p-value and effect size labels (based on SSMD),
#' presented either as symbols or numeric values, to highlight significant
#' differences and the magnitude of effects.
#' When an output filename is provided the plot is saved
#' to disk; otherwise the ggplot object is returned and drawn on the
#' current graphics device.
#'
#' @param data A data frame containing at least one numeric column
#'   and a grouping column.
#' @param group_col Character string naming the column that defines
#'   groups.  This column will be coerced to a factor.
#' @param p_lab Logical.  If `TRUE` (default) p‑value labels
#'   are displayed for group comparisons.
#' @param es_lab Logical.  If `TRUE` (default) effect‑size
#'   labels are displayed.
#' @param class_symbol Logical.  If `TRUE`, p‑value and
#'   effect‑size labels are encoded using symbols (e.g., `*`, `>>>`).
#'   If `FALSE`(default), numeric values are shown instead.
#' @param x_lab Character string for the x-axis label.  If empty a
#'   default label is generated.
#' @param y_lab Character string for the y-axis label.  If empty a
#'   default label is generated.
#' @param title Character string for the plot title.  If empty a
#'   default title is generated.
#' @param stat Character.  Central tendency statistic to use.  Choices
#'   are "mean" or "median"; default is "mean".  Added to
#'   support non‑mean summaries.
#' @param error Character.  Error measure visualized around the
#'   statistic.  Options are "se" (standard error; default),
#'   "sd" (standard deviation), "mad" (median absolute
#'   deviation) or "ci" (approximate 95 % confidence interval).
#' @param scale Character controlling data transformation before
#'   analysis.  Accepts "none" (default), "log2", "log10",
#'   "zscore" or "custom".
#' @param custom_fn A user‑supplied function applied to numeric
#'   columns when `scale = "custom"`.
#' @param method Character controlling the statistical test used for
#'   pairwise comparisons.  Options are "auto" (default; choose
#'   between t‑test and Wilcoxon based on a normality test),
#'   "ttest" or "wilcox".
#' @param p_adjust_method Character.  If non‑NULL, specifies the
#'   method used by `p.adjust()` to correct p‑values across all
#'   comparisons (e.g., "BH" for Benjamini–Hochberg).  If
#'   `NULL` (default) no adjustment is performed.
#' @param label_size Numeric.  Font size for p‑value and effect‑size
#'   labels.  Default is 4.
#' @param output_file Optional file path.  If provided, the plot is
#'   saved using `ggsave()`; otherwise the plot is returned and
#'   automatically printed.
#'
#' @return An error-bar plot (a `ggplot` object) is produced and optionally
#'   saved as a PDF. If `output_file` is specified, the function returns
#'   returns the `ggplot` object.
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @importFrom stats t.test wilcox.test shapiro.test mad sd median
#' @author Xiaohua Douglas Zhang and Shubh Saraswat
#' @export
#' @examples
#' # Basic usage with default settings
#' df <- ExampleData1[, c("Group", "CCL-20/MIP-3A", "IL-10")]
#' cyt_errbp(df, group_col = "Group")
#' # Use mean and SD, log2 transform and show significance
#' cyt_errbp(df, group_col = "Group", stat = "mean", error = "sd",
#'           scale = "log2", class_symbol = TRUE, method = "ttest")
cyt_errbp <- function(
  data,
  group_col = NULL,
  p_lab = TRUE,
  es_lab = TRUE,
  class_symbol = FALSE,
  x_lab = "",
  y_lab = "",
  title = "",
  stat = c("mean", "median"),
  error = c("se", "sd", "mad", "ci"),
  scale = c("none", "log2", "log10", "zscore", "custom"),
  custom_fn = NULL,
  method = c("auto", "ttest", "wilcox"),
  p_adjust_method = NULL,
  output_file = NULL,
  label_size = 4
) {
  names(data) <- make.names(names(data), unique = TRUE)
  # Validate grouping column
  if (
    is.null(group_col) || !is.character(group_col) || length(group_col) != 1
  ) {
    stop(
      "'group_col' must be provided as a single character string naming a column in 'data'.",
      call. = FALSE
    )
  }
  if (!group_col %in% names(data)) {
    stop(
      paste0(
        "The specified group_col '",
        group_col,
        "' was not found in the data frame."
      ),
      call. = FALSE
    )
  }
  # Match arguments
  stat <- match.arg(stat)
  error <- match.arg(error)
  scale <- match.arg(scale)
  method <- match.arg(method)
  data <- as.data.frame(data)
  # Convert grouping column to factor
  data[[group_col]] <- as.factor(data[[group_col]])
  # Identify numeric columns
  num_vars <- names(data)[sapply(data, is.numeric) & names(data) != group_col]
  if (length(num_vars) == 0) {
    stop("No numeric columns found in the data.")
  }
  # Apply scaling
  if (scale != "none") {
    data <- apply_scale(
      data,
      columns = num_vars,
      scale = scale,
      custom_fn = custom_fn
    )
  }
  # Reshape to long format
  long_df <- data %>%
    dplyr::select(dplyr::all_of(c(group_col, num_vars))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(num_vars),
      names_to = "Measure",
      values_to = "Value"
    )
  # Summary statistics
  summarised <- long_df %>%
    dplyr::group_by(.data[[group_col]], Measure) %>%
    dplyr::summarise(
      n = sum(!is.na(Value)),
      center = if (stat == "mean") {
        mean(Value, na.rm = TRUE)
      } else {
        median(Value, na.rm = TRUE)
      },
      sd = sd(Value, na.rm = TRUE),
      .groups = "drop"
    )
  # Compute spread based on error metric
  summarised <- summarised %>%
    dplyr::mutate(
      spread = dplyr::case_when(
        error == "sd" ~ sd,
        error == "se" ~ sd / sqrt(n),
        error == "mad" ~ {
          # compute MAD per group and measure
          mapply(
            function(g, m) {
              mad(
                long_df$Value[long_df$Measure == m & long_df[[group_col]] == g],
                na.rm = TRUE
              )
            },
            .data[[group_col]],
            Measure
          )
        },
        error == "ci" ~ (sd / sqrt(n)) * 1.96,
        TRUE ~ sd / sqrt(n)
      )
    )
  # Determine baseline group
  group_levels <- levels(data[[group_col]])
  baseline <- group_levels[1]
  # Initialise p‑value and effect size columns
  summarised <- summarised %>%
    dplyr::mutate(P_value = NA_real_, EffectSize = NA_real_)
  # Compute p‑values and effect sizes
  if (p_lab || es_lab) {
    for (m in unique(summarised$Measure)) {
      baseline_values <- long_df$Value[
        long_df$Measure == m & long_df[[group_col]] == baseline
      ]
      base_mean <- mean(baseline_values, na.rm = TRUE)
      base_sd <- sd(baseline_values, na.rm = TRUE)
      for (g in setdiff(group_levels, baseline)) {
        idx <- summarised$Measure == m & summarised[[group_col]] == g
        grp_values <- long_df$Value[
          long_df$Measure == m & long_df[[group_col]] == g
        ]
        # Choose test
        use_test <- method
        if (method == "auto") {
          combined <- c(baseline_values, grp_values)
          p_norm <- tryCatch(
            shapiro.test(combined)$p.value,
            error = function(e) NA_real_
          )
          use_test <- if (!is.na(p_norm) && p_norm < 0.05) "wilcox" else "ttest"
        }
        if (use_test == "ttest") {
          tt <- t.test(grp_values, baseline_values)
          p_val <- tt$p.value
          grp_mean <- mean(grp_values, na.rm = TRUE)
          grp_sd <- sd(grp_values, na.rm = TRUE)
          pooled_sd <- sqrt(
            ((length(grp_values) - 1) *
              grp_sd^2 +
              (length(baseline_values) - 1) * base_sd^2) /
              (length(grp_values) + length(baseline_values) - 2)
          )
          eff <- ifelse(
            pooled_sd == 0,
            NA_real_,
            (grp_mean - base_mean) / pooled_sd
          )
        } else {
          wt <- wilcox.test(grp_values, baseline_values, exact = FALSE)
          p_val <- wt$p.value
          u <- wt$statistic
          eff <- (u / (length(grp_values) * length(baseline_values))) * 2 - 1
        }
        summarised$P_value[idx] <- p_val
        summarised$EffectSize[idx] <- eff
      }
    }
    # Adjust p‑values if a method is specified
    if (!is.null(p_adjust_method)) {
      summarised$P_adj <- adjust_p(summarised$P_value, method = p_adjust_method)
    } else {
      summarised$P_adj <- summarised$P_value
    }
  } else {
    summarised$P_adj <- summarised$P_value
  }
  # Compute label positions
  y_range <- diff(range(summarised$center + summarised$spread, na.rm = TRUE))
  summarised <- summarised %>%
    dplyr::mutate(
      p_y = center + spread + 0.05 * y_range,
      es_y = center + spread + 0.15 * y_range
    )
  # Create labels for p-values and effect sizes.  We compute both symbol and
  # numeric representations using the original significance and effect size
  # classification functions and combine them based on the values of
  # `p_lab`, `es_lab` and `class_symbol`.  When both a label flag and
  # class_symbol are TRUE, both the symbol and the numeric value are
  # displayed on separate lines.  When only one of p_lab or es_lab is TRUE,
  # the corresponding numeric label is shown; when only class_symbol is
  # TRUE, only the symbol is shown.  When neither a label flag nor
  # class_symbol applies, no label is drawn.
  # Initialise label columns
  summarised$p_label <- NA_character_
  summarised$es_label <- NA_character_
  # Define significance and effect size classification functions (copied
  # from the original implementation).
  significance_mark_fn <- function(p_value) {
    if (is.na(p_value)) {
      return(NA_character_)
    }
    if (p_value <= 0.00001) {
      return("*****")
    }
    if (p_value <= 0.0001) {
      return("****")
    }
    if (p_value <= 0.001) {
      return("***")
    }
    if (p_value <= 0.01) {
      return("**")
    }
    if (p_value <= 0.05) {
      return("*")
    }
    return("")
  }
  effect_size_mark_fn <- function(es) {
    if (is.na(es)) {
      return(NA_character_)
    }
    if (es >= 5) {
      return(">>>>>")
    }
    if (es >= 3) {
      return(">>>>")
    }
    if (es >= 1.645) {
      return(">>>")
    }
    if (es >= 1) {
      return(">>")
    }
    if (es > 0.25) {
      return(">")
    }
    if (es >= -0.25) {
      return(" ")
    }
    if (es > -1) {
      return("<")
    }
    if (es > -1.645) {
      return("<<")
    }
    if (es > -3) {
      return("<<<")
    }
    if (es > -5) {
      return("<<<<")
    }
    return("<<<<<")
  }
  # Compute star-coded significance and arrow-coded effect size symbols
  p_symbol <- sapply(summarised$P_adj, significance_mark_fn)
  es_symbol <- sapply(summarised$EffectSize, effect_size_mark_fn)
  # Compute numeric labels
  p_numeric <- ifelse(
    is.na(summarised$P_adj),
    NA_character_,
    ifelse(
      summarised$P_adj > 0.001,
      sprintf("p=%.3f", summarised$P_adj),
      paste0("p=", formatC(summarised$P_adj, format = "e", digits = 1))
    )
  )
  es_numeric <- ifelse(
    is.na(summarised$EffectSize),
    NA_character_,
    sprintf("es=%.3f", summarised$EffectSize)
  )
  # Combine according to user options for p-values
  summarised$p_label <- mapply(
    function(sym, num) {
      # When both p_lab and class_symbol are TRUE, display the symbol before
      # the numeric p-value, separated by a space.  When only one option is
      # enabled, show the corresponding label.  Empty symbols ("" or " ")
      # are treated as absent.
      if (p_lab && class_symbol) {
        # Both requested
        if (!is.na(sym) && sym != "" && sym != " " && !is.na(num)) {
          return(paste(sym, num))
        }
        # If symbol is blank or NA but numeric exists
        if (!is.na(num)) {
          return(num)
        }
        # No numeric or symbol
        return(NA_character_)
      } else if (p_lab) {
        # Only numeric requested
        return(ifelse(is.na(num), NA_character_, num))
      } else if (class_symbol) {
        # Only symbol requested
        return(ifelse(
          is.na(sym) || sym == "" || sym == " ",
          NA_character_,
          sym
        ))
      }
      return(NA_character_)
    },
    p_symbol,
    p_numeric,
    SIMPLIFY = FALSE
  )
  summarised$p_label <- unlist(summarised$p_label)
  # Combine according to user options for effect sizes
  summarised$es_label <- mapply(
    function(sym, num) {
      if (es_lab && class_symbol) {
        if (!is.na(sym) && sym != "" && sym != " " && !is.na(num)) {
          return(paste(sym, num))
        }
        if (!is.na(num)) {
          return(num)
        }
        return(NA_character_)
      } else if (es_lab) {
        return(ifelse(is.na(num), NA_character_, num))
      } else if (class_symbol) {
        return(ifelse(
          is.na(sym) || sym == "" || sym == " ",
          NA_character_,
          sym
        ))
      }
      return(NA_character_)
    },
    es_symbol,
    es_numeric,
    SIMPLIFY = FALSE
  )
  summarised$es_label <- unlist(summarised$es_label)
  # Default labels
  if (x_lab == "") {
    x_lab <- group_col
  }
  if (y_lab == "") {
    y_lab <- paste(stat, "value")
  }
  if (title == "") {
    title <- paste(
      "Error Bar Plots for",
      paste(unique(summarised$Measure), collapse = ", ")
    )
  }
  # Build plot
  p <- ggplot2::ggplot(
    summarised,
    ggplot2::aes(x = .data[[group_col]], y = center)
  ) +
    ggplot2::geom_col(fill = "grey80") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = center - spread, ymax = center + spread),
      width = 0.2
    ) +
    ggplot2::facet_wrap(~Measure, scales = "free_y") +
    ggplot2::labs(x = x_lab, y = y_lab, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  # Add p‑value labels
  if (p_lab) {
    p <- p +
      ggplot2::geom_text(
        data = summarised %>%
          dplyr::filter(.data[[group_col]] != baseline & !is.na(p_label)),
        ggplot2::aes(x = .data[[group_col]], y = p_y, label = p_label),
        size = label_size,
        colour = "black",
        vjust = 0,
        na.rm = TRUE
      )
  }
  # Add effect size labels
  if (es_lab) {
    p <- p +
      ggplot2::geom_text(
        data = summarised %>%
          dplyr::filter(.data[[group_col]] != baseline & !is.na(es_label)),
        ggplot2::aes(x = .data[[group_col]], y = es_y, label = es_label),
        size = label_size,
        colour = "black",
        vjust = 0,
        na.rm = TRUE
      )
  }
  # Save or display
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 8, height = 5)
    return(invisible(p))
  }
  print(p)
  invisible(p)
}
