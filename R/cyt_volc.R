#' Volcano plot
#'
#' @param data A data frame containing numeric variables and a
#'   grouping column.
#' @param group_col Character.  Name of the grouping column.
#' @param cond1 Character strings specifying the levels of
#'   `group_col` to compare.  If either is `NULL`, all pairwise
#'   combinations of conditions are used.
#' @param cond2 Character strings specifying the levels of
#'   `group_col` to compare.  If either is `NULL`, all pairwise
#'   combinations of conditions are used.
#' @param fold_change_thresh Numeric.  Threshold for absolute fold
#'   change (in original scale).  Default is 2.
#' @param p_value_thresh Numeric.  Threshold for the p‑value (raw or
#'   adjusted).  Default is 0.05.
#' @param top_labels Integer.  Number of top points to label in each
#'   plot.  Default is 10.
#' @param method Character.  Statistical test to use.  "ttest" (default)
#'   uses two‑sample t‑tests; "wilcox" uses Wilcoxon rank‑sum tests.
#' @param p_adjust_method Character or `NULL`.  Method to adjust
#'   p‑values across variables within each comparison (e.g., "BH").
#'   If `NULL` (default) no adjustment is performed. See \code{\link[stats]{p.adjust}}
#'   for details.
#' @param add_effect Logical.  If `TRUE`, effect sizes are computed and
#'   returned in the results (Cohen's d for t‑tests; rank‑biserial for
#'   Wilcoxon).  Default is `FALSE`.
#' @param verbose Logical.  If `TRUE`, prints the data frame used for
#'   the final comparison without the label column.  Default is
#'   `FALSE`.
#'
#' @description This function subsets the numeric columns from the input data
#'   and compares them based on a selected grouping column. It computes the fold
#'   changes (as the ratio of means) and associated p-values for each numeric variable
#'   between two groups. The results are log2-transformed (for fold change) and
#'   -log10-transformed (for p-values) to generate a volcano plot. Additionally,
#'   there is a choice between t‑tests and Wilcoxon rank‑sum tests and adjusting
#'   p‑values for multiple comparisons
#'
#' @note If \code{cond1} and \code{cond2} are not provided, the function
#'   automatically generates all possible pairwise combinations of groups from
#'   the specified \code{group_col} for comparisons.
#'
#' @return A list of ggplot objects (one per comparison).  Each plot
#'   visualizes log2 fold change on the x‑axis and –log10 of the
#'   (adjusted) p‑value on the y‑axis.  The underlying data used to
#'   construct the final plot are printed when `verbose = TRUE`.
#'
#' @author Xiaohua Douglas Zhang and Shubh Saraswat
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr arrange mutate desc row_number
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats t.test wilcox.test var na.omit
#'
#' @examples
#' # Loading data
#' data_df <- ExampleData1[,-c(2:3)]
#'
#' cyt_volc(data_df, "Group", cond1 = "T2D", cond2 = "ND", fold_change_thresh = 2.0, top_labels= 15)

cyt_volc <- function(
  data,
  group_col,
  cond1 = NULL,
  cond2 = NULL,
  fold_change_thresh = 2,
  p_value_thresh = 0.05,
  top_labels = 10,
  method = c("ttest", "wilcox"),
  p_adjust_method = NULL,
  add_effect = FALSE,
  verbose = FALSE
) {
  method <- match.arg(method)
  data <- as.data.frame(data)
  if (!group_col %in% names(data)) {
    stop(paste0("group_col '", group_col, "' not found."))
  }
  group <- factor(data[[group_col]])
  # Determine comparison pairs
  if (!is.null(cond1) && !is.null(cond2)) {
    condition_pairs <- list(c(cond1, cond2))
  } else {
    conds <- unique(group)
    if (length(conds) < 2) {
      stop("At least two groups are required for comparison.")
    }
    condition_pairs <- combn(as.character(conds), 2, simplify = FALSE)
  }
  # Identify numeric variables
  numeric_columns <- sapply(data, is.numeric)
  num_vars <- names(data)[numeric_columns & names(data) != group_col]
  if (length(num_vars) == 0) {
    stop("No numeric variables found for volcano plot.")
  }
  plots <- list()
  # Loop through pairs
  for (pair in condition_pairs) {
    c1 <- pair[1]
    c2 <- pair[2]
    # Subset data
    d1 <- data[data[[group_col]] == c1, , drop = FALSE]
    d2 <- data[data[[group_col]] == c2, , drop = FALSE]
    # Compute statistics for each numeric variable
    res <- lapply(num_vars, function(var) {
      x1 <- d1[[var]]
      x2 <- d2[[var]]
      # Fold change (ratio of means) and log2FC
      fc <- mean(x2, na.rm = TRUE) / (mean(x1, na.rm = TRUE) + 1e-12)
      log2fc <- log2(fc + 1e-12)
      # Choose test
      if (method == "ttest") {
        tst <- t.test(x1, x2)
        pval <- tst$p.value
        eff <- if (add_effect) {
          m1 <- mean(x1, na.rm = TRUE)
          m2 <- mean(x2, na.rm = TRUE)
          s1 <- var(x1, na.rm = TRUE)
          s2 <- var(x2, na.rm = TRUE)
          n1 <- length(na.omit(x1))
          n2 <- length(na.omit(x2))
          pooled_sd <- sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
          if (is.na(pooled_sd) || pooled_sd == 0) {
            NA_real_
          } else {
            (m2 - m1) / pooled_sd
          }
        } else {
          NA_real_
        }
      } else {
        tst <- wilcox.test(x1, x2, exact = FALSE)
        pval <- tst$p.value
        eff <- if (add_effect) {
          u <- tst$statistic
          n1 <- length(na.omit(x1))
          n2 <- length(na.omit(x2))
          (u / (n1 * n2)) * 2 - 1
        } else {
          NA_real_
        }
      }
      data.frame(
        cytokine = var,
        log2FC = log2fc,
        P_value = pval,
        EffectSize = eff
      )
    })
    df <- do.call(rbind, res)
    # Adjust p‑values if requested
    if (!is.null(p_adjust_method)) {
      df$P_adj <- adjust_p(df$P_value, method = p_adjust_method)
    } else {
      df$P_adj <- df$P_value
    }
    # Compute –log10 p
    df$minusLog10P <- -log10(df$P_adj + 1e-300)
    # Determine significance by thresholds (fold change on original scale)
    df$Significant <- with(
      df,
      P_adj <= p_value_thresh & abs(log2FC) >= log2(fold_change_thresh)
    )
    # Sort and label top points
    df <- df %>%
      dplyr::arrange(dplyr::desc(Significant), dplyr::desc(minusLog10P)) %>%
      dplyr::mutate(
        label = ifelse(
          dplyr::row_number() <= top_labels,
          as.character(cytokine),
          ""
        )
      )
    # Create plot
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = log2FC,
        y = minusLog10P,
        label = label,
        colour = Significant
      )
    ) +
      ggplot2::geom_point(alpha = 1, size = 2) +
      ggplot2::geom_vline(
        xintercept = c(log2(fold_change_thresh), -log2(fold_change_thresh)),
        linetype = "dashed",
        color = "blue"
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(p_value_thresh),
        linetype = "dashed",
        color = "blue"
      ) +
      ggrepel::geom_text_repel(
        aes(label = label),
        size = 3,
        vjust = 1.5,
        hjust = 0.5,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = c("TRUE" = "red", "FALSE" = "grey50")
      ) +
      ggplot2::labs(
        x = "Log2 Fold Change",
        y = "-Log10 P-Value",
        title = paste("Volcano Plot of Cytokine Levels:", c1, "vs", c2),
        colour = paste0(
          "Significant (p<",
          p_value_thresh,
          ", |FC|>",
          fold_change_thresh,
          ")"
        )
      ) +
      ggplot2::theme_minimal()
    plots[[paste(c1, "vs", c2)]] <- p
    # Optionally print final data frame
    if (verbose) {
      print(df[, setdiff(names(df), "label")], row.names = FALSE)
    }
  }
  return(plots)
}
