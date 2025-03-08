##############################################################################
## Function to draw an error-bar plot enriched with p-value                  #
## and/or effect size for the comparison of each of selected                 #
## groups to the baseline group                                              #
## Author: Xiaohua Douglas Zhang, January 2022                               #
## Arguments:                                                                #
##   data: a data frame containing the following columns for each group:     #
##            "name" for group names                                         #
##            "center" for mean or median,                                   #
##            "spread" for standard deviation, MAD or standard error,        #
##            "p.value" for p-value,                                         #
##            "effect.size" for effect size based on SSMD                    #
##            note: the first column of data must be for the baseline group  #
##   pLab: whether to label the p-value                                      #
##   esLab: whether to label the effect size usually in SSMD                 #
## Output:                                                                   #
##   None                                                                    #
##############################################################################

#' Error-bar Plot.
#'
#' This function draws an error-bar plot for comparing groups to a
#' baseline group. It creates a barplot of the central tendency
#' (mean or median) and overlays error bars representing the
#' spread (e.g., standard deviation, MAD, or standard error).
#' Optionally, p-value and effect size labels (based on SSMD) are added,
#' either as symbols or numeric values.
#'
#' @param data A data frame containing the following columns for each group:
#'   \itemize{
#'     \item \code{name}: Group names.
#'     \item \code{center}: Mean or median values.
#'     \item \code{spread}: Standard deviation, MAD, or standard error.
#'     \item \code{p.value}: P-value for the comparison.
#'     \item \code{effect.size}: Effect size based on SSMD.
#'   }
#'   Note: The first row of \code{data} must correspond to the baseline group.
#' @param p_lab Logical. Whether to label the p-values on the plot.
#' Default is \code{TRUE}.
#' @param es_lab Logical. Whether to label the effect sizes on the plot.
#' Default is \code{TRUE}.
#' @param class_symbol Logical. Whether to use symbolic notation for
#' significance and effect size.
#'   Default is \code{TRUE}.
#' @param x_lab Character. Label for the x-axis.
#' @param y_lab Character. Label for the y-axis.
#' @param main Character. Title of the graph.
#'
#' @return An error-bar plot is produced.
#'
#' @export
#' @examples
#' # Load sample data
#' data_df <- cytodata[, -1]
#' cyt_mat <- log2(data_df[, -c(1:3)])
#' data_df1 <- data.frame(data_df[, 1:3], cyt_mat)
#' cytokineNames <- colnames(cyt_mat)
#' nCytokine <- length(cytokineNames)
#' condt <- !is.na(cyt_mat) & (cyt_mat > 0)
#' Cutoff <- min(cyt_mat[condt], na.rm = TRUE) / 10
#' # Create matrices for ANOVA and Tukey results
#' p_aov_mat <- matrix(NA, nrow = nCytokine, ncol = 3)
#' dimnames(p_aov_mat) <- list(cytokineNames,
#'                          c("Group", "Treatment", "Interaction"))
#' p_groupComp_mat <- matrix(NA, nrow = nCytokine, ncol = 3)
#' dimnames(p_groupComp_mat) <- list(cytokineNames,
#'                                  c("2-1", "3-1", "3-2"))
#' ssmd_groupComp_stm_mat <- mD_groupComp_stm_mat <- p_groupComp_stm_mat <-
#' p_groupComp_mat
#'
#' for (i in 1:nCytokine) {
#' Cytokine <- (cyt_mat[, i] + Cutoff)
#' cytokine_aov <- aov(Cytokine ~ Group * Treatment, data = data_df)
#' aov_table <- summary(cytokine_aov)[[1]]
#' p_aov_mat[i, ] <- aov_table[1:3, 5]
#' p_groupComp_mat[i, ] <- TukeyHSD(cytokine_aov)$Group[1:3, 4]
#' p_groupComp_stm_mat[i, ] <- TukeyHSD(cytokine_aov)$`Group:Treatment`[1:3, 4]
#' mD_groupComp_stm_mat[i, ] <- TukeyHSD(cytokine_aov)$`Group:Treatment`[1:3, 1]
#' ssmd_groupComp_stm_mat[i, ] <- mD_groupComp_stm_mat[i, ] / sqrt(2 *
#' aov_table["Residuals", "Mean Sq"])
#' }
#'
#' results <- cyt_skku(cytodata[, -c(1,4)], print_res_log = TRUE,
#'                  group_cols = c("Group", "Treatment"))
#' pdf("bar_error_plot_enriched.pdf")
#' par(mfrow = c(2,3), mar = c(8.1, 4.1, 4.1, 2.1))
#' for (k in 1:nCytokine) {
#' result_mat <- results[1:9, , k]
#'center_df <- data.frame(
#'  name = rownames(result_mat),
#'  result_mat[, c("center", "spread")],
#'  p.value = c(1, p_groupComp_stm_mat[k, 1:2]),
#'  effect.size = c(0, ssmd_groupComp_stm_mat[k, 1:2])
#'  )
#' cyt_errbp(center_df, p_lab = TRUE, es_lab = TRUE,
#'          class_symbol = TRUE,
#'          y_lab = "Concentration in log2 scale",
#'          main = cytokineNames[k])
#' }
#' dev.off()
cyt_errbp <- function(data, p_lab = TRUE, es_lab = TRUE, class_symbol = TRUE,
                      x_lab = "", y_lab = "", main = "") {
  # Significance mark function
  significance_mark_fn <- function(p_value) {
    ###########################################################################
    ## function to mark the significant level based on p-value
    ## Input:
    ##   p_value: a vector of p-values
    ## Author: Xiaohua Douglas Zhang, 2022
    ## Reference:
    ##   Zhang XHD, 2011. Optimal High-Throughput Screening: Practical
    ##   Experimental Design and Data Analysis for Genome-scale RNAi Research.
    ##   Cambridge University Press, Cambridge, UK
    ###########################################################################
    significance_class_vec <- rep(NA, length(p_value))
    significance_class_vec[p_value <= 0.05 & p_value > 0.01] <- "*"
    significance_class_vec[p_value <= 0.01 & p_value > 0.001] <- "**"
    significance_class_vec[p_value <= 0.001 & p_value > 0.0001] <- "***"
    significance_class_vec[p_value <= 0.0001 & p_value > 0.00001] <- "****"
    significance_class_vec[p_value <= 0.00001] <- "*****"
    if (sum(is.na(p_value)) > 0)
      significance_class_vec[is.na(p_value)] <- NA
    significance_class_vec
  }

  # Effect size mark function
  effect_size_mark_fn <- function(ssmd) {
    ##########################################################################
    ## function to mark SSMD-based effect size being large or beyond
    ## Input:
    ##   ssmd: a vector of SSMD values
    ## Author: Xiaohua Douglas Zhang, 2022
    ## Reference:
    ##   Zhang XHD, 2011. Optimal High-Throughput Screening: Practical
    ##   Experimental Design and Data Analysis for Genome-scale RNAi Research.
    ##   Cambridge University Press, Cambridge, UK
    ##########################################################################
    effect_class_vec <- rep(NA, length(ssmd))
    effect_class_vec[ssmd >= 5] <- ">>>>>"
    effect_class_vec[ssmd >= 3 & ssmd < 5] <- ">>>>"
    effect_class_vec[ssmd >= 1.645 & ssmd < 3] <- ">>>"
    effect_class_vec[ssmd >= 1 & ssmd < 1.645] <- ">>"
    effect_class_vec[ssmd > 0.25 & ssmd < 1] <- ">"
    effect_class_vec[ssmd < 0.25 & ssmd > -0.25] <- " "
    effect_class_vec[ssmd <= -0.25 & ssmd > -1] <- "<"
    effect_class_vec[ssmd <= -1 & ssmd > -1.645] <- "<<"
    effect_class_vec[ssmd <= -1.645 & ssmd > -3] <- "<<<"
    effect_class_vec[ssmd <= -3 & ssmd > -5] <- "<<<<"
    effect_class_vec[ssmd <= -5] <- "<<<<<"
    if (sum(is.na(ssmd)) > 0)
      effect_class_vec[is.na(ssmd)] <- NA
    effect_class_vec
  }

  center <- data[, "center"]
  names(center) <- data[, "name"]
  spread <- data[, "spread"]
  mid <- barplot(center, plot = FALSE)
  p_value <- data[-1, "p.value"]
  effect_size <- data[-1, "effect.size"]

  # Plot barplot
  y_lim0 <- c(
    min(0, min(center - spread)),
    max(0, max(center + spread))
  ) * 1.4
  y_range <- diff(y_lim0)
  y_lim <- range(c(y_lim0, center[-1] +
                     sign(center[-1]) * (spread[-1] + y_range / 8)))
  barplot(center, ylim = y_lim, las = 2, xlab = x_lab,
          ylab = y_lab, main = main)

  # Add error bars
  arrows(x0 = mid, y0 = center, x1 = mid, y1 = center +
           sign(center) * spread, code = 3, angle = 90, length = 0.1)

  # Add text for p-value
  if (p_lab) {
    if (class_symbol) {
      text(
        x = mid[-1],
        y = center[-1] + sign(center[-1]) * (spread[-1] + y_range / 20),
        labels = significance_mark_fn(p_value), cex = 1
      )
    } else {
      text(
        x = mid[-1],
        y = center[-1] + sign(center[-1]) * (spread[-1] + y_range / 20),
        labels = paste("p=", ifelse(p_value > 0.001, round(p_value, 3),
                                    formatC(p_value, format = "e",
                                            digits = 1)), sep = ""),
        cex = 0.75
      )
    }
  }

  # Add text for effect size
  if (es_lab) {
    if (class_symbol) {
      text(
        x = mid[-1],
        y = center[-1] + sign(center[-1]) * (spread[-1] + y_range / 8),
        labels = effect_size_mark_fn(effect_size), cex = 1
      )
    } else {
      text(
        x = mid[-1],
        y = center[-1] + sign(center[-1]) * (spread[-1] + y_range / 8),
        labels = paste("s=", round(effect_size, 3), sep = ""),
        cex = 0.75
      )
    }
  }
}
