#' Sample Cytokine Profiling Data.
#'
#' Contains observed values of cytokines and their respective treatment
#' and groups, derived from:
#'
#' @format A data frame with 297 rows and 29 columns:
#' \describe{
#'   \item{Group}{Group assigned to the subjects.}
#'   \item{Stimulation}{Stimulation received by subjects.}
#'   \item{Treatment}{Treatment received by subjects.}
#'   \item{...}{Additional numeric columns representing measured cytokines.}
#' }
#' @source Example data compiled for cytokine profiling.
#' @references
#' Pugh GH, Fouladvand S, SantaCruz-Calvo S, Agrawal M, Zhang XD, Chen J,
#' Kern PA, Nikolajczyk BS.
#' T cells dominate peripheral inflammation in a cross-sectional analysis of
#' obesity-associated diabetes. \emph{Obesity (Silver Spring)}. 2022;30(10):
#' 1983–1994. doi:10.1002/oby.23528.
#'
#' @examples
#' data(cytodata)
"cytodata"
