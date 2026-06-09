#' URT upper respiratory tract microbiome data (filtered)
#'
#' A small example dataset derived from Charlson et al. (2010) as distributed
#' by the **LOCOM** package and filtered to left oropharyngeal samples.
#' Provided for illustrating CAFT workflows.
#'
#' @format A list with three components:
#' \describe{
#'   \item{otu}{Integer matrix of counts, samples in rows and taxa in columns.}
#'   \item{meta}{Data frame of sample metadata (e.g., SmokingStatus, Sex).}
#'   \item{tax}{Data frame of taxonomy annotations for the taxa.}
#' }
#' @details Zeros in \code{otu} correspond to under-detection and are handled
#' as left-censored measurements by CAFT.
#'
#' @source Charlson ES et al. (2010) PLoS One 5(12):e15216; dataset accessed via
#' the **LOCOM** R package.
#'
#' @examples
#' data("URT", package = "CAFT")
#' str(URT, max.level = 1)
#' # counts:
#' dim(URT$otu)
#' # metadata variables:
#' names(URT$meta)
#'
#' @usage data("URT")
#' @keywords datasets
"URT"
