#' Extract Empirical Quantiles
#'
#' This function returns empirical quantiles from a numeric sample.
#' It uses a partial sort to avoid sorting the full vector when only
#' one order statistic is needed.
#'
#' @details
#' The function computes
#' 
#' q = `floor(perc * length(smpl))`
#' 
#' and extracts the q-th order statistic using sort with the partial
#' argument, which performs only the work needed to find the selected
#' element.
#' 
#' @param smpl Numeric vector containing the sample values.
#' 
#' @param perc Numeric value between 0 and 1 giving the desired
#' quantile level. For example, 0.5 returns the median.
#'
#' @return
#' A numeric value equal to the sample quantile at the specified
#' probability level.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' get_quantile(x, 0.5)
#' get_quantile(x, 0.9)
#' 
#' @keywords internal
#' 
get_quantile <- function(smpl, perc) {
	qnts <- floor(perc * length(smpl))
	return(sort(smpl, na.last = FALSE, partial = qnts)[qnts])
}
