#' Johansen (1995) Test
#'
#' Performs a Johansen (1995) test
#' 
#' @param sigl Significance level.
#' 
#' @param tsrs (Multivariate) time series to be tested for the cointegrating
#' rank.
#' 
#' @param ordr Order of the VECM. Equal to one plus the number of lags.
#' 
#' Default is `1`.
#' 
#' @param dpow Highest power of polynomial time trends (e.g., -1: 
#' no deterministic components, 0: constant only, 1: constant and linear trend)
#' 
#' Default is `0`.
#' 
#' @param excl Logical value indicating whether the highest power 
#' of polynomial time trends is restricted to the long-run 
#' relationship.
#' 
#' Default is `FALSE`.
#' 
#' @param smpl Sample size to consider in case critical values 
#' have to be simulated.
#' 
#' @param simu Number of repetitions in case critical values have to 
#' be simulated.
#' 
#' @param tolr Tolerance limit for lookup tables of critical values.
#'
#' @return A list containing:
#' 
#'   - `stat`: Test statistics
#'   - `crit`: Critical values
#'   - `rjct`: Logical value indicating the rejection of null hypotheses
#' 
#' @references
#' Johansen, S. (1995). Likelihood-Based Inference in Cointegrated Vector
#' Autoregressive Models. Oxford University Press, Oxford.
#' 
#' @export
#' 
johi1test <- function(
	sigl, tsrs, ordr = 1, dpow = 0, excl = FALSE,
    smpl = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    excl <- ifelse(dpow == -1, FALSE, excl)
    stat <- johi1stat(tsrs, ordr, dpow, excl)
    crit <- getjoh95qntl(1 - sigl, dpow, NCOL(tsrs), excl, smpl, simu, tolr)
    if (is.matrix(crit)) {
        rjct <- t(apply(crit, 1, function(x) {return(stat > x)}))
    } else {
        rjct <- (stat > crit)
    }
    return(list(stat = stat, crit = crit, rjct = rjct))
}
