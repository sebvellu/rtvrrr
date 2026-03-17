#' Estimation of I(1) VECMs Under Restrictions
#'
#' Estimates an I(1) VECM under restrictions.
#' 
#' @param tsrs Matrix of values for multivariate I(1) process.
#' 
#' @param rrst Specification of the restrictions on the long-run 
#' coefficients matrix. Either `NULL` or a list of length 3. If list of 
#' length 3, then first element restriction matrix, second element restriction
#' vector, third element vector of initial values for the free parameters.
#' 
#' @param lrst Specification of the restrictions on the adjustment
#' coefficients matrix. Either `NULL` or a list. Must have length 4 if 
#' `rrst` is `NULL`, and length 3 otherwise.
#' 
#'   - If list of length 4, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters, and fourth element initial value 
#'   for the long-run coefficients matrix.
#'   - If list of length 3, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters.
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
#' @param fish Should fisher information matrix be returned?
#' 
#' Default is `FALSE`.
#' 
#' @param bima Order of the Chebyshev time polynomials for a time-varying
#' VECM as described in Bierens and Martins (2010). Set to 0 for the
#' standard time-invariant VECM case.
#'
#' Default is `0`.
#' 
#' @param tolr Tolerance limit for optimization. 
#' 
#' Default is 100 * sqrt(.Machine$double.eps). 
#' 
#' @param mxit Maximum number of iterations.
#' 
#' Default is 10000.
#' 
#' @return Returns a list of the following values:
#' 
#'   - `lfct`: Estimated adjument coefficients matrix
#'   - `rfct`: Estimated long-run coefficients matrix
#'   - `lfun`: Estimated unrestricted adjustment coefficients vector
#'   - `rfun`: Estimated unrestricted long-run coefficients vector
#'   - `evar`: Estimated variance-covariance matrix of the error term
#'   - `augm`: Estimated coefficient matrix of the short-run dynamics
#'   - `lrco`: Cointegrating rank
#'   - `fish`: Fisher information marix (if argument `fish` is set to
#'     `TRUE`, otherwise `NULL`)
#'   - `unvr`: Inverse of the Fisher information matrix (if argument
#'     `fish` is set to `TRUE`, otherwise `NULL`)
#'   - `idnt`: Logical indicating whether identification is achieved 
#'     (if argument `fish` is set to `TRUE`, otherwise `NULL`)
#'   - `nlik`: Negative of the log-likelihood (apart from constants) evaluated 
#'     at the optimal solutions found
#'   - `ordr`: Order of the VECM
#'   - `dpow`: Highest power of polynomial time trends
#'   - `excl`: Logical indicating whether the highest power of polynomial
#'     time trends should only be captured in the long-run relationship
#'   - `bima`: Order of the Chebyshev time polynomials for a time-varying
#'     VECM as described in Bierens and Martins (2010)
#' 
#' @references
#' Bierens, H. J. and Martins, L. F. (2010). Time-Varying Cointegration.
#' Econometric Theory 26, 1453-1490.
#' 
#' Boswijk, H. P. and Doornik, J. A. (2004). Identifying, Estimating and
#' Testing Restricted Cointegrated Systems: An Overview. Statistica
#' Neerlandica 58, 440-465.
#' 
#' Doornik, J. A. (1995). Testing General Restrictions on the Cointegrating
#' Space. Working Paper, Nuffield College, Oxford.
#' 
#' Hansen, P. R. (2003). Generalized Reduced Rank Regression. Working Paper,
#' Brown University.
#' 
#' Johansen, S. (1995). Likelihood-Based Inference in Cointegrated Vector
#' Autoregressive Models. Oxford University Press, Oxford.
#' 
#' Pesaran, M. H. and Shin, Y. (2002). Long-Run Structural Modelling.
#' Econometric Reviews 21, 49-87.
#' 
#' @export
#' 
vecm <- function(
	tsrs, rrst, lrst, ordr = 1, dpow = 0, fish = FALSE, bima = 0,
	tolr = 100 * sqrt(.Machine$double.eps), mxit = 10000
) {
	nums <- ncol(tsrs)
	if (is.null(lrst)) {
		if (is.null(rrst)) {
			return(NULL)
		} else {
			if (nrow(rrst[[1]]) == (length(rrst[[4]]) * (bima + 1))) {
				excl <- FALSE
			} else {
				excl <- TRUE
			}
		}
	} else {
		if (is.null(rrst)) {
			if (nrow(lrst[[1]]) == length(lrst[[4]])) {
				excl <- FALSE
			} else {
				excl <- TRUE
			}
		} else {
			if (nrow(lrst[[1]]) == nrow(rrst[[1]])) {
				excl <- FALSE
			} else {
				excl <- TRUE
			}
		} 
	}
	rslt <- johi1setup(tsrs, ordr, dpow, excl, 0, 0, bima) #dtrn, bima)
	#
	rslt <- rrr(
		yvls = rslt$yvls,
		xvls = rslt$xvls,
		uvls = rslt$uvls,
		rrst = rrst,
		lrst = lrst,
		arst = NULL,
		fish = fish,
		tolr = tolr,
		mxit = mxit
	)
	#
	return(c(
		rslt,
		list(
			ordr = ordr,
			dpow = dpow,
			excl = excl,
			bima = bima
		)
	))
}