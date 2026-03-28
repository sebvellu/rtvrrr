#' Estimation of I(1) RS-TV-VECMs Under Restrictions
#'
#' Estimates an I(1) RS-TV-VECM under restrictions.
#' 
#' @param tsrs Matrix of values for the measured I(1) output.
#' 
#' @param rtrs Specification of the restrictions on the time-varying
#' long-run coefficients matrix.  Either `NULL` or a list. Must have 
#' length 4 if `lrst` is `NULL`, and length 3 otherwise.
#' 
#'   - If list of length 4, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters, and fourth element initial value 
#'   for the adjustment coefficients matrix.
#'   - If list of length 3, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters.
#' 
#' @param rnrs Specification of the restrictions on the non-time-varying
#' long-run coefficients matrix. Either `NULL` or a list. Must have 
#' length 4 only if `rtrs` is `NULL` and `lrst` is `NULL`; 
#' otherwise must have length 3.
#' 
#'   - If list of length 4, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters, and fourth element initial value 
#'   for the adjustment coefficients matrix.
#'   - If list of length 3, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters.
#' 
#' @param lrst Specification of the restrictions on the adjustment coefficients
#' matrix. Either `NULL` or a list of length 3. If list of 
#' length 3, then first element restriction matrix, second element restriction
#' vector, third element vector of initial values for the free parameters.
#' 
#' @param evar Initial value for the covariance matrix of the error terms
#' in the measurement equation.
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
#' @param bsbv Logical indication whether the state equation considered 
#' by Brune, Scherrer and Bura (2022) should be used.
#' 
#' Default is `FALSE`.
#' 
#' @param tolr Tolerance limit for optimization. 
#' 
#' Default is 100 * sqrt(.Machine$double.eps). 
#' 
#' @param mxit Maximum number of iterations (outer loop).
#' 
#' Default is 10000.
#' 
#' @param mxis Maximum number of iterations (inner loop, switching).
#' 
#' Default is 1.
#' 
#' @return Returns `NULL` if `arst` is not `NULL`, otherwise 
#' a list of the following values:
#' 
#'   - `kalm`: Results obtained from Kalman smoother
#'   - `sime`: Estimated mean of initial state
#'   - `sicv`: Estimated covariance matrix of initial state
#'   - `evar`: Estimated covariance matrix of measurement equation errors
#'   - `sscf`: Estimated state transition matrix
#'   - `smea`: Estimated drift term in state equation
#'   - `secv`: Estimated covariance matrix of state equation errors
#'   - `augm`: Estimated coefficient matrix of input vectors in measurement 
#'     equation
#'   - `lfct`: Estimated adjumstment coefficients matrix
#'   - `rnfc`: Estimated non-time-varying part of long-run coefficients matrix
#'   - `rtfc`: Estimated time-varying part of long-run coefficients matrix
#'   - `lfur`: Estimated unrestricted adjustment coefficients vector
#'   - `rnur`: Estimated unrestricted non-time-varying long-run coefficients 
#'      vector
#'   - `rtur`: Estimated unrestricted time-varying long-run coefficients vector
#'   - `nlik`: Negative of the log-likelihood (apart from constants) evaluated 
#'     at the optimal solutions found
#'   - `rtrm`: Restriction matrix on time-varying long-run coefficients matrix
#'   - `rtrv`: Restriction vector on time-varying long-run coefficients matrix
#'   - `ordr`: Order of the VECM
#'   - `dpow`: Highest power of polynomial time trends
#'   - `excl`: Logical indicating whether the highest power of polynomial
#'     time trends should only be captured in the long-run relationship
#' 
#' @references
#' Brune, B., Scherrer, W. and Bura, E. (2022). A State-Space Approach to
#' Time-Varying Reduced-Rank Regression. Econometric Reviews 41, 895-917.
#' 
#' Veldhuis, S. and Wagner, M. (2026). Restrictions in a State Space Approach 
#' to Time-Varying Vector Error Correction Models: Modelling Instabilities 
#' in Long-Run Money Demand. Mimeo.
#' 
#' @export

rstvvecm <- function(
	tsrs, rtrs, rnrs, lrst, evar, ordr = 1, dpow = 0, bsbv = FALSE,
    tolr = 100 * sqrt(.Machine$double.eps), mxit = 10000, mxis = 1
) {
	if (is.null(rtrs)) {
        if (is.null(rnrs)) {
            if (is.null(lrst)) {
                return(NULL)
            } else {
                if (nrow(lrst[[1]]) == length(lrst[[4]])) {
                    excl <- FALSE
                } else {
                    excl <- TRUE
                }
            }
        } else {
            if (is.null(lrst)) {
                if (nrow(rnrs[[1]]) == length(rnrs[[4]])) {
                    excl <- FALSE
                } else {
                    excl <- TRUE
                }
            } else {
                if (nrow(rnrs[[1]]) == nrow(lrst[[1]])) {
                    excl <- FALSE
                } else {
                    excl <- TRUE
                }
            } 
        }
    } else {
        if (is.null(lrst)) {
            if (nrow(rtrs[[1]]) == length(rtrs[[4]])) {
                excl <- FALSE
            } else {
                excl <- TRUE
            }
        } else {
            if (nrow(rtrs[[1]]) == nrow(lrst[[1]])) {
                excl <- FALSE
            } else {
                excl <- TRUE
            }
        }
    }
	rslt <- johi1setup(tsrs, ordr, dpow, excl, 0, 0, 0)
	rslt <- rstvrrr(
		yvls = rslt$yvls,
		xvls = rslt$xvls,
		uvls = rslt$uvls,
		rtrs = rtrs,
		rnrs = rnrs,
		lrst = lrst,
		evar = evar,
		arst = NULL,
		bsbv = bsbv,
        tolr = tolr,
        mxit = mxit,
        mxis = mxis
	)
    rslt <- c(
        rslt,
        list(
			ordr = ordr,
			dpow = dpow,
			excl = excl
		)
    )
	#
	return(rslt)
}