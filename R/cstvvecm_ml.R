#' Estimation of I(1) CS-TV-VECMs Under Restrictions Using Maximum Likelihood
#'
#' Estimates an I(1) CS-TV-VECM under restrictions using maximum likelihood.
#' 
#' @param tsrs Matrix of values for the measured I(1) output.
#' 
#' @param rrst Specification of the restrictions on the long-run 
#' coefficients matrix. Either `NULL` or a list of length 3. If list of 
#' length 3, then first element restriction matrix, second element restriction
#' vector, third element vector of initial values for the free parameters.
#' 
#' @param lnrs Specification of the restrictions on the non-time-varying
#' adjustment coefficients matrix. Either `NULL` or a list. Must have
#' length 4 only if `ltrs` is `NULL` and `rrst` is `NULL`; 
#' otherwise must have length 3.
#' 
#'   - If list of length 4, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters, and fourth element initial value 
#'   for the long-run coefficients matrix.
#'   - If list of length 3, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters.
#' 
#' @param ltrs Specification of the restrictions on the time-varying
#' adjustment coefficients matrix. Either `NULL` or a list. Must have 
#' length 4 if `rrst` is `NULL`, and length 3 otherwise.
#' 
#'   - If list of length 4, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters, and fourth element initial value 
#'   for the long-run coefficients matrix.
#'   - If list of length 3, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters.
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
#' @param meth Method to be used in ``stats::optim''.
#' 
#' Default is `Nelder-Mead`.
#' 
#' @param lowr Lower bound for certain methods in ``stats::optim''.
#' 
#' Default is `-Inf'.
#' 
#' @param uppr Upper bound for certain methods in ``stats::optim''.
#' 
#' Default is `Inf'.
#' 
#' @param ctrl List of control parameters for ``stats::optim''.
#' 
#' Default is `list()'.
#' 
#' @param hess Logical. Should ``stats::optim'' return a numerically 
#' differentiated Hessian matrix?
#' 
#' Default is `FALSE'.
#' 
#' @return Returns a list of the following values:
#' 
#'   - `kalm`: Results obtained from Kalman smoother
#'   - `sime`: Estimated mean of initial state
#'   - `sicv`: Estimated covariance matrix of initial state
#'   - `evar`: Estimated covariance matrix of measurement equation errors
#'   - `sscf`: Estimated state transition matrix
#'   - `smea`: Estimated drift term in state equation
#'   - `secv`: Estimated covariance matrix of state equation errors
#'   - `augm`: Estimated coefficient matrix of input vectors in measurement 
#'      equation
#'   - `rfct`: Estimated long-run coefficients matrix
#'   - `lnfc`: Estimated non-time-varying part of adjustment coefficients matrix
#'   - `ltfc`: Estimated time-varying part of adjustment coefficient matrix
#'   - `rfur`: Estimated unrestricted long-run coefficients vector
#'   - `lnur`: Estimated unrestricted non-time-varying adjustment coefficients
#'     vector
#'   - `ltur`: Estimated unrestricted time-varying adjustment coefficients
#'     vector
#'   - `nlik`: Negative of the log-likelihood (apart from constants) evaluated 
#'     at the optimal solutions found
#'   - `ltrm`: Restriction matrix on time-varying adjustment coefficients matrix
#'   - `ltrv`: Restriction vector on time-varying adjustment coefficients matrix
#'   - `conv`: Integer code indicating convergece. `0` indicates successful
#'     completion.
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
#' 
cstvvecm_ml <- function(
    tsrs, rrst, lnrs, ltrs, evar, ordr = 1, dpow = 0,
    meth = "Nelder-Mead", lowr = -Inf, uppr = Inf, ctrl = list(), hess = FALSE
) {
  if (is.null(ltrs)) {
    if (is.null(lnrs)) {
      if (is.null(rrst)) {
        return(NULL)
      } else {
        if (nrow(rrst[[1]]) == length(rrst[[4]])) {
          excl <- FALSE
        } else {
          excl <- TRUE
        }
      }
    } else {
      if (is.null(rrst)) {
        if (nrow(lnrs[[1]]) == length(lnrs[[4]])) {
          excl <- FALSE
        } else {
          excl <- TRUE
        }
      } else {
        if (nrow(lnrs[[1]]) == nrow(rrst[[1]])) {
          excl <- FALSE
        } else {
          excl <- TRUE
        }
      } 
    }
  } else {
    if (is.null(rrst)) {
      if (nrow(ltrs[[1]]) == length(ltrs[[4]])) {
        excl <- FALSE
      } else {
        excl <- TRUE
      }
    } else {
      if (nrow(ltrs[[1]]) == nrow(rrst[[1]])) {
        excl <- FALSE
      } else {
        excl <- TRUE
      }
    }
  }
  rslt <- johi1setup(tsrs, ordr, dpow, excl, 0, 0, 0)
  rslt <- cstvrrr_ml(
    yvls = rslt$yvls,
    xvls = rslt$xvls,
    uvls = rslt$uvls,
    rrst = rrst,
    lnrs = lnrs,
    ltrs = ltrs,
    evar = evar,
    arst = NULL,
    meth = meth,
    lowr = lowr,
    uppr = uppr,
    ctrl = ctrl,
    hess = hess
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