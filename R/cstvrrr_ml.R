#' Estimation of CS-TV-RRRs Under Restrictions Using Maximum Likelihood
#'
#' Estimates an CS-TV-RRR under restrictions using maximum likelihood.
#' 
#' @param yvls Matrix of dependent variables.
#' 
#' @param xvls Matrix of (reduced rank coefficients related) regressors.
#' 
#' @param uvls Matrix of (not reduced rank coefficients related) regressors.
#' 
#' @param rrst Specification of the restrictions on the row space 
#' parameter matrix. Either `NULL` or a list of length 3. If list of 
#' length 3, then first element restriction matrix, second element restriction
#' vector, third element vector of initial values for the free parameters.
#' 
#' @param lnrs Specification of the restrictions on the non-time-varying part
#' of the column space parameter matrix. Either `NULL` or a list. Must have
#' length 4 only if `ltrs` is `NULL` and `rrst` is `NULL`; 
#' otherwise must have length 3.
#' 
#'   - If list of length 4, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters, and fourth element initial value 
#'   for the row space parameter matrix.
#'   - If list of length 3, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters.
#' 
#' @param ltrs Specification of the restrictions on the time-varying part
#' of the column space parameter matrix. Either `NULL` or a list. Must have 
#' length 4 if `rrst` is `NULL`, and length 3 otherwise.
#' 
#'   - If list of length 4, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters, and fourth element initial value 
#'   for the row space parameter matrix.
#'   - If list of length 3, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters.
#' 
#' @param evar Initial value for the covariance matrix of the error terms
#' in the measurement equation.
#' 
#' @param arst Specification of the restrictions on the coefficients of 
#' the input vectors in measurement equation. If not `NULL`, `NULL` is returned.
#' 
#' Default is `NULL`.
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
#'   - `rfct`: Estimated row space parameter matrix
#'   - `lnfc`: Estimated non-time-varying part of column space parameter matrix
#'   - `ltfc`: Estimated time-varying part of column space parameter matrix
#'   - `rfur`: Estimated unrestricted row space parameter vector
#'   - `lnur`: Estimated unrestricted non-time-varying column space parameter 
#'     vector
#'   - `ltur`: Estimated unrestricted time-varying column space parameter vector
#'   - `nlik`: Negative of the log-likelihood (apart from constants) evaluated 
#'     at the optimal solutions found
#'   - `ltrm`: Restriction matrix on time-varying column space parameter matrix
#'   - `ltrv`: Restriction vector on time-varying column space parameter matrix
#'   - `conv`: Integer code indicating convergece. `0` indicates successful
#'     completion.
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
cstvrrr_ml <- function(
    yvls, xvls, uvls, rrst, lnrs, ltrs, evar, arst = NULL, 
    meth = "Nelder-Mead", lowr = -Inf, uppr = Inf, ctrl = list(), hess = FALSE
) {
  if (!is.null(arst)) {
    return(NULL)
  }
  #
  lgth <- nrow(yvls)
  #
  secs <- sum(diag(stats::var(yvls)))/ncol(yvls)
  secs <- secs^2/10^2 # state error covariance matrix scale
  sics <- 10^3 # initial state covariance matrix scale
  #
  rfro <- ncol(xvls)
  lfro <- ncol(yvls)
  #
  if (is.null(ltrs)) {
    if (is.null(lnrs)) {
      if (is.null(rrst)) {
        return(NULL)
      } else {
        rfrm <- rrst[[1]]
        rfrv <- rrst[[2]]
        rfur <- rrst[[3]] # kappa
        #lnfc <- rrst[[4]]
        #
        lrco <- as.integer(length(rfrv)/rfro) # ncol(lnfc)
        #
        lnrm <- matrix(nrow = lfro * lrco, ncol = 0) # diag(1, lfro * lrco)
        lnrv <- numeric(lfro * lrco)
        lnur <- numeric(0) #c(lnfc) # gamma
        #
        ltrm <- diag(1, lfro * lrco)
        ltrv <- numeric(lfro * lrco)
        ltur <- numeric(lfro * lrco) # phi_t (initial)
      }
    } else {
      lnrm <- lnrs[[1]]
      lnrv <- lnrs[[2]]
      lnur <- lnrs[[3]] # gamma
      #
      lrco <- as.integer(length(lnrv)/lfro)
      #
      ltrm <- diag(1, lfro * lrco)
      ltrv <- numeric(lfro * lrco)
      ltur <- numeric(lfro * lrco) # phi_t (initial)
      if (is.null(rrst)) {
        rfct <- lnrs[[4]]
        #
        rfrm <- diag(1, rfro * lrco)
        rfrv <- numeric(rfro * lrco)
        rfur <- c(t(rfct)) # kappa
      } else {
        rfrm <- rrst[[1]]
        rfrv <- rrst[[2]]
        rfur <- rrst[[3]] # kappa
      }
    }
  } else {
    ltrm <- ltrs[[1]]
    ltrv <- ltrs[[2]]
    ltur <- ltrs[[3]] # phi_t (initial)
    #ltur <- numeric(lfro * lrco) 
    #
    lrco <- as.integer(length(ltrv)/lfro)
    if (is.null(lnrs)) {
      #lnfc <- ltrs[[4]][[1]]
      #
      lnrm <- matrix(nrow = lfro * lrco, ncol = 0) #diag(1, lfro * lrco)
      lnrv <- numeric(lfro * lrco)
      lnur <- numeric(0) #c(lnfc) # gamma
      #
      if (is.null(rrst)) {
        rfct <- ltrs[[4]]#[[2]]
        #
        rfrm <- diag(1, rfro * lrco)
        rfrv <- numeric(rfro * lrco)
        rfur <- c(t(rfct)) # kappa
      } else {
        rfrm <- rrst[[1]]
        rfrv <- rrst[[2]]
        rfur <- rrst[[3]] # kappa
      }
    } else {
      lnrm <- lnrs[[1]]
      lnrv <- lnrs[[2]]
      lnur <- lnrs[[3]] # gamma
      #
      if (is.null(rrst)) {
        rfct <- ltrs[[4]]#[[2]]
        #
        rfrm <- diag(1, rfro * lrco)
        rfrv <- numeric(rfro * lrco)
        rfur <- c(t(rfct)) # kappa
      } else {
        rfrm <- rrst[[1]]
        rfrv <- rrst[[2]]
        rfur <- rrst[[3]] # kappa
      }
    }
  }
  #
  args <- new.env()
  negloglike <- function(rfur, lnur, augm, evar, sscf, smea, secv, sime, sicv) { #sprd, cprd) {
    rfct <- rfrm %*% rfur + rfrv
    rfct <- matrix(rfct, rfro, lrco, TRUE) # beta update
    if (ncol(lnrm) > 0) {
      lnfc <- lnrm %*% lnur + lnrv
      lnfc <- matrix(lnfc, lfro, lrco) # alpha update
    }
    #
    xrfv <- xvls %*% rfct
    xarr <- array(t(xrfv), c(1, lrco, lgth))
    xidn <- xarr %x% diag(1, lfro)
    #
    mscf <- helperkit::arraymultmat(xidn, ltrm)
    mmea <- helperkit::arraymultvec(xidn, ltrv)
    mmea <- mmea + tcrossprod(xrfv, lnfc)
    mmea <- mmea + tcrossprod(uvls, augm)
    #
    kalm <- kalman_filter(yvls, mscf, mmea, evar, sscf, smea, secv, sime, sicv)
    sprd <- kalm$sprd
    cprd <- kalm$cprd
    #
    ltfc <- array(NA_real_, c(lfro, lrco, lgth))
    rslt <- 0
    for (time in 1:lgth) {
      msct <- ltrm %*% sprd[time, ] + ltrv
      tmpa <- (xvls[time, , drop = FALSE] %*% rfct) %x% diag(1, lfro)
      tmpb <- tmpa %*% msct
      ltfc[, , time] <- matrix(msct, lfro, lrco)
      if (ncol(lnrm) > 0) {
        tmpb <- tmpb + lnfc %*% crossprod(rfct, xvls[time, ])
      }
      tmpb <- yvls[time, ] - tmpb - augm %*% uvls[time, ]
      tmpc <- tmpa %*% ltrm
      tmpd <- helperkit::array3tomat(cprd, time)
      tmpc <- tcrossprod(tmpc %*% tmpd, tmpc) + evar
      #tmpb <- crossprod(tmpb, helperkit::safesolve(tmpc) %*% tmpb)
      tmpb <- sum(tmpb * (helperkit::safesolve(tmpc) %*% tmpb))
      rslt <- rslt + helperkit::ldet(tmpc) + tmpb
    }
    args$last <- list(kalm = kalm, rfct = rfct, lnfc = lnfc, ltfc = ltfc)
    return(drop(rslt)/lgth)
  }
  ofun <- function(rfur, lnur, augm, evar, sscf, smea, secv, sime, sicv) {
    npar <- list(
      rfur = rfur, #kappa
      lnur = lnur, #gamma_nu
      augm = augm, #Famma
      evar = evar, #Omega
      sscf = sscf, #Pi
      smea = smea, #pi
      secv = secv, #Sigma
      sime = sime, #mu
      sicv = sicv #Psi
    )
    tmpl <- pack::make_template(npar)
    opar <- pack::pack(npar)
    objf <- function(opar) {
      npar <- pack::unpack(opar, tmpl)
      return(negloglike(#(rfur, lnur, augm, evar, sscf, smea, secv, sime, sicv
        rfur = npar[[1]],
        lnur = npar[[2]],
        augm = npar[[3]],
        evar = npar[[4]],
        sscf = npar[[5]],
        smea = npar[[6]],
        secv = npar[[7]],
        sime = npar[[8]],
        sicv = npar[[9]]
      ))
    }
    opti <- stats::optim(
      par = opar, 
      fn = objf, 
      gr = NULL, 
      method = meth,
      lower = lowr,
      upper = uppr,
      control = ctrl,
      hessian = hess
    )
    para <- pack::unpack(opti$par, tmpl)
    return(list(
      kalm = args$last$kalm,
      sime = para[[8]],
      sicv = para[[9]],
      evar = para[[4]],
      sscf = para[[5]],
      smea = para[[6]],
      secv = para[[7]],
      augm = para[[3]],
      rfct = args$last$rfct,
      lnfc = args$last$lnfc,
      ltfc = args$last$ltfc,
      rfur = para[[1]],
      lnur = para[[2]],
      ltur = args$last$kalm$pred,
      nlik = opti$fn,
      ltrm = ltrm,
      ltrv = ltrv,
      conv = opti$convergence
    ))
  }
  #
  rfct <- rfrm %*% rfur + rfrv
  rfct <- matrix(rfct, rfro, lrco, TRUE) # beta
  #
  lnfc <- lnrm %*% lnur + lnrv
  lnfc <- matrix(lnfc, lfro, lrco) # nu
  #
  ltfc <- ltrm %*% ltur + ltrv
  ltfc <- matrix(ltfc, lfro, lrco) # alpha (initial)
  #
  secv <- diag(secs, length(ltur))
  sicv <- sics * secv # initial state covariance matrix
  sime <- ltur #numeric(length(ltur)) # initial state mean
  #
  smea <- numeric(length(ltur)) # pi (initial)
  sscf <- diag(1, length(ltur)) # Pi (initial)
  #
  if (is.null(uvls)) {
    uvls <- matrix(nrow = lgth, ncol = 0)
    uinv <- matrix(nrow = 0, ncol = 0)
  } else {
    uinv <- helperkit::safesolve(crossprod(uvls))
  }
  xrfv <- xvls %*% rfct
  augm <- tcrossprod(xrfv, ltfc + lnfc)
  augm <- crossprod(yvls - augm, uvls) %*% uinv
  #
  return(ofun(rfur, lnur, augm, evar, sscf, smea, secv, sime, sicv))
}
