#' Estimation of RS-TV-RRRs Under Restrictions
#'
#' Estimates an RS-TV-RRR under restrictions.
#' 
#' @param yvls Matrix of dependent variables.
#' 
#' @param xvls Matrix of (reduced rank coefficients related) regressors.
#' 
#' @param uvls Matrix of (not reduced rank coefficients related) regressors.
#' 
#' @param rtrs Specification of the restrictions on the time-varying
#' row space parameter matrix. Either `NULL` or a list. Must have 
#' length 4 if `lrst` is `NULL`, and length 3 otherwise.
#' 
#'   - If list of length 4, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters, and fourth element initial value 
#'   for the column space parameter matrix.
#'   - If list of length 3, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters.
#' 
#' @param rnrs Specification of the restrictions on the non-time-varying
#' row space parameter matrix. Either `NULL` or a list. Must have 
#' length 4 only if `rtrs` is `NULL` and `lrst` is `NULL`; 
#' otherwise must have length 3.
#' 
#'   - If list of length 4, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters, and fourth element initial value 
#'   for the column space parameter matrix.
#'   - If list of length 3, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters.
#' 
#' @param lrst Specification of the restrictions on the column space 
#' parameter matrix. Either `NULL` or a list of length 3. If list of 
#' length 3, then first element restriction matrix, second element restriction
#' vector, third element vector of initial values for the free parameters.
#' 
#' @param evar Initial value for the covariance matrix of the error terms
#' in the measurement equation.
#' 
#' @param arst Specification of the restrictions on the coefficients of 
#' the input vectors in measurement equation. If not `NULL`, `NULL` is returned.
#' 
#' Default is `NULL`.
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
#' Default is Inf.
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
#'   - `lfct`: Estimated column space parameter matrix
#'   - `rnfc`: Estimated non-time-varying part of row space parameter matrix
#'   - `rtfc`: Estimated time-varying part of row space parameter matrix
#'   - `lfur`: Estimated unrestricted column space parameter vector
#'   - `rnur`: Estimated unrestricted non-time-varying row space parameter 
#'     vector
#'   - `rtur`: Estimated unrestricted time-varying row space parameter vector
#'   - `nlik`: Negative of the log-likelihood (apart from constants) evaluated 
#'     at the optimal solutions found
#'   - `rtrm`: Restriction matrix on time-varying row space parameter matrix
#'   - `rtrv`: Restriction vector on time-varying row space parameter matrix
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

rstvrrr <- function(
	yvls, xvls, uvls, rtrs, rnrs, lrst, evar, arst = NULL, bsbv = FALSE,
	tolr = 100 * sqrt(.Machine$double.eps), mxit = Inf, mxis = 1
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
	negloglike <- function(mscf, mmea, mecv, sprd, cprd) {
		rslt <- 0
		for (time in 1:lgth) {
			msct <- helperkit::array3tomat(mscf, time)
			tmpa <- yvls[time, ] - msct %*% sprd[time, ] - mmea[time, ]
			tmpb <- helperkit::array3tomat(cprd, time)
			tmpb <- tcrossprod(msct %*% tmpb, msct) + mecv
			tmpa <- crossprod(tmpa, helperkit::safesolve(tmpb) %*% tmpa)
			rslt <- rslt + helperkit::ldet(tmpb) + tmpa
		}
		return(drop(rslt)/lgth)
	}
	stopcrit <- function(fnew, fold) {
		return(abs(fnew - fold) < tolr * max(1, abs(fold)))
	}
	#
	rfro <- ncol(xvls)
    lfro <- ncol(yvls)
    #
    if (is.null(rtrs)) {
        if (is.null(rnrs)) {
			if (is.null(lrst)) {
                return(NULL)
            } else {
                lfrm <- lrst[[1]]
                lfrv <- lrst[[2]]
                lfur <- lrst[[3]] # gamma
                #rnfc <- lrst[[4]]
                #
                lrco <- as.integer(length(lfrv)/lfro) # ncol(rnfc)
                #
                rnrm <- matrix(nrow = rfro * lrco, ncol = 0) # diag(1, rfro * lrco)
                rnrv <- numeric(rfro * lrco)
                rnur <- numeric(0) #c(rnfc) # kappa
                #
                rtrm <- diag(1, rfro * lrco)
                rtrv <- numeric(rfro * lrco)
                rtur <- numeric(rfro * lrco) # phi_t (initial)
            }
        } else {
            rnrm <- rnrs[[1]]
            rnrv <- rnrs[[2]]
            rnur <- rnrs[[3]] # kappa
            #
            lrco <- as.integer(length(rnrv)/lfro)
            #
            rtrm <- diag(1, rfro * lrco)
            rtrv <- numeric(rfro * lrco)
            rtur <- numeric(rfro * lrco) # phi_t (initial)
			if (is.null(lrst)) {
                lfct <- rnrs[[4]]
                #
                lfrm <- diag(1, lfro * lrco)
                lfrv <- numeric(lfro * lrco)
                lfur <- c(t(lfct))  # gamma
            } else {
                lfrm <- lrst[[1]]
                lfrv <- lrst[[2]]
                lfur <- lrst[[3]] # kappa
            }
        }
    } else {
        rtrm <- rtrs[[1]]
        rtrv <- rtrs[[2]]
        rtur <- rtrs[[3]] # phi_t (initial)
        #ltur <- numeric(rfro * lrco) 
        #
        lrco <- as.integer(length(rtrv)/rfro)
        if (is.null(rnrs)) {
            #rnfc <- rtrs[[4]][[1]]
            #
            rnrm <- matrix(nrow = rfro * lrco, ncol = 0) #diag(1, rfro * lrco)
            rnrv <- numeric(rfro * lrco)
            rnur <- numeric(0) #c(rnfc) # kappa
            #
			if (is.null(lrst)) {
                lfct <- rtrs[[4]]#[[2]]
                #
                lfrm <- diag(1, lfro * lrco)
                lfrv <- numeric(lfro * lrco)
                lfur <- c(t(lfct)) # gamma
            } else {
                lfrm <- lrst[[1]]
                lfrv <- lrst[[2]]
                lfur <- lrst[[3]] # gamma
            }
        } else {
            rnrm <- rnrs[[1]]
            rnrv <- rnrs[[2]]
            rnur <- rnrs[[3]] # kappa
            #
			if (is.null(lrst)) {
                lfct <- rtrs[[4]]#[[2]]
                #
                lfrm <- diag(1, lfro * lrco)
                lfrv <- numeric(lfro * lrco)
                lfur <- c(t(lfct)) # gamma
            } else {
                lfrm <- lrst[[1]]
                lfrv <- lrst[[2]]
                lfur <- lrst[[3]] # gamma
            }
        }
    }
    #
    rnfc <- rnrm %*% rnur + rnrv
    rnfc <- matrix(rnfc, rfro, lrco, TRUE) # varrho
	#
    rtfc <- rtrm %*% rtur + rtrv
    rtfc <- matrix(rtfc, rfro, lrco, TRUE) # beta (initial)
    #
    lfct <- lfrm %*% lfur + lfrv
    lfct <- matrix(lfct, lfro, lrco) # alpha
    #
    secv <- diag(secs, length(rtur))
    sicv <- sics * secv # initial state covariance matrix
    sime <- rtur #numeric(length(rtur)) # initial state mean
    #
    smea <- numeric(length(rtur)) # pi (initial)
    sscf <- diag(1, length(rtur)) # Pi (initial)
	#
    #if (!is.null(arst)) {
    #    return(NULL)
    #} else {
        if (is.null(uvls)) {
            uvls <- matrix(nrow = lgth, ncol = 0)
            uinv <- matrix(nrow = 0, ncol = 0)
        } else {
            uinv <- helperkit::safesolve(crossprod(uvls))
        }
		xdet <- xvls - tcrossprod(uvls, crossprod(xvls, uvls) %*% uinv)
		ydet <- yvls - tcrossprod(uvls, crossprod(yvls, uvls) %*% uinv)
		xxdt <- crossprod(xdet)
		#
		augm <- tcrossprod(xvls %*% (rtfc + rnfc), lfct)
		augm <- crossprod(yvls - augm, uvls) %*% uinv
		#
		#rtfc <- array(rtfc, c(rfro, lrco, lgth)) #beta
		#
		xarr <- array(t(xvls), c(1, rfro, lgth))
		xidn <- xarr %x% diag(1, lrco)
		#
		tmpa <- xarr %x% lfct
		mscf <- helperkit::arraymultmat(tmpa, rtrm)
		mmea <- helperkit::arraymultvec(tmpa, rtrv)
		mmea <- mmea + tcrossprod(xvls %*% rnfc, lfct)
		mmea <- mmea + tcrossprod(uvls, augm)
		#
		kalm <- kalman_filter(
			yvls, mscf, mmea, evar, sscf, smea, secv, sime, sicv
		)
		fnew <- negloglike(mscf, mmea, evar, kalm$sprd, kalm$cprd)
		#
		iter <- 1
		#
		repeat {
			fold <- fnew
			#
			kalm <- kalman_smoother(
				sscf, kalm$sprd, kalm$supd, kalm$cprd, kalm$cupd
			)
			#
			sime <- kalm$ssmo[1, ]
			sicv <- helperkit::array3tomat(kalm$csmo, 1)
			sicv <- helperkit::forcesym(sicv)
			#
			ssd0 <- kalm$ssmo[1:lgth, , drop = FALSE]
			ssm0 <- colMeans(ssd0)
			ssd0 <- sweep(ssd0, 2, ssm0) # lambda0
			ssd1 <- kalm$ssmo[2:(lgth + 1), , drop = FALSE]
			ssm1 <- colMeans(ssd1)
			ssd1 <- sweep(ssd1, 2, ssm1) #lambda1
			#
			s00v <- 0
			s10v <- 0
			s11v <- 0
			scov <- 0
			#
			rtfc <- array(NA_real_, c(rfro, lrco, lgth))
			xrtv <- matrix(0, lgth, lrco)
			#
			for (time in 1:lgth) {
				tmpa <- tcrossprod(ssd0[time, ])
				s00v <- s00v + tmpa + helperkit::array3tomat(kalm$csmo, time)
				tmpa <- tcrossprod(ssd1[time, ], ssd0[time, ]) 
				s10v <- s10v + tmpa + helperkit::array3tomat(kalm$clag, time)
				tmpa <- tcrossprod(ssd1[time, ])
				#
				csmt <- helperkit::array3tomat(kalm$csmo, time + 1)
				#
				s11v <- s11v + tmpa + csmt
				#
				tmpa <- helperkit::array3tomat(xidn, time) %*% rtrm
				tmpa <- tcrossprod(tmpa %*% csmt, tmpa)
				scov <- scov + tmpa # (x' o I) HPH' (x' o I)'
				#
				rtur <- kalm$ssmo[time + 1, ]
				tmpa <- rtrm %*% rtur + rtrv #beta
				tmpa <- matrix(tmpa, rfro, lrco, TRUE)
				xrtv[time, ] <- xvls[time, , drop = FALSE] %*% tmpa
				rtfc[, , time] <- matrix(tmpa, rfro, lrco, TRUE)
			}
			#
			xrtd <- xrtv - tcrossprod(uvls, crossprod(xrtv, uvls) %*% uinv)
			#
			if (bsbv) {
				secv <- (s11v - s10v - t(s10v) + s00v)/lgth # Sigma
			} else {
				sscf <- s10v %*% helperkit::safesolve(s00v) #Pi
				smea <- drop(ssm1 - sscf %*% ssm0) #pi
				secv <- (s11v - tcrossprod(sscf, s10v))/lgth #Sigma
			}
			secv <- helperkit::forcesym(secv)
			#
			submin <- function(lfct, rnfc) {
				tmpa <- ydet - tcrossprod(xrtd + xdet %*% rnfc, lfct)
				evar <- tcrossprod(lfct %*% scov, lfct) + crossprod(tmpa)
				evar <- helperkit::forcesym(evar/lgth)
				nllk <- lgth * helperkit::ldet(evar)
				#
				return(list(nllk = nllk/lgth, evar = evar))
			}
			gnew <- submin(lfct, rnfc)
			evai <- helperkit::safesolve(gnew$evar)
			#
			ites <- 1
			#
			repeat {
				gold <- gnew
				#
				tmpa <- xrtd + xdet %*% rnfc
				tmpb <- crossprod(lfrm, (crossprod(tmpa) + scov) %x% evai)
				tmpd <- crossprod(lfrm, c(evai %*% crossprod(ydet, tmpa)))
				lfur <- helperkit::safesolve(tmpb %*% lfrm) %*% (tmpd - tmpb %*% lfrv) # gamma update
				#
				lfct <- lfrm %*% lfur + lfrv
				lfct <- matrix(lfct, lfro, lrco) # alpha update
				#
				if (ncol(rnrm) != 0) {
					tmpa <- crossprod(lfct, evai)
					tmpb <- helperkit::safesolve(tmpa %*% lfct)
					tmpb <- crossprod(rnrm, xxdt %x% tmpb)
					tmpc <- ydet - tcrossprod(xrtd, lfct)
					tmpd <- crossprod(rnrm, c(tmpa %*% crossprod(tmpc, xdet)))
					rnur <- helperkit::safesolve(tmpb %*% rnrm) %*% (tmpd - tmpb %*% rnrv) # kappa update
					#
					rnfc <- rnrm %*% rnur + rnrv
					rnfc <- matrix(rnfc, rfro, lrco, TRUE) # varrho update
				}
				#
				gnew <- submin(lfct, rnfc)
				evai <- helperkit::safesolve(gnew$evar)
		        #print(c(gnew, helperkit::ldet(evar) + ncol(yvls))))
				# == ldet(evar) ##?
				#
				if (stopcrit(gnew$nllk, gold$nllk)|| (ites > mxis)) {
					break
				}
				#
				ites <- ites + 1
			}
			#
			augm <- tcrossprod(xrtv + xvls %*% rnfc, lfct)
			augm <- crossprod(yvls - augm, uvls) %*% uinv
			#
			tmpa <- xarr %x% lfct
			mscf <- helperkit::arraymultmat(tmpa, rtrm)
			mmea <- helperkit::arraymultvec(tmpa, rtrv)
			mmea <- mmea + tcrossprod(xvls %*% rnfc, lfct)
			mmea <- mmea + tcrossprod(uvls, augm)
			#
			kalm <- kalman_filter(
				yvls, mscf, mmea, evar, sscf, smea, secv, sime, sicv
			)
			fnew <- negloglike(mscf, mmea, evar, kalm$sprd, kalm$cprd)
			#print(fnew)
			#
			if (stopcrit(fnew, fold) || (iter > mxit)) {
				break
			}
			#
			iter <- iter + 1
		}
		kalm <- kalman_smoother(
			sscf, kalm$sprd, kalm$supd, kalm$cprd, kalm$cupd
		)
		rtur <- kalm$ssmo[-1, , drop = FALSE]
	#}
	return(list(
		kalm = kalm,
		sime = sime,
		sicv = sicv,
		evar = gnew$evar, 
		sscf = sscf,
		smea = smea,
		secv = secv,
		augm = augm,
        lfct = lfct,
		rnfc = rnfc,
        rtfc = rtfc,
        lfur = lfur,
        rnur = rnur,
        rtur = rtur,
		nlik = fnew,
        rtrm = rtrm,
        rtrv = rtrv
	))
}
