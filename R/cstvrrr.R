#' Estimation of CS-TV-RRRs Under Restrictions
#'
#' Estimates an CS-TV-RRR under restrictions.
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
cstvrrr <- function(
	yvls, xvls, uvls, rrst, lnrs, ltrs, evar, arst = NULL, bsbv = FALSE,
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
			msct <- array3tomat(mscf, time)
			tmpa <- yvls[time, ] - msct %*% sprd[time, ] - mmea[time, ]
			tmpb <- tcrossprod(msct %*% array3tomat(cprd, time), msct) + mecv
			tmpa <- crossprod(tmpa, safesolve(tmpb) %*% tmpa)
			rslt <- rslt + ldet(tmpb) + tmpa
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
    #if (!is.null(arst)) {
    #    return(NULL)
    #} else {
        if (is.null(uvls)) {
            uvls <- matrix(nrow = lgth, ncol = 0)
            uinv <- matrix(nrow = 0, ncol = 0)
        } else {
            uinv <- safesolve(crossprod(uvls))
        }
        #
        # ALGORITHM:
        #
        xrfv <- xvls %*% rfct
        xdet <- xvls - tcrossprod(uvls, crossprod(xvls, uvls) %*% uinv)
        ydet <- yvls - tcrossprod(uvls, crossprod(yvls, uvls) %*% uinv)
		xxdt <- crossprod(xdet)
        #
        augm <- tcrossprod(xrfv, ltfc + lnfc)
        augm <- crossprod(yvls - augm, uvls) %*% uinv
        #
        xarr <- array(t(xrfv), c(1, lrco, lgth))
        xidn <- xarr %x% diag(1, lfro)
        #
        mscf <- arraymultmat(xidn, ltrm)
        mmea <- arraymultvec(xidn, ltrv)
        mmea <- mmea + tcrossprod(xrfv, lnfc)
        mmea <- mmea + tcrossprod(uvls, augm)
        #
        kalm <- kalman_filter(
            yvls, mscf, mmea, evar, sscf, smea, secv, sime, sicv
        )
        #print(list(mscf, secv, sicv))
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
            sicv <- array3tomat(kalm$csmo, 1)
            sicv <- forcesym(sicv)
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
            #scov <- 0
			#
            ltfc <- array(NA_real_, c(lfro, lrco, lgth))
            for (time in 1:lgth) {
                tmpa <- tcrossprod(ssd0[time, ])
                s00v <- s00v + tmpa + array3tomat(kalm$csmo, time)
                tmpa <- tcrossprod(ssd1[time, ], ssd0[time, ]) 
                s10v <- s10v + tmpa + array3tomat(kalm$clag, time)
                tmpa <- tcrossprod(ssd1[time, ]) 
                #
                csmt <- array3tomat(kalm$csmo, time + 1)
                #
                s11v <- s11v + tmpa + csmt
                #
				#tmpa <- (xrfv[time, ] %x% diag(1, lfro)) %*% ltrm
				#tmpa <- tcrossprod(tmpa %*% csmt, tmpa)
				#scov <- scov + tmpa # (x'b o I) HPH' (x'b o I)'
                #
                ltur <- kalm$ssmo[time + 1, ]
                tmpa <- ltrm %*% ltur + ltrv # alpha
                ltfc[, , time] <- matrix(tmpa, lfro, lrco)
            }
			xrlvfun <- function(xrfv) {
				rslt <- matrix(0, lgth, lfro)
				for (time in seq_len(lgth)) {
					tmpb <- array3tomat(ltfc, time)
					rslt[time, ] <- drop(tmpb %*% xrfv[time, ])
				}
				return(rslt)
			}
			xrlv <- xrlvfun(xrfv)
            xrld <- xrlv - tcrossprod(uvls, crossprod(xrlv, uvls) %*% uinv)
			#
            if (bsbv) {
                secv <- (s11v - s10v - t(s10v) + s00v)/lgth # Sigma
            } else {
                sscf <- s10v %*% safesolve(s00v) #Pi
                smea <- drop(ssm1 - sscf %*% ssm0) #pi
                secv <- (s11v - tcrossprod(sscf, s10v))/lgth #Sigma
            }
            secv <- forcesym(secv)
            #
            submin <- function(xrfv, xrld, rfct, lnfc) {
	        	tmpb <- ydet - xrld - tcrossprod(xdet %*% rfct, lnfc)
                scov <- 0
				for (time in 1:lgth) {
				    tmpa <- (xrfv[time, ] %x% diag(1, lfro)) %*% ltrm
                    csmt <- array3tomat(kalm$csmo, time + 1)
				    tmpa <- tcrossprod(tmpa %*% csmt, tmpa)
				    scov <- scov + tmpa # (x'b o I) HPH' (x'b o I)'
				}
                evar <- scov + crossprod(tmpb)
                evar <- forcesym(evar/lgth)
                nllk <- lgth * ldet(evar)
                #
                return(list(nllk = nllk/lgth, evar = evar))
            }
            gnew <- submin(xrfv, xrld, rfct, lnfc)
            evai <- safesolve(gnew$evar)
            #
            ites <- 1
            #
            repeat {
                gold <- gnew
                #
			    pmat <- array(NA_real_, c(lrco, lrco, lgth))
			    qnta <- 0
			    qntb <- 0
                for (time in 1:lgth) {
                    for(iidx in seq_len(lrco)) {
                        aidx <- (iidx - 1) * lfro + 1
                        bidx <- iidx * lfro
                        tmpa <- ltrm[aidx:bidx, , drop = FALSE] 
                        tmpa <- tmpa %*% kalm$csmo[, , time + 1]
                        for (jidx in seq_len(iidx)) {
                            cidx <- (jidx - 1) * lfro + 1
                            didx <- jidx * lfro
                            tmpb <- ltrm[cidx:didx, , drop = FALSE]
                            # Note that: sum(A * B) == sum(diag(A %*% B)) # if A and/or B are/is symmetric
                            tmpb <- sum(evai * tcrossprod(tmpa, tmpb)) # == sum(diag(evai %*% crossprod(tmpa, tmpb)))
                            pmat[iidx, jidx, time] <- tmpb
                            pmat[jidx, iidx, time] <- tmpb
                        }
                    }
                    tmpa <- array3tomat(ltfc, time) + lnfc
                    tmpb <- crossprod(tmpa, evai)
                    qnta <- qnta + (tmpb %*% tcrossprod(ydet[time, ], xdet[time, ]))
                    qntb <- qntb + (tmpb %*% tmpa + array3tomat(pmat, time))
                }
				#
                tmpb <- crossprod(rfrm, xxdt %x% qntb)
                tmpd <- crossprod(rfrm, c(qnta))
                rfur <- safesolve(tmpb %*% rfrm) %*% (tmpd - tmpb %*% rfrv) # kappa update
				#
    			rfct <- rfrm %*% rfur + rfrv
    			rfct <- matrix(rfct, rfro, lrco, TRUE) # beta update
				#
                xrfv <- xvls %*% rfct
			    xrlv <- xrlvfun(xrfv)
                xrld <- xrlv - tcrossprod(uvls, crossprod(xrlv, uvls) %*% uinv)
				#
                if (ncol(lnrm) != 0) {
                    tmpb <- crossprod(rfct, xxdt %*% rfct)
                    tmpb <- crossprod(lnrm, tmpb %x% evai)
                    tmpd <- evai %*% crossprod(ydet - xrld, xdet %*% rfct)
                    tmpd <- crossprod(lnrm, c(tmpd))
                    lnur <- safesolve(tmpb %*% lnrm) %*% (tmpd - tmpb %*% lnrv) # gamma update
                    #
                    lnfc <- lnrm %*% lnur + lnrv
                    lnfc <- matrix(lnfc, lfro, lrco) # alpha update
                }
                #
                gnew <- submin(xrfv, xrld, rfct, lnfc)
                evai <- safesolve(gnew$evar)
                #print(c(gnew, ldet(evar) + ncol(yvls))))
				# == ldet(evar) ##?
                #
                if (stopcrit(gnew$nllk, gold$nllk) || (ites > mxis)) {
                    break
                }
                #
                ites <- ites + 1
			}
            #
			augm <- xrlv - tcrossprod(xrfv, lnfc)
			augm <- crossprod(yvls - augm, uvls) %*% uinv
			#
			xarr <- array(t(xrfv), c(1, lrco, lgth))
            xidn <- xarr %x% diag(1, lfro)
			#
			mscf <- arraymultmat(xidn, ltrm)
			mmea <- arraymultvec(xidn, ltrv)
			mmea <- mmea + tcrossprod(xrfv, lnfc)
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
        ltur <- kalm$ssmo[-1, , drop = FALSE]
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
        rfct = rfct,
		lnfc = lnfc,
        ltfc = ltfc,
        rfur = rfur,
        lnur = lnur,
        ltur = ltur,
		nlik = fnew,
        ltrm = ltrm,
        ltrv = ltrv
	))
}
