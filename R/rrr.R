#' Estimation of a Reduced Rank Regressions Under Restrictions
#'
#' Estimates a reduced rank regression under restrictions.
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
#' @param lrst Specification of the restrictions on the column space
#' parameter matrix. Either `NULL` or a list. Must have length 4 if 
#' `rrst` is `NULL`, and length 3 otherwise.
#' 
#'   - If list of length 4, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters, and fourth element initial value 
#'   for the row space parameter matrix.
#'   - If list of length 3, then first element restriction matrix, 
#'   second element restriction vector, third element vector of initial values
#'   for the free parameters.
#' 
#' @param arst Specification of the restrictions on the coefficients of 
#' the input vectors in measurement equation. If not `NULL`, `NULL` is returned.
#' 
#' Default is `NULL`.
#' 
#' @param fish Should fisher information matrix be returned?
#' 
#' Default is `FALSE`.
#' 
#' @param tolr Tolerance limit for optimization. 
#' 
#' Default is 100 * sqrt(.Machine$double.eps). 
#' 
#' @param mxit Maximum number of iterations.
#' 
#' Default is Inf.
#' 
#' @return Returns `NULL` if `arst` is not `NULL`, otherwise 
#' a list of the following values:
#' 
#'   - `lfct`: Estimated column space parameter matrix
#'   - `rfct`: Estimated rows space parameter matrix
#'   - `lfun`: Estimated unrestricted column space parameter vector
#'   - `rfun`: Estimated unrestricted row space parameter vector
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
#' 
#' @references
#' Hansen, P. R. (2003). Generalized Reduced Rank Regression. Working Paper,
#' Brown University.
#' 
#' @export
#' 
rrr <- function(#restricted reduced rank regression estimation
	yvls, xvls, uvls, rrst, lrst, arst = NULL, fish = FALSE,
	tolr = 100 * sqrt(.Machine$double.eps), mxit = Inf
) {
    if (!is.null(arst)) {
        return(NULL)
    }
	#
	lgth <- nrow(yvls)
	#
	# objtfun <- function(lfct, rfct) {
	# 	tmpa <- tcrossprod(lfct, s01v %*% rfct)
	# 	tmpb <- crossprod(rfct, s11v %*% rfct)
	# 	tmpa <- s00v - tmpa - t(tmpa) + tcrossprod(lfct %*% tmpb, lfct)
	# 	evar <- forcesym(tmpa/lgth)
	# 	nllk <- lgth * ldet(evar) #+ sum(diag(safesolve(evar) %*% tmpa))
	# 	return(list(nllk = nllk/lgth, evar = evar))
	# }
	stopcrit <- function(fnew, fold) {
		return(abs(fnew - fold) < tolr * max(1, abs(fold)))
	}
	#
	rfro <- ncol(xvls)
	lfro <- ncol(yvls)
	#
	if (is.null(lrst)) {
		if (is.null(rrst)) {
			return(NULL)
		} else {
			rfrm <- rrst[[1]]
			rfrv <- rrst[[2]]
			rfur <- rrst[[3]] # kappa
			#lfct <- rrst[[4]]
			#
			lrco <- as.integer(length(rfrv)/rfro) # ncol(lfct)
			#
			lfrm <- diag(1, lfro * lrco)
			lfrv <- numeric(lfro * lrco)
			lfur <- numeric(lfro * lrco) # gamma
		}
	} else {
		lfrm <- lrst[[1]]
		lfrv <- lrst[[2]]
		lfur <- lrst[[3]] # gamma
		#
		lrco <- as.integer(length(lfrv)/lfro)
		if (is.null(rrst)) {
			rfct <- lrst[[4]]
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
	#
	rfct <- rfrm %*% rfur + rfrv
	rfct <- matrix(rfct, rfro, lrco, TRUE)
	#
	lfct <- lfrm %*% lfur + lfrv
	lfct <- matrix(lfct, lfro, lrco)
	#
    #if (!is.null(arst)) {
    #    return(NULL)
    #} else {
		if (is.null(uvls)) {
			m00v <- crossprod(yvls)
			m01v <- crossprod(yvls, xvls)
			m11v <- crossprod(xvls)
			m02v <- NULL
			m12v <- NULL
			m22v <- NULL
			m22i <- NULL
			#
			s00v <- m00v
			s01v <- m01v
			s11v <- m11v
		} else {
			m00v <- crossprod(yvls)
			m01v <- crossprod(yvls, xvls)
			m02v <- crossprod(yvls, uvls)
			m11v <- crossprod(xvls)
			m12v <- crossprod(xvls, uvls)
			m22v <- crossprod(uvls)
			#
			m22i <- safesolve(m22v)
			#
			s00v <- m00v - tcrossprod(m02v %*% m22i, m02v)
			s01v <- m01v - tcrossprod(m02v %*% m22i, m12v)
			s11v <- m11v - tcrossprod(m12v %*% m22i, m12v)
		}
		#
		objtfun <- function(lfct, rfct) {
			tmpa <- tcrossprod(lfct, s01v %*% rfct)
			tmpb <- crossprod(rfct, s11v %*% rfct)
			tmpa <- s00v - tmpa - t(tmpa) + tcrossprod(lfct %*% tmpb, lfct)
			evar <- forcesym(tmpa/lgth)
			nllk <- lgth * ldet(evar) #+ sum(diag(safesolve(evar) %*% tmpa))
			return(list(nllk = nllk/lgth, evar = evar))
		}
		#
		fnew <- objtfun(lfct, rfct)
		evai <- safesolve(fnew$evar)
		#
		iter <- 1
		repeat {
			fold <- fnew
			#
			tmpb <- crossprod(lfct, evai)
			tmpa <- crossprod(rfrm, s11v %x% (tmpb %*% lfct))
			tmpb <- crossprod(rfrm, as.vector(tmpb %*% s01v))
			rfur <- safesolve(tmpa %*% rfrm) %*% (tmpb - tmpa %*% rfrv)
			#
			rfct <- rfrm %*% rfur + rfrv
			rfct <- matrix(rfct, rfro, lrco, TRUE)
			#
			tmpa <- crossprod(lfrm, crossprod(rfct, s11v %*% rfct) %x% evai)
			tmpb <- crossprod(lfrm, as.vector(evai %*% s01v %*% rfct))
			lfur <- safesolve(tmpa %*% lfrm) %*% (tmpb - tmpa %*% lfrv)
			#
			lfct <- lfrm %*% lfur + lfrv
			lfct <- matrix(lfct, lfro, lrco)
			#
			fnew <- objtfun(lfct, rfct)
			evai <- safesolve(fnew$evar)
			#
			if (stopcrit(fnew$nllk, fold$nllk) || (iter > mxit)) {
				break
			}
			#
			iter <- iter + 1
		}
		rsds <- yvls - tcrossprod(xvls %*% rfct, lfct)
		if (!is.null(m22i)) {
			augm <- (m02v - tcrossprod(lfct, rfct) %*% m12v) %*% m22i
			rsds <- rsds - tcrossprod(uvls, augm)
		} else {
			augm <- NULL
			uvls <- matrix(nrow = nrow(xvls), ncol = 0)
			m12v <- matrix(nrow = ncol(xvls), ncol = 0)
			m22v <- matrix(nrow = 0, ncol = 0)
		}
		#
		#evai <- safesolve(fnew$evar)
		flvr <- m11v %x% crossprod(lfct, evai %*% lfct)
		flvr <- crossprod(rfrm, flvr %*% rfrm)
		#
		#Fisher matrix
		if (fish) {
			fish <- flvr
			tmpa <- diag(1, ncol(rfct)) %x% (crossprod(xvls, rsds) %*% evai)
			tmpa <- commat(nrow(rfct), ncol(rfct)) %*% tmpa
			tmpa <- ((m11v %*% rfct) %x% crossprod(lfct, evai)) - tmpa
			tmpa <- crossprod(rfrm, tmpa %*% lfrm)
			tmpb <- crossprod(rfrm, m12v %x% crossprod(lfct, evai))
			tmpa <- evai %x% evai
			tmpc <- t(lfct) %x% crossprod(xvls, rsds)
			tmpc <- commat(nrow(lfct), ncol(lfct)) %*% tmpc
			tmpd <- crossprod(xvls, rsds) %x% t(lfct)
			tmpc <- (crossprod(rfrm, tmpc + tmpd) %*% tmpa)/2
			fish <- cbind(fish, tmpa, tmpb, tmpc)
			tmpd <- crossprod(rfct, m11v %*% rfct) %x% evai
			tmpd <- crossprod(lfrm, tmpd %*% lfrm)
			tmpe <- crossprod(lfrm, crossprod(rfct, m12v) %x% evai)
			tmpf <- diag(1, nrow(lfct)) %x% crossprod(rfct, crossprod(xvls, rsds))
			tmpf <- crossprod(commat(nrow(lfct), ncol(lfct)), tmpf)
			tmpg <- crossprod(rfct, crossprod(xvls, rsds)) %x% diag(1, nrow(lfct))
			tmpf <- (crossprod(lfrm, tmpf + tmpg) %*% tmpa)/2
			fish <- rbind(fish, cbind(t(tmpa), tmpd, tmpe, tmpf))
			tmpg <- m22v %x% evai
			tmph <- diag(1, nrow(lfct)) %x% crossprod(uvls, rsds)
			tmph <- crossprod(commat(nrow(lfct), ncol(uvls)), tmph)
			tmpi <- crossprod(uvls, rsds) %x% diag(1, nrow(lfct))
			tmph <- ((tmph + tmpi) %*% tmpa)/2
			fish <- rbind(fish, cbind(t(tmpb), t(tmpe), tmpg, tmph))
			tmpi <- evai %*% crossprod(rsds)
			tmpi <- (tmpi %x% diag(1, nrow(lfct))) + (diag(1, nrow(lfct)) %x% tmpi)
			tmpi <- ((tmpi - diag(1, nrow(lfct)^2)) %*% tmpa)/2
			fish <- rbind(fish, cbind(t(tmpc), t(tmpf), t(tmph), tmpi))
			#
			fish <- forcesym(fish)
			#
			unvr <- safesolve(fish)
			idnt <- (qr(fish)$rank == ncol(fish))
		} else {
			fish <- NULL
			unvr <- NULL
			idnt <- NULL
		}
		#
		#lvar <- safesolve(flvr)
	#}
	#
	return(list(
		lfct = lfct,
		rfct = rfct,
		lfun = lfur,
		rfun = rfur,
		evar = fnew$evar,
		augm = augm,
		lrco = lrco,
		fish = fish,
		unvr = unvr,
		idnt = idnt,
		nllk = fnew$nllk#,
		#lvar = lvar
	))
}