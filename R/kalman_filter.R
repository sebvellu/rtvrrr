kalman_filter <- function(
	mvls, mscf, mmea, mecv, sscf, smea, secv, sime, sicv, nllk = FALSE
) {
	mrow <- nrow(mvls)
	mcol <- ncol(mvls)
	scol <- nrow(sscf)
	#
	sprd <- matrix(NA_real_, mrow, scol)     #x_1^0,x_2^1,...,x_T^T-1
	supd <- matrix(NA_real_, mrow + 1, scol) #x_0^0,x_1^1,...,x_T^T
	#
	cprd <- array(NA_real_, c(scol, scol, mrow))     #P_1^0,P_2^1,...,P_T^T-1
	cupd <- array(NA_real_, c(scol, scol, mrow + 1)) #P_0^0,P_1^1,...,P_T^T
	#
	gain <- array(NA_real_, c(scol, mcol, mrow))
	#
	supd[1, ] <- sime
	cupd[, , 1] <- sicv
	#
	if (nllk) {
		nlik <- 0
	} else {
		nlik <- NULL
	}
	#
	for (time in 1:mrow) {
		temp <- sscf %*% supd[time, ]
		temp <- temp + smea
		sprd[time, ] <- temp
		#
		temp <- sscf %*% cupd[, , time]
		temp <- tcrossprod(temp, sscf)
		temq <- temp + secv
		cprd[, , time] <- temq
		#
		msct <- helperkit::array3tomat(mscf, time)
		cprt <- helperkit::array3tomat(cprd, time)
		#
		temp <- tcrossprod(cprt, msct)
		temq <- msct %*% temp + mecv
		temq <- helperkit::forcesym(temq)
		temr <- helperkit::safesolve(temq)
		gain[, , time] <- temp %*% temr
		#
		gait <- helperkit::array3tomat(gain, time)
		#
		inov <- mvls[time, ] - (msct %*% sprd[time, ]) - mmea[time, ]
		if (nllk) {
			#choq <- try(chol(temq), silent = TRUE)
			#if (!inherits(choq, "try-error")) {
			#	quad <- sum(backsolve(choq, inov, transpose = TRUE)^2)
			#	ldet <- 2 * sum(log(diag(choq)))
			#} else {
				quad <- as.numeric(t(inov) %*% temr %*% inov)
				ldet <- helperkit::ldet(temq)
			#}
			nlik <- nlik + 0.5 * (mcol * log(2 * pi) + ldet + quad)
		}
		temq <- sprd[time, ] + gait %*% inov
		supd[time + 1, ] <- temq
		#
		temp <- cprt - tcrossprod(gait, temp)
		cupd[, , time + 1] <- helperkit::forcesym(temp)
	}
	#
	return(list(sprd = sprd, supd = supd, cprd = cprd, cupd = cupd, nlik = nlik))
}