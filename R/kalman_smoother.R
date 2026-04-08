kalman_smoother <- function(sscf, sprd, supd, cprd, cupd) {
    mrow <- nrow(sprd)
	scol <- nrow(sscf)
	#
	ssmo <- matrix(NA_real_, mrow + 1, scol) #x_0^T,x_1^T,...,x_T^T
    #
	csmo <- array(NA_real_, c(scol, scol, mrow + 1)) #P_0^T,P_1^T,...,P_T^T
	#
	clag <- array(NA_real_, c(scol, scol, mrow)) #P_1,0^T,...,P_T,T-1^T
	#
	ssmo[mrow + 1, ] <- supd[mrow + 1, ]
	csmo[, , mrow + 1] <- cupd[, , mrow + 1]
	#
	temr <- matrix(0, scol, scol)
	#
	for (time in mrow:1) {
		cupt <- helperkit::array3tomat(cupd, time)
		cprt <- helperkit::array3tomat(cprd, time)
		#
		temp <- tcrossprod(cupt, sscf)
		temp <- temp %*% helperkit::safesolve(cprt)
		#
		temq <- ssmo[time + 1, ] - sprd[time, ]
		temq <- supd[time, ] + temp %*% temq
		ssmo[time, ] <- temq
		#
		temq <- helperkit::array3tomat(csmo, time + 1) - cprt
		temq <- cupt + tcrossprod(temp %*% temq, temp)
		csmo[, , time] <- temq
		#
		tems <- tcrossprod(helperkit::array3tomat(cupd, time + 1) + temr, temp)
		clag[, , time] <- tems
		#
		temr <- sscf %*% cupt
		temr <- temp %*% (helperkit::array3tomat(clag, time) - temr)
	}
	#
	return(list(ssmo = ssmo, csmo = csmo, clag = clag))
}