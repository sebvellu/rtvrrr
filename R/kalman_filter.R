kalman_filter <- function(mvls, mscf, mmea, mecv, sscf, smea, secv, sime, sicv) {
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
		msct <- array3tomat(mscf, time)
		cprt <- array3tomat(cprd, time)
		#
		temp <- tcrossprod(cprt, msct)
		temq <- msct %*% temp + mecv
		gain[, , time] <- temp %*% safesolve(temq)
		#
		gait <- array3tomat(gain, time)
		#
		temq <- mvls[time, ] - (msct %*% sprd[time, ]) - mmea[time, ]
		temq <- sprd[time, ] + gait %*% temq
		supd[time + 1, ] <- temq
		#
		temp <- cprt - tcrossprod(gait, temp)
		cupd[, , time + 1] <- forcesym(temp)
	}
	#
	return(list(sprd = sprd, supd = supd, cprd = cprd, cupd = cupd))
}