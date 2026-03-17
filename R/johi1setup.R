johi1setup <- function(
		tsrs, ordr = 1, dpow = -1, excl = FALSE, citv = 0, dtrn = 0, bima = 0
) {
	lgth <- nrow(tsrs)
	if (dpow > -1) {
		#idxa <- rep(2:lgth, dpow + 1)
		#idxb <- rep(0:dpow, each = lgth - 1)
		dvls <- polydet(lgth, dpow)[2:lgth, , drop = FALSE] #matrix(idxa^idxb, lgth - 1)
		dcol <- ncol(dvls)
		if (excl) {
			vvls <- dvls[, dcol, drop = FALSE]
			#wvls <- dvls[, seq_len(dcol - 1), drop = FALSE]
			if (dcol != 1) {
				wvls <- dvls[, 1:(dcol - 1), drop = FALSE]
			} else {
				wvls <- NULL
			}
		} else {
			vvls <- NULL
			wvls <- dvls
		}
	} else {
		vvls <- NULL
		wvls <- NULL
	}
	yvls <- diff(tsrs)
	if (citv == 0) {
		xvls <- tsrs[-lgth, , drop = FALSE]
		if (bima != 0) {
			chsh <- sapply(0:bima, function(indx) {
				return(chebyshev(seq_len(lgth - 1) + 1, indx, lgth))
			})
			xvls <- t(sapply(1:(lgth - 1), function(irow) {
				return(chsh[irow, ] %x% xvls[irow, ])
			}))
		}
		rslt <- getrrr(yvls, xvls, vvls, wvls, ordr - 1, dtrn)
		return(list(
			yvls = rslt$yvls,
			xvls = cbind(rslt$xvls, rslt$vvls),
			uvls = rslt$uvls
		))
	} else {
		vals <- ncol(tsrs)
		if (citv < vals) {
			xvls <- tsrs[-lgth, 1:citv, drop = FALSE]
			vvls <- cbind(tsrs[-lgth, (citv + 1):vals, drop = FALSE], vvls)
		} else { #if (citv == ncol(tsrs)) {
			xvls <- tsrs[-lgth, , drop = FALSE]
		}
		rslt <- getrrr(yvls, xvls, vvls, wvls, ordr - 1, dtrn)
		return(list(
			yvls = rslt$yvls,
			xvls = rslt$xvls,
			vvls = rslt$vvls,
			uvls = rslt$uvls,
			lgth = lgth
		))
	}
}
