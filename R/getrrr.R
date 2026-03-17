getrrr <- function(
		yvls, xvls = NULL, vvls = NULL, wvls = NULL, lags = 0, dtrn = 0 #dtrn = 0, 1, 2, dtrn = 1 = pesaran and shin 2001, dtrn = 2 full detrending, dtrn = 0, no detrending
) {
	lgth <- nrow(yvls)
	vals <- ncol(yvls)
	#
	augm <- matrix(NA_real_, lgth - lags, lags * vals)
	#
	strt <- 1
	for (indx in seq_len(lags)) {
		fnsh <- indx * vals
		augm[, strt:fnsh] <- yvls[(lags + 1 - indx):(lgth - indx), ]
		strt <- fnsh + 1
	}
	yvls <- yvls[(lags + 1):lgth, , drop = FALSE]
	if (!is.null(xvls)) {
		xvls <- xvls[(lags + 1):lgth, , drop = FALSE]
	}
	if (!is.null(vvls)) {
		vvls <- vvls[(lags + 1):lgth, , drop = FALSE]
	}
	if (is.null(wvls)) {
		if (lags == 0) {
			uvls <- NULL
		} else {
			uvls <- augm
			if (dtrn == 2) {
				qruv <- qr(uvls[(lags + 1):lgth, ])
				yvls <- qr.resid(qruv, yvls)
				if (!is.null(xvls)) {
					xvls <- qr.resid(qruv, xvls)
				}
				if (!is.null(vvls)) {
					vvls <- qr.resid(qruv, vvls)
				}
				uvls <- NULL
			}
		}
	} else {
		if (dtrn == 1) {
			qrwv <- qr(wvls[(lags + 1):lgth, ])
			yvls <- qr.resid(qrwv, yvls)
			if (!is.null(xvls)) {
				xvls <- qr.resid(qrwv, xvls)
			}
			if (!is.null(vvls)) {
				vvls <- qr.resid(qrwv, vvls)
			}
			if (lags == 0) {
				uvls <- NULL
			} else {
				uvls <- qr.resid(qrwv, augm)
			}
		} else if (dtrn == 2) {
			qruv <- qr(cbind(augm, wvls[(lags + 1):lgth, ]))
			yvls <- qr.resid(qruv, yvls)
			if (!is.null(xvls)) {
				xvls <- qr.resid(qruv, xvls)
			}
			if (!is.null(vvls)) {
				vvls <- qr.resid(qruv, vvls)
			}
			uvls <- NULL
		} else {
			uvls <- cbind(augm, wvls[(lags + 1):lgth, ])
		}
	}
	return(list(yvls = yvls, xvls = xvls, vvls = vvls, uvls = uvls))
}
