#' Forecasting of VECM for I(1) Processes Under Restrictions
#'
#' Forecasts a VECM for I(1) processes under restrictions. 
#' 
#' @param tsrs Time series to forecast.
#'
#' @param objt Result of a call to `vecm`.
#' 
#' @param hriz Positive integer indicating the forecasting horizon.
#' 
#' Default is `1`.
#' 
#' @return Returns the time series to forecast 
#' extended by the forecasted values.
#' 
#' @references
#' Johansen, S. (1995). Likelihood-Based Inference in Cointegrated Vector
#' Autoregressive Models. Oxford University Press, Oxford.
#' 
#' @export
#' 
fvecm <- function(tsrs, objt, hriz = 1) {
	#tsrs <- objt$tsrs
	ordr <- objt$ordr
	dpow <- objt$dpow
	excl <- objt$excl
	bima <- objt$bima
	lfct <- objt$lfct
	rfct <- objt$rfct
	augm <- objt$augm
	#
	lgth <- nrow(tsrs)
	vals <- ncol(tsrs)
	lags <- ordr - 1
	#
	yvls <- diff(tsrs)
	#
	tsrs <- rbind(tsrs, matrix(NA_real_, hriz, ncol(tsrs)))
	yvls <- rbind(yvls, matrix(NA_real_, hriz, ncol(yvls)))
	#
	for (tidx in (lgth + 1:hriz)) {
		if (dpow > -1) {
			dvls <- tidx^(0:dpow)
			dcol <- length(dvls)
			if (excl) {
				vvls <- dvls[dcol]
				if (dcol != 1) {
					wvls <- dvls[1:(dcol - 1)]
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
		#
		chsh <- sapply(0:bima, function(indx) {
			return(chebyshev(tidx, indx, lgth))
		})
		xvls <- chsh %x% tsrs[tidx - 1, ]
		xvls <- c(xvls, vvls)
		temp <- tcrossprod(lfct, rfct) %*% xvls
		#
		if (!is.null(augm)) {
			ufor <- numeric(lags * vals)
			strt <- 1
			for (indx in seq_len(lags)) {
				fnsh <- indx * vals
				ufor[strt:fnsh] <- yvls[tidx - indx - 1, ]
				strt <- fnsh + 1
			}
			ufor <- c(ufor, wvls)
			temp <- temp + augm %*% ufor
		}
		#
		yvls[tidx - 1, ] <- temp
		tsrs[tidx, ] <- tsrs[tidx - 1, ] + yvls[tidx - 1, ]
	}
	return(tsrs) #list(tsrs = tsrs, dtsr = yvls))
}