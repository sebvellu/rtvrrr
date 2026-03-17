#' Forecasting of I(1) RS-TV-VECMs Under Restrictions
#'
#' Forecasts an I(1) RS-TV-VECM under restrictions. 
#' 
#' @param tsrs Time series to forecast.
#'
#' @param objt Result of a call to `rstvvecm`.
#' 
#' @param hriz Positive integer indicating the forecasting horizon.
#' 
#' Default is `1`.
#' 
#' @return Returns the time series to forecast 
#' extended by the forecasted values.
#' 
#' @export
#' 
frstvvecm <- function(tsrs, objt, hriz = 1) {
  #tsrs <- objt$tsrs
	ordr <- objt$ordr
	dpow <- objt$dpow
	excl <- objt$excl
	rnfc <- objt$rnfc
	lfct <- objt$lfct
	augm <- objt$augm
  rtrm <- objt$rtrm 
  rtrv <- objt$rtrv
  sscf <- objt$sscf
  smea <- objt$smea
  kalm <- objt$kalm
  #
  supd <- kalm$ssmo[nrow(kalm$ssmo), ]
	#cupd <- kalm$csmo[, , nrow(kalm$ssmo), drop = FALSE]
  #cupd <- matrix(cupd, dim(kalm$csmo)[1], dim(kalm$csmo)[2])
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
  for (tidx in lgth + 1:hriz) {
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
		xvls <- c(tsrs[tidx - 1, ], vvls)
		tmpa <- tcrossprod(lfct, rnfc) %*% xvls
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
			tmpa <- tmpa + augm %*% ufor
		}
    #
    tmpb <- t(xvls) %x% lfct
    mmea <- (tmpb %*% rtrv) + tmpa
    mscf <- tmpb %*% rtrm
    #
    sprd <- (sscf %*% supd) + smea
    #
    yvls[tidx - 1, ] <- mscf %*% sprd + mmea
    tsrs[tidx, ] <- tsrs[tidx - 1, ] + yvls[tidx - 1, ]
		#
		supd <- sprd
  }
  return(tsrs)
}