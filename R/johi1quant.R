#' @importFrom cireg get_quantile
#' 
#' @keywords internal

johi1quant <- function(
	prob, hypo, dpow = -1, excl = FALSE, lgth = 1000, simu = 10000 #hypo = nums - cint
) {
	rslt <- matrix(NA_real_, simu, hypo)
	if (dpow > - 1) {
		detr <- polydet(lgth, dpow + 1 - excl)
		#dprt <- detr %*% matrix(1/lgth^dpow, ncol(detr), hypo)
		dprt <- rep(rowSums(detr/lgth^dpow), hypo)
		dim(dprt) <- c(nrow(detr), hypo)
	} else {
		dprt <- 0
	}
	for (indx in seq_len(simu)) {
		tsrs <- dprt + apply(matrix(stats::rnorm(hypo * lgth), lgth), 2, cumsum)
		#rslt[indx, ] <- johi1stat(tsrs, 1, dpow, excl)
		for (hidx in seq_len(hypo)) {
			temp <- tsrs[, seq_len(hidx), drop = FALSE]
			rslt[indx, hypo - hidx + 1] <- johi1stat(temp, 1, dpow, excl)[1]
		}
	}
	#
	return(apply(rslt, 2, cireg::get_quantile, prob))
}
