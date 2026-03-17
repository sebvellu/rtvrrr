polydet <- function(lgth, dpow, optn = 1) {
	if (dpow == -1) {
		dvls <- NULL
	} else if (optn == 1) {
		idxa <- rep(1:lgth, dpow + 1)
		idxb <- rep(0:dpow, each = lgth)
		dvls <- matrix(idxa^idxb, lgth)
	} else if (optn == 2) {
		idxa <- rep(1:lgth/lgth, dpow + 1)
		idxb <- rep(0:dpow, each = lgth)
		dvls <- matrix(idxa^idxb, lgth)
	} else { #if (optn == 3) {
		idxa <- rep((1:lgth - (lgth + 1)/2)/lgth, dpow + 1)
		idxb <- rep(0:dpow, each = lgth)
		dvls <- matrix(idxa^idxb, lgth)
	}
	return(dvls)
}
