mpinv <- function(matr) {#, tolr = 100 * .Machine$double.eps) { #tolr # tolerance for detection zero singular value #sqrt(.Machine$double.eps)
	svdm <- svd(matr)
	#pstv <- (svdm$d > tolr * svdm$d[1])
	tolr <- max(dim(matr)) * .Machine$double.eps
    pstv <- (svdm$d > tolr * svdm$d[1])
    # if (all(pstv)) {
	# 	return(tcrossprod(svdm$v %*% diag(1/svdm$d, ncol(svdm$v)), svdm$u))
	# } else if (!any(pstv)) {
	# 	return(matrix(0, ncol(matr), nrow(matr)))
	# } else {
	# 	temp <- svdm$v[, pstv, drop = FALSE] %*% diag(1/svdm$d[pstv], sum(pstv))
	# 	return(tcrossprod(temp, svdm$u[, pstv, drop = FALSE]))
	# }
    if (!any(pstv)) {
        return(matrix(0, ncol(matr), nrow(matr)))
    } else {
		temp <- svdm$v[, pstv, drop = FALSE] %*% diag(1/svdm$d[pstv], sum(pstv))
		return(tcrossprod(temp, svdm$u[, pstv, drop = FALSE]))
	}
}