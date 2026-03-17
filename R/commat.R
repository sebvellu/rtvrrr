commat <- function(m, n) {
	if (m == 0 | n == 0) {
		return(matrix(nrow = 0, ncol = 0))
	}
	rslt <- matrix(1:(m * n), n, m, byrow = TRUE)
	return(diag(1, m * n)[c(rslt), , drop = FALSE])
}
