sympow <- function(matr, powr = 1/2) {
	eign <- eigen(matr, TRUE)
	temp <- diag(eign$values^powr, ncol(matr))
	return(tcrossprod(eign$vectors %*% temp, eign$vectors))
}
