arraymultmat <- function(arry, matr) {
	rslt <- array(NA_real_, c(dim(arry)[1], ncol(matr), dim(arry)[3]))
	rslt[] <- apply(arry, 3, function(x) {return(x %*% matr)})
	return(rslt)
}