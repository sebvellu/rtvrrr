arraymultvec <- function(arry, vect) {
	rslt <- apply(arry, 3, function(x) {return(x %*% vect)})
	rslt <- t(rslt)
	return(rslt)
}