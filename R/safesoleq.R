safesoleq <- function(mata, matb) {
	tsol <- try(solve(mata, matb), TRUE)
	if (inherits(tsol, "try-error")) {
		return(mpinv(mata) %*% matb)
	} else {
		return(tsol)
	}
}