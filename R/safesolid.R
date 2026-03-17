safesolid <- function(matr) {
	tsol <- try(solve(matr), TRUE)
	if (inherits(tsol, "try-error")) {
		return(mpinv(matr))
	} else {
		return(tsol)
	}
}