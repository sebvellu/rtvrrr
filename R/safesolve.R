safesolve <- function(mata, matb = NULL) {
	if (is.null(matb)) {
		# #return(solve(mata))
		# tsol <- try(solve(mata), TRUE)
		# if (inherits(tsol, "try-error")) {
		# 	return(mpinv(mata))
		# } else {
		# 	return(tsol)
		# }
		return(safesolid(mata))
	} else {
		# #return(solve(mata, matb))
		# tsol <- try(solve(mata, matb), TRUE)
		# if (inherits(tsol, "try-error")) {
		# 	return(mpinv(mata) %*% matb)
		# } else {
		# 	return(tsol)
		# }
		return(safesoleq(mata, matb))
	}
}
