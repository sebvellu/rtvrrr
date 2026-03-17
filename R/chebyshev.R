chebyshev <- function(time, ordr = 0, lgth = NULL) {
	if (is.null(lgth)) {
		lgth <- length(time)
	}
	if (ordr == 0) {
		return(rep(1, length(time)))
	} else {
		return(sqrt(2) * cos(ordr * pi * (time - 0.5)/lgth))
	}
}