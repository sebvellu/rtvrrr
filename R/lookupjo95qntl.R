#' Lookup critical values for Johansen (1995)
#'
#' Uses an internal lookup table of Johansen (1995) critical values
#' stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupjo95qntl <- function(
		prob, dpow, hypo, excl, tolr = .Machine$double.eps
) {
	tabl <- crit_jo95
	cndn <- (tabl$dpow == dpow) & (tabl$hypo == hypo)
	cndn <- cndn & (tabl$excl == excl)
	rslt <- numeric(length(prob))
	for (pidx in seq_along(prob)) {
		temp <- cndn & (abs(tabl$prob - prob[pidx]) < tolr)
		if (any(temp)) {
			rslt[pidx] <- tabl[temp, 2]
		} else {
			rslt[pidx] <- NA_real_
		}
	}
	#names(rslt) <- prob
	return(rslt)
}