johi1stat <- function(
	tsrs, ordr = 1, dpow = 0, excl = FALSE#, dtrn = 0
) { #, alth = NULL ##### cirk = 1
	# johansen 1991, 1988, 1995 I(1)
	#lgth <- nrow(tsrs) #OPTION 1
	rslt <- johi1setup(tsrs, ordr, dpow, excl, 0, 0, 0) #0, 1, 0
	lgth <- nrow(rslt$yvls) #OPTION 2
	rslt <- rrrstuff(rslt$yvls, rslt$xvls, rslt$uvls)
	nums <- length(rslt$eigo$values)
	#
	rslt$eigo$values <- sort(
		x = pmin(pmax(Re(rslt$eigo$values), 0), 1 - 1e-12),
		decreasing = TRUE
	)
	#
	alth <- nums - excl # length(rslt$eigo$values)
	stat <- function(crnk) {
		return(-lgth * sum(log(1 - rslt$eigo$values[(crnk + 1):alth])))
	}
	return(sapply(0:(alth - 1), stat))
}
