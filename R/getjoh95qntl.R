getjoh95qntl <- function(
	prob, dpow, hypo, excl, lgth = 1000, simu = 10000,
	tolr = .Machine$double.eps
) {
	temp <- try(sapply(rev(seq_len(hypo)), function(hypo) {
		return(lookupjo95qntl(
			prob = prob,
			dpow = dpow,
			hypo = hypo,
			excl = ifelse(dpow == -1, FALSE, excl),
			tolr = tolr
		))
	}), TRUE)
	if (!inherits(temp, "try-error")) {
        if (!any(is.na(temp))) {
            return(temp)
        }
	}
	return(johi1quant(prob, hypo, dpow, excl, lgth, simu))
}