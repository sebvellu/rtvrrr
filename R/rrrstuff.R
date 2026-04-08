rrrstuff <- function(yvls, xvls, uvls = NULL, lgth) {
	if (is.null(uvls)) {
		m00v <- crossprod(yvls)
		m01v <- crossprod(yvls, xvls)
		m11v <- crossprod(xvls)
		m02v <- NULL
		m12v <- NULL
		m22v <- NULL
		m22i <- NULL
		#
		s00v <- m00v
		s01v <- m01v
		s11v <- m11v
		#
		s00i <- helperkit::safesolve(s00v)
		s11i <- helperkit::safesolve(s11v)
		s11s <- helperkit::sympow(s11i) #sympow(s11v, -1/2)
	} else {
		m00v <- crossprod(yvls)
		m01v <- crossprod(yvls, xvls)
		m02v <- crossprod(yvls, uvls)
		m11v <- crossprod(xvls)
		m12v <- crossprod(xvls, uvls)
		m22v <- crossprod(uvls)
		#
		m22i <- helperkit::safesolve(m22v)
		#
		s00v <- m00v - tcrossprod(m02v %*% m22i, m02v)
		s01v <- m01v - tcrossprod(m02v %*% m22i, m12v)
		s11v <- m11v - tcrossprod(m12v %*% m22i, m12v)
		#
		s00i <- helperkit::safesolve(s00v)
		s11i <- helperkit::safesolve(s11v)
		s11s <- helperkit::sympow(s11i)
	}
	temp <- crossprod(s01v, s00i %*% s01v)
	eigo <- eigen(tcrossprod(s11s %*% temp, s11s))
	#print(sapply(eigo$values, function(x) {det(x * s11v - crossprod(s01v, s00i %*% s01v))}))
	#
	rslt <- list(
		yvls = yvls,
		xvls = xvls,
		uvls = uvls,
		s00v = s00v,
		s00i = s00i,
		s01v = s01v,
		s11v = s11v,
		s11i = s11i,
		s11s = s11s,
		m02v = m02v,
		m12v = m12v,
		m22v = m22v,
		m22i = m22i,
		eigo = eigo
	)
	return(rslt)
}