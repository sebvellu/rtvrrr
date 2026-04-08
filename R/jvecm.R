#' Estimation of VECM for I(1) Processes Using Johansen (1995) Methodology
#'
#' Estimates a VECM for I(1) processes using the Johansen (1995) methodology.
#' 
#' @param tsrs Matrix of values for multivariate I(1) process.
#' 
#' @param cirk A positive integer specifying the cointegrating rank.
#' 
#' @param ordr Order of the VECM. Equal to one plus the number of lags.
#' 
#' Default is `1`.
#' 
#' @param dpow Highest power of polynomial time trends (e.g., -1: 
#' no deterministic components, 0: constant only, 1: constant and linear trend)
#' 
#' Default is `0`.
#' 
#' @param excl Logical value indicating whether the highest power 
#' of polynomial time trends is restricted to the long-run 
#' relationship.
#' 
#' Default is `FALSE`.
#' 
#' @return Returns a list of the following values:
#' 
#'   - `cint`: Estimated long-run coefficient matrix
#'   - `ldng`: Estimated loading/adjument coefficient matrix
#'   - `rvar`: Estimated variance-covariance matrix of the error term
#'   - `rrkc`: Matrix product `ldng` * `cint`'
#'   - `agmt`: Estimated coefficient matrix of the short-run dynamics
#'   - `cior`: Orthogonal complement of `cint`
#'   - `ldor`: Orthogonal complement of `ldng`
#' 
#' @references
#' Johansen, S. (1995). Likelihood-Based Inference in Cointegrated Vector
#' Autoregressive Models. Oxford University Press, Oxford.
#' 
#' @export
#' 
jvecm <- function(tsrs, cirk, ordr = 1, dpow = 0, excl = FALSE) {
	excl <- ifelse(dpow == -1, FALSE, excl)
	# johansen 1991, 1988, 1995 I(1)
	#lgth <- nrow(tsrs) #OPTION 1
	rslt <- johi1setup(tsrs, ordr, dpow, excl, 0, 0, 0)#FALSE) #dtrn)
	lgth <- nrow(rslt$yvls) #OPTION 2
	rslt <- rrrstuff(rslt$yvls, rslt$xvls, rslt$uvls)
	nums <- length(rslt$eigo$values)
	#
	eign <- crossprod(sqrt(lgth) * rslt$s11s, rslt$eigo$vectors) # Johansen normalized eigenvectors
	#eign <- crossprod(rslt$eigo$vectors, (rslt$s11v/lgth) %*% rslt$eigo$vectors)
	#eign <- rslt$eigo$vectors %*% solve(roth(ncol(eign)) %*% chol(eign)) #roth(ncol(eign))
	#print(t(eign) %*% (rslt$s11v/lgth) %*% eign)
	cint <- eign[, 1:cirk, drop = FALSE]
	#
	temp <- crossprod(cint, rslt$s11v %*% cint)
	ldng <- (rslt$s01v %*% cint) %*% helperkit::safesolve(temp) #rslt$s01v %*% cint/sqrt(lgth) #*sqrt(lgth)
	#print(crossprod(ldng, solve(rslt$s00v/lgth) %*% ldng))
	#print(rslt$eigo$values)
	#print(((rslt$s01v %*% rslt$eigo$vectors %*% safesolve(crossprod(rslt$eigo$vectors , rslt$s11v %*% rslt$eigo$vectors))) %*% t(rslt$eigo$vectors)))
	#print(ldng %*% t(cint))
	#stop()
	#
	rvar <- (rslt$s00v - tcrossprod(ldng %*% temp, ldng))/lgth
	#temp <- tcrossprod(ldng, rslt$s01v %*% cint)
	#temq <- crossprod(cint, rslt$s11v %*% cint)
	#rvar <- (rslt$s00v - temp - t(temp) + tcrossprod(ldng %*% temq, ldng))/lgth
	#
	rrkc <- tcrossprod(ldng, cint)
	#
	if (!is.null(rslt$m22v)) { #(dpow != 0) & (order != 1)
		agmt <- rslt$m02v %*% rslt$m22i - rrkc %*% rslt$m12v %*% rslt$m22i
	} else {
		agmt <- NULL
	}
	#
	eior <- eign[, seq_len(nums - cirk) + cirk, drop = FALSE]#/sqrt(lgth)
	cior <- (rslt$s11v %*% eior)/lgth
	ldor <- (rslt$s00i %*% rslt$s01v) %*% eior
	#print(crossprod(ldng, ldor))
	#print(crossprod(cint, cior))
	#stop()
	#
	return(list(
		cint = cint, # beta
		ldng = ldng, # alpha
		rvar = rvar, # Omega
		rrkc = rrkc, # alpha * beta'
		agmt = agmt, # Psi
		cior = cior, # beta_perp
		ldor = ldor  # alpha_perp
	))
}