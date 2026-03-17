ldet <- function(x) {
  rslt <- determinant(x, TRUE)$modulus
  attributes(rslt) <- NULL
  return(rslt)
}