array3tomat <- function(arry, indx) {
    return(matrix(arry[, , indx, drop = FALSE], dim(arry)[1], dim(arry)[2]))
}