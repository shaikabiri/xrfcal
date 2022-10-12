#' Convert a matrix to additive logratio matrix
#'
#' This functions take a matrix X of compositions of D elements, close it and convert it
#' to additive logratio space with D-1 elements with column d as the denominator.
#' 
#' @param X A matrix of compositions 
#' @param d Index of column to be used as denominator
#' @return An alr matrix of D-1 elements
#' @export


#' @export

alr <- function(X,d){
  X <- X/rowSums(X)
  X <- X/X[,d]
  X <- X[,-d]
  X <- log(X)
  return(X)
}