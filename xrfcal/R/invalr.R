#' Convert back a matrix from additive logratio to simplex space
#'
#' This functions take a matrix X of compositions of D-1 elements in alr space and converts it
#' to simplex space with D elements with column d as the denominator.
#' 
#' @param X A matrix of compositions in alr
#' @param d Index of column to be used as denominator
#' @return A matrix in simplex space with D columns
#' @export


#' @export

invalr <- function(X,d){
  X <- exp(X)
  X <- tibble::add_column(X,1,before=d)
  X <- X/X[,i]
}