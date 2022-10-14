#' Simplex matrix to additive logratio matrix
#'
#' This function takes a composition matrix X with D elements in simplex space, closes it 1 and converts it
#' to a matrix in additive logratio space with D-1 elements with column d as the denominator.
#' @references 
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data. Journal of the Royal Statistical Society. Series B (Methodological), 44(2), 139-177. http://www.jstor.org/stable/2345821 
#' @param X A matrix of compositions in simplex space
#' @param d Index of the column to be used as denominator
#' @return An alr matrix of D-1 elements
#' @export
alr <- function(X,d){
  X <- X/rowSums(X)
  X <- X/X[,d]
  X <- X[,-d]
  X <- log(X)
  return(X)
}

