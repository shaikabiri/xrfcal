#' Convert back a matrix from additive logratio to simplex space
#'
#' This functions take a matrix X of compositions of D-1 elements in alr space and converts it
#' to simplex space with D elements with column d as the denominator.
#' @references 
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data. Journal of the Royal Statistical Society. Series B (Methodological), 44(2), 139-177. http://www.jstor.org/stable/2345821 
#' @param X A matrix of compositions in alr
#' @param d Index of column to be used as denominator
#' @return A matrix in simplex space with D columns
#' @export

invalr <- function(X,d){
  X <- exp(X)
  X <- tibble::add_column(as.data.frame(X),1,.before=d)
  X <- X/rowSums(X)
  return(X)
}