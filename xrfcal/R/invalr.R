#' Additive logratio matrix to simplex matrix 
#'
#' This function takes a composition matrix X with D-1 elements in additive logratio space converts it
#' back to a matrix in simplex space with D elements with column d as the denominator of alr transformation. 
#' @examples
#' #Transform elemetal counts to alr with calcium as denominator and transform it back to simplex. 
#' Xalr <- alr(xrf$X,6) 
#' Xsimplex <- invalr(Xalr,6)
#' @references 
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data. Journal of the Royal Statistical Society. Series B (Methodological), 44(2), 139-177. http://www.jstor.org/stable/2345821 
#' @param X A matrix of compositions in alr space with D-1 columns
#' @param d Index of the column that was used as denominator for original alr transformation
#' @return A matrix of compositions in simplex space with D columns
#' @export

invalr <- function(X,d){
  X <- exp(X)
  X <- tibble::add_column(as.data.frame(X),1,.before=d)
  X <- X/rowSums(X)
  return(X)
}