% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alr.R
\name{alr}
\alias{alr}
\title{Simplex matrix to additive logratio matrix}
\usage{
alr(X, d)
}
\arguments{
\item{X}{A matrix of compositions in simplex space with D columns. Doesn't need to be closed.}

\item{d}{Index of the column to be used as denominator}
}
\value{
A matrix of compositions in alr space with D-1 columns
}
\description{
This function takes a composition matrix X with D elements in simplex space, closes it to 1 and converts it
to a matrix in additive logratio space with D-1 elements with column d as the denominator.
}
\examples{
#Transform elemetal counts to alr with calcium as denominator 
Xalr <- alr(xrf$X,6)
}
\references{
Aitchison, J. (1982). The Statistical Analysis of Compositional Data. Journal of the Royal Statistical Society. Series B (Methodological), 44(2), 139-177. http://www.jstor.org/stable/2345821
}
