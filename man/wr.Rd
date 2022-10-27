\name{wr}
\title{Find the Response Vector and Design Matrix for a W&R Model Formula}
\alias{wr}
\usage{
wr(formula, data=NULL, expand=TRUE)
}
\description{
\code{wr} gives the response vector and design matrix for a formula in
Wilkinson and Rogers notation.
}
\arguments{
\item{formula}{A model formula.}
\item{data}{A data object or environment.}
\item{expand}{If FALSE, the covariates are read from the \code{tccov}
object without expanding to the length of the response variable.}
}
\value{
\code{wr} returns a list containing the response vector
(\code{z$response}), if included in the formula, and the design matrix
(\code{z$design}) from the data object or environment supplied or from
the global environment for the formula supplied.
}
\author{J.K. Lindsey}
\examples{
y <- rnorm(20)
x <- gl(4,5)
z <- rpois(20,2)
wr(y~x+z)
}
\keyword{programming}
