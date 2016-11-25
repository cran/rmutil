\name{mexp}
\title{Matrix Exponentiation}
\alias{mexp}
\usage{
mexp(x, t=1, type="spectral decomposition", n=20, k=3)
}
\description{
\code{mexp} calculates \code{exp(t*x)} for the square matrix, \code{x}, by
spectral decomposition or series expansion.
}
\arguments{
\item{x}{A square matrix.}
\item{t}{Constant multiplying the matrix.}
\item{type}{Algorithm used: spectral decomposition or series approximation.}
\item{n}{Number of terms in the series expansion.}
\item{k}{Constant divisor to avoid over- or underflow (series
approximation only).}
}
\value{
\code{mexp} returns the exponential of a matrix.
}
\author{J.K. Lindsey}
\examples{
x <- matrix(c(1,2,3,4),nrow=2)
mexp(x)
}
\keyword{array}
