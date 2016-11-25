\name{int2}
\title{Vectorized Two-dimensional Numerical Integration}
\alias{int2}
\usage{
int2(f, a=c(-Inf,-Inf), b=c(Inf,Inf), eps=1.0e-6, max=16, d=5)
}
\description{
\code{int} performs vectorized numerical integration of a given
two-dimensional function.
}
\arguments{
\item{f}{The function (of two variables) to integrate, returning either
a scalar or a vector.}
\item{a}{A two-element vector or a two-column matrix giving the lower
bounds. It cannot contain both -Inf and finite values.}
\item{b}{A two-element vector or a two-column matrix giving the upper
bounds. It cannot contain both Inf and finite values.}
\item{eps}{Precision.}
\item{max}{The maximum number of steps, by default set to 16.}
\item{d}{The number of extrapolation points so that 2k is the order of
integration, by default set to 5; d=2 is Simpson's rule.}
}
\value{
The vector of values of the integrals of the function supplied.
}
\author{J.K. Lindsey}
\examples{
f <- function(x,y) sin(x)+cos(y)-x^2
int2(f, a=c(0,1), b=c(2,4))
#
fn1 <- function(x, y) x^2+y^2
fn2 <- function(x, y) (1:4)*x^2+(2:5)*y^2
int2(fn1, c(1,2), c(2,4))
int2(fn2, c(1,2), c(2,4))
int2(fn1, matrix(c(1:4,1:4),ncol=2), matrix(c(2:5,2:5),ncol=2))
int2(fn2, matrix(c(1:4,1:4),ncol=2), matrix(c(2:5,2:5),ncol=2))
}
\keyword{math}
