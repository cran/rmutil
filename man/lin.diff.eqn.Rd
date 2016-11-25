\name{lin.diff.eqn}
\title{Solution of Autonomous Linear Differential Equations}
\alias{lin.diff.eqn}
\usage{
lin.diff.eqn(A, initial, t=1)
}
\description{
\code{lin.diff.eqn} numerically solves a system of autonomous linear
differential equations with given initial conditions by matrix
exponentiation.
}
\arguments{
\item{A}{A square matrix giving the coefficients of the equations.}
\item{initial}{The vector of initial values of the system.}
\item{t}{A scalar or vector of values of the independent variable for
which solutions are sought.}
}
\value{
A matrix of solutions with one row for each value of \code{t}.
}
\author{J.K. Lindsey}
\examples{
a <- matrix(c(1,0,1,0,0,0,0,0,-1),ncol=3,byrow=TRUE)
x <- c(5,7,6)
lin.diff.eqn(a,x,1)
# function giving the exact solution
exact <- function(t) c(8*exp(t)-3*exp(-t),7,6*exp(-t))
exact(1)
}
\keyword{math}
