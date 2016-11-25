\name{fmobj}
\title{Object Finder}
\alias{fmobj}
\usage{
fmobj(z, envir=parent.frame())
}
\description{
\code{fmobj} inspects a formula and returns a list containing the
objects referred to, with indicators as to which are unknown parameters,
covariates, factor variables, and functions.
}
\arguments{
\item{z}{A model formula beginning with ~, either in Wilkinson and
Rogers notation or containing unknown parameters.}
\item{envir}{The environment in which the formula is to be
interpreted.}
}
\value{
A list, of class \code{fmobj}, containing a character vector
(\code{objects}) with the names of the objects used in a formula, and
logical vectors indicating which are unknown parameters
(\code{parameters}), covariates (\code{covariates}), factor variables
(\code{factors}), and functions (\code{functions}).
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{finterp}}
}
\examples{
x1 <- rpois(20,2)
x2 <- rnorm(20)
x3 <- gl(2,10)
#
# W&R formula
fmobj(~x1+x2+x3)
#
# formula with unknowns
fmobj(~b0+b1*x1+b2*x2)
#
# nonlinear formulae with unknowns
# log link
fmobj(~exp(b0+b1*x1+b2*x2))
}
\keyword{programming}
