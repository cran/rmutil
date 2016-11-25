\name{runge.kutta}
\title{Runge-Kutta Method for Solving Differential Equations}
\alias{runge.kutta}
\usage{
runge.kutta(f, initial, x)
}
\description{
\code{runge.kutta} numerically solves a differential equation by the
fourth-order Runge-Kutta method.
}
\arguments{
\item{f}{A function \code{dy/dx=func(y,x)}.}
\item{initial}{The initial value of \code{y}.}
\item{x}{A vector of values of \code{x} for which the values or
\code{y} are required.}
}
\value{
A vector of values of \code{y} as solution of the function \code{f}
corresponding to the values in \code{x}.
}
\author{J.K. Lindsey}
\examples{
fn <- function(y,x) (x*y-y^2)/x^2
soln <- runge.kutta(fn,2,seq(1,3,by=1/128))
## exact solution
exact <- seq(1,3,by=1/128)/(0.5+log(seq(1,3,by=1/128)))
rbind(soln, exact)
}
\keyword{math}
