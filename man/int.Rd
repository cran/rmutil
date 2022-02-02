\name{int}
\title{Vectorized Numerical Integration}
\alias{int}
\usage{
int(f, a=-Inf, b=Inf, type="Romberg", eps=0.0001, max=NULL, d=NULL, p=0)
}
\description{
\code{int} performs numerical integration of a given function using
either Romberg integration or algorithm 614 of the collected
algorithms from ACM. Only the former is vectorized. The latter uses
formulae optimal in certain Hardy spaces h(p,d).

Functions may have singularities at one or both end-points of the
interval (a,b).
}
\arguments{
\item{f}{The function (of one variable) to integrate, returning either
a scalar or a vector.}
\item{a}{A scalar or vector (only Romberg) giving the lower bound(s).
A vector cannot contain both -Inf and finite values.}
\item{b}{A scalar or vector (only Romberg) giving the upper bound(s).
A vector cannot contain both Inf and finite values.}
\item{type}{The algorithm to be used, by default Romberg integration.
Otherwise, it uses the TOMS614 algorithm.}
\item{eps}{Precision.}
\item{max}{For Romberg, the maximum number of steps, by default set
to 16. For TOMS614, the maximum number of function evaluations, by
default set to 100.}
\item{d}{For Romberg, the number of extrapolation points so that
2d is the order of integration, by default set to 5; d=2 is Simpson's
rule. For TOMS614, heuristic termination = any real number;
deterministic termination = a number in the range 0 < d < pi/2
by default, set to 1.}
\item{p}{For TOMS614, p = 0: heuristic termination, p = 1:
deterministic termination with the infinity norm, p > 1: deterministic
termination with the p-th norm.}
}
\value{
The vector of values of the integrals of the function supplied.
}
\author{J.K. Lindsey}
\references{
ACM algorithm 614 appeared in

ACM-Trans. Math. Software, Vol.10, No. 2, Jun., 1984, p. 152-160.

See also

Sikorski,K., Optimal quadrature algorithms in HP spaces, Num. Math.,
39, 405-410 (1982).
}
\examples{
f <- function(x) sin(x)+cos(x)-x^2
int(f, a=0, b=2)
int(f, a=0, b=2, type="TOMS614")
#
f <- function(x) exp(-(x-2)^2/2)/sqrt(2*pi)
int(f, a=0:3)
int(f, a=0:3, d=2)
1-pnorm(0:3, 2)
#
f <- function(x) dnorm(x)
int(f, a=-Inf, b=qnorm(0.975))
int(f, a=-Inf, b=qnorm(0.975), type="TOMS614", max=1e2)
}
\keyword{math}
