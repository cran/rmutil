\name{int}
\title{Numerical Integration of a Function}
\usage{
int(f, a="-infty", b="infty", type="Romberg", eps=1.0e-6, max, d, p=0)
}
\alias{int}
\arguments{
\item{f}{The function (of one variable) to integrate, returning either
a scalar or a vector.}
\item{a}{A scalar or vector giving the lower bound. If non-numeric,
taken to be -infty, in which case it must be the same for the whole
vector.}
\item{b}{A scalar or vector giving the upper bound. If non-numeric,
taken to be infty, in which case it must be the same for the whole
vector.}
\item{type}{The algorithm to be used, by default Romberg integration.
Otherwise, it uses the TOMS614 algorithm.}
\item{eps}{Precision.}
\item{max}{For Romberg, the maximum number of steps, by default set
to 16. For TOMS614, the maximum number of function evaluations, by
default set to 100.}
\item{d}{For Romberg, the number of extrapolation points so that
2k is the order of integration, by default set to 5; d=2 is Simpson's
rule. For TOMS614, heuristic termination = any real number;
deterministic termination = a number in the range 0 < d < pi/2
by default, set to 1.}
\item{p}{ For TOMS614, p = 0: heuristic termination, p = 1:
deterministic termination with the infinity norm, p > 1: deterministic
termination with the p-th norm.}
}
\description{
\code{int} performs numerical integration of a given function using
either Romberg integration or algorithm 614 of the collected
algorithms from ACM. Only the former is vectorized. The latter
appeared in ACM-Trans. Math. Software, Vol.10, No. 2, Jun., 1984,
p. 152-160 and uses formulae optimal in certain Hardy spaces h(p,d);
see Sikorski,K., Optimal quadrature algorithms in HP spaces,
Num. Math., 39, 405-410 (1982).

Functions may have singularities at one or both end-points of the
interval (a,b).
}
\author{J.K. Lindsey}
\examples{
f <- function(x) sin(x)+cos(x)-x^2
int(f, a=0, b=2)
#
f <- function(x) exp(-(x-2)^2/2)/sqrt(2*pi)
int(f, a=0:3)
1-pnorm(0:3, 2)
}
