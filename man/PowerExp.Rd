\name{power exponential}
\title{The power exponential Distribution}
\usage{
dpowexp(y, m=0, s=1, f=1)
ppowexp(q, m=0, s=1, f=1)
}
\alias{dpowexp}
\alias{ppowexp}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of means.}
  \item{s}{vector of dispersion parameters.}
  \item{f}{vector of family parameters.}
}
\description{
  These functions provide information about the power exponential
distribution with mean parameter equal to \code{m}, dispersion equal
to \code{s}, and family parameter equal to \code{f}.
  \code{dpowexp} gives the density, \code{ppowexp} gives the distribution
  function.

  The power exponential distribution has density
  \deqn{
    f(y) = \frac{\exp(-(\abs{y-\mu}/\sqrt{\sigma})^{2 \nu}/2)}{
    \sqrt{\sigma} Gamma(1+1/(2 \nu)) 2^{1+1/(2 \nu)}}
    }{
    f(y) = exp(-(abs(y-m)/sqrt(s))^(2 f)/2)/
    (sqrt(s) Gamma(1+1/(2 f)) 2^(1+1/(2 f)))
  }
	
  where \eqn{\mu}{m} is the mean of the distribution,
  \eqn{\sigma}{s} is the dispersion, and \eqn{\nu}{f} is the family
  parameter. \eqn{\nu=1}{f=1} yields a normal distribution,
  \eqn{\nu=0.5}{f=0.5} a Laplace distribution, and
  \eqn{\nu=\infty}{f=Inf}  a uniform distribution.
}
\author{J.K. Lindsey}
\examples{
dpowexp(5, 5, 1, 2)
ppowexp(5, 5, 1, 2)
}
\keyword{distribution}
