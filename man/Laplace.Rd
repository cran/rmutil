\name{Laplace}
\title{The Laplace Distribution}
\usage{
dlaplace(y, m, s)
plaplace(q, m, s)
hlaplace(y, m, s)
}
\alias{dlaplace}
\alias{plaplace}
\alias{hlaplace}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of location parameters.}
  \item{s}{vector of dispersion parameters.}
}
\description{
  These functions provide information about the Laplace distribution
with location parameter equal to \code{m} and dispersion equal to
\code{s}.  \code{dlaplace} gives the density, \code{plaplace} gives
the distribution function, \code{hlaplace} gives the hazard function.

  The Laplace distribution has density
  \deqn{
    f(y) = \frac{\exp(-\abs(y-\mu)/\sigma)}{(2*\sigma)}}{
    f(y) = exp(-abs(y-m)/s)/(2*s)
  }
  where \eqn{\mu}{m} is the location parameter of the distribution and
  \eqn{\sigma}{s} is the dispersion.
}
\author{J.K. Lindsey}
\seealso{
  \code{\link{dexp}} for the exponential distribution and
  \code{\link{dcauchy}} for the Cauchy distribution.
}
\examples{
dlaplace(5, 2, 1)
plaplace(5, 2, 1)
}
\keyword{distribution}
