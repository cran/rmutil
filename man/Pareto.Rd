\name{Pareto}
\title{The Pareto Distribution}
\usage{
dpareto(y, m, s)
ppareto(q, m, s)
hpareto(y, m, s)
}
\alias{dpareto}
\alias{ppareto}
\alias{hpareto}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of location parameters.}
  \item{s}{vector of dispersion parameters.}
}
\description{
  These functions provide information about the Pareto distribution
with location parameter equal to \code{m} and dispersion equal to
\code{s}.  \code{dpareto} gives the density, \code{ppareto} gives
the distribution function, \code{hpareto} gives the hazard function.

  The Pareto distribution has density
\deqn{
  f(y) = \frac{\sigma }{\mu (\sigma-1)(1 + y/(\mu (\sigma-1)))^{\sigma+1}}}{
  f(y) = s (1 + y/(m (s-1)))^(-s-1)/(m (s-1))
}
where \eqn{\mu}{m} is the mean parameter of the distribution and
\eqn{\sigma}{s} is the dispersion.

This distribution can be obtained as a mixture distribution from the
exponential distribution using a gamma mixing distribution.
}
\author{J.K. Lindsey}
\seealso{
  \code{\link{dexp}} for the exponential distribution.
}
\examples{
dpareto(5, 2, 1)
ppareto(5, 2, 1)
}
\keyword{distribution}
