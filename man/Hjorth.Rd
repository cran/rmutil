\name{Hjorth}
\title{The Hjorth Distribution}
\usage{
dhjorth(y, m, s, f)
phjorth(q, m, s, f)
hhjorth(y, m, s, f)
}
\alias{dhjorth}
\alias{phjorth}
\alias{hhjorth}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of location parameters.}
  \item{s}{vector of dispersion parameters.}
  \item{f}{vector of family parameters.}
}
\description{
  These functions provide information about the Hjorth
distribution with location parameter equal to \code{m}, dispersion equal
to \code{s}, and family parameter equal to \code{f}.
  \code{dhjorth} gives the density, \code{phjorth} gives the distribution
  function, \code{hhjorth} gives the hazard function.

  The Hjorth distribution has density
  \deqn{
    f(y) = (1+\sigma y)^{-\nu/\sigma} \exp(-(y/\mu)^2/2)
    (\frac{y}{\mu^2}+\frac{\nu}{1+\sigma y})
    }{
    f(y) = (1+s y)^(-f/s) exp(-(y/m)^2/2) (y/m^2+f/(1+s y))
  }
	
  where \eqn{\mu}{m} is the location parameter of the distribution,
  \eqn{\sigma}{s} is the dispersion, and \eqn{\nu}{f} is the family
  parameter.
}
\author{J.K. Lindsey}
\examples{
dhjorth(5, 5, 5, 2)
phjorth(5, 5, 5, 2)
}
\keyword{distribution}
