\name{Burr}
\title{The Burr Distribution}
\usage{
dburr(y, m, s, f)
pburr(q, m, s, f)
hburr(y, m, s, f)
}
\alias{dburr}
\alias{pburr}
\alias{hburr}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of location parameters.}
  \item{s}{vector of dispersion parameters.}
  \item{f}{vector of family parameters.}
}
\description{
  These functions provide information about the Burr distribution with
location parameter equal to \code{m}, dispersion equal to
\code{s}, and family parameter equal to \code{f}.
  \code{dburr} gives the density, \code{pburr} gives the distribution
  function, \code{hburr} gives the hazard function.

  The Burr distribution has density
  \deqn{
    f(y) = \frac{\sigma (y / \mu)^{\sigma-1}}
    {\mu^\sigma (1+(y/\mu)^\sigma/\nu)^{\nu+1}}}{
    f(y) = s (y/m)^(s-1)/(m^s (1+(y/m)^s/f)^(f+1))
  }
  where \eqn{\mu}{m} is the location parameter of the distribution,
  \eqn{\sigma}{s} is the dispersion, and \eqn{\nu}{f} is the family
  parameter.
}
\author{J.K. Lindsey}
\examples{
dburr(2, 5, 1, 2)
pburr(2, 5, 1, 2)
}
\keyword{distribution}
