\name{Generalized Weibull}
\title{The Generalized Weibull Distribution}
\usage{
dgweibull(y, m, s, f)
pgweibull(q, m, s, f)
hgweibull(y, m, s, f)
}
\alias{dgweibull}
\alias{pgweibull}
\alias{hgweibull}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of location parameters.}
  \item{s}{vector of dispersion parameters.}
  \item{f}{vector of family parameters.}
}
\description{
  These functions provide information about the generalized Weibull
distribution with location parameter equal to \code{m}, dispersion equal
to \code{s}, and family parameter equal to \code{f}.
  \code{dgweibull} gives the density, \code{pgweibull} gives the distribution
  function, \code{hgweibull} gives the hazard function.

  The generalized Weibull distribution has density
  \deqn{
    f(y) = \frac{\sigma \nu y^{\sigma-1} (1-\exp(-(y/\mu)^\sigma))^{\nu-1}
    \exp(-(y/\mu)^\sigma)}{\mu^\sigma}
    }{
    f(y) = s f y^(s-1) (1-exp(-(y/m)^s))^(f-1) exp(-(y/m)^s)/m^s
  }
	
  where \eqn{\mu}{m} is the location parameter of the distribution,
  \eqn{\sigma}{s} is the dispersion, and \eqn{\nu}{f} is the family
  parameter.

\eqn{\nu=\infty}{f=infinity} gives a Weibull distribution, for
\eqn{\sigma=1}{s=1}, \eqn{\nu<0}{f<0} a generalized F distribution,
and for \eqn{\sigma>0}{s>0}, \eqn{\nu<0}{f<0} a Burr type XII distribution.
}
\author{J.K. Lindsey}
\seealso{
  \code{\link{dweibull}} for the Weibull distribution,
\code{\link{df}} for the F distribution,
\code{\link{dburr}} for the Burr distribution.
}
\examples{
dgweibull(5, 1, 3, 2)
pgweibull(5, 1, 3, 2)
}
\keyword{distribution}
