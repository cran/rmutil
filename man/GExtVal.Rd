\name{Generalized Extreme Value}
\title{The Generalized Extreme Value Distribution}
\usage{
dgextval(y, m, s, f)
pgextval(q, m, s, f)
hgextval(y, m, s, f)
}
\alias{dgextval}
\alias{pgextval}
\alias{hgextval}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of location parameters.}
  \item{s}{vector of dispersion parameters.}
  \item{f}{vector of family parameters.}
}
\description{
  These functions provide information about the generalized extreme
value distribution with location parameter equal to \code{m}, dispersion
equal to \code{s}, and family parameter equal to \code{f}.
  \code{dgextval} gives the density, \code{pgextval} gives the distribution
  function, \code{hgextval} gives the hazard function.

  The generalized extreme value distribution has density
  \deqn{
    f(y) =
    \nu y^{\nu-1} \exp(y^\nu) \frac{\sigma}{\mu}
    \frac{\exp(y^\nu)}{\mu^{\sigma-1}/(1-I(\nu>0)+sign(\nu) exp(-\mu^-\sigma))}
    \exp(-(\exp(y^\nu)/\mu)^\sigma)}{
    f(y) = f y^(f-1) exp(y^f) (s/m) (exp(y^f)/m)^(s-1)
  exp(-(exp(y^f)/m)^s)/(1-I(f>0)+sign(f) exp(-m^-s))
  }
	
  where \eqn{\mu}{m} is the location parameter of the distribution,
  \eqn{\sigma}{s} is the dispersion, \eqn{\nu}{f} is the family
  parameter, \eqn{I()} is the indicator function, and \eqn{y>0}.

\eqn{\nu=0}{f=0} yields a Weibull distribution and
  \eqn{\nu=1}{f=1} a truncated extreme value distribution.
\author{J.K. Lindsey}
}
\seealso{
  \code{\link{dweibull}} for the Weibull distribution.
}
\examples{
dgextval(1, 2, 1, 2)
pgextval(1, 2, 1, 2)
}
\keyword{distribution}
