\name{Generalized Gamma}
\title{The Generalized Gamma Distribution}
\usage{
dggamma(y, m, s, f)
pggamma(q, m, s, f)
hggamma(y, m, s, f)
}
\alias{dggamma}
\alias{pggamma}
\alias{hggamma}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of location parameters.}
  \item{s}{vector of dispersion parameters.}
  \item{f}{vector of family parameters.}
}
\description{
  These functions provide information about the generalized gamma
distribution with location parameter equal to \code{m}, dispersion equal
to \code{s}, and family parameter equal to \code{f}.
  \code{dggamma} gives the density, \code{pggamma} gives the distribution
  function, \code{hggamma} gives the hazard function.

  The generalized gamma distribution has density
  \deqn{
    f(y) =
    \frac{1}{(\mu/\sigma)^{\nu\sigma} Gamma(\sigma)} y^{\nu(\sigma-1)}
    \exp(-(y \sigma/\mu)^\nu)}{
    f(y) = 1/((m/s)^(fs) Gamma(s)) y^(f(s-1)) exp(-(y s/m)^f)
  }
	
  where \eqn{\mu}{m} is the location parameter of the distribution,
  \eqn{\sigma}{s} is the dispersion, and \eqn{\nu}{f} is the family
  parameter.

\eqn{\nu=1}{f=1} yields a gamma distribution, \eqn{\sigma=1}{s=1} a
Weibull distribution, and \eqn{\sigma=\infty}{s=infinity} a
log normal distribution.
}
\author{J.K. Lindsey}
\seealso{
  \code{\link{dgamma}} for the gamma distribution,
\code{\link{dweibull}} for the Weibull distribution, \code{\link{dlnorm}}
for the log normal distribution.
}
\examples{
dggamma(2, 5, 4, 2)
pggamma(2, 5, 4, 2)
}
\keyword{distribution}
