\name{Generalized Inverse Gaussian}
\title{The Generalized Inverse Gaussian Distribution}
\usage{
dginvgauss(y, m, s, f)
pginvgauss(q, m, s, f)
}
\alias{dginvgauss}
\alias{pginvgauss}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of means.}
  \item{s}{vector of dispersion parameters.}
  \item{f}{vector of family parameters.}
}
\description{
  These functions provide information about the generalized inverse
Gaussian distribution with mean equal to \code{m}, dispersion equal to
\code{s}, and family parameter equal to \code{f}.
  \code{dginvgauss} gives the density, \code{pginvgauss} gives the distribution
  function.

  The generalized inverse Gaussian distribution has density
  \deqn{
    f(y) =
    \frac{y^{\nu-1}}{2 \mu^\nu K(1/(\sigma \mu),abs(\nu))}
	\exp(-(y-\mu)^2/(2 y \sigma^2 \mu^2))}{
    f(y) = y^(f-1)/(2 m^f K(1/(s m),abs(f)))
	 exp(-((y - m)^2/(2 y s^2 m^2)))
  }
  where \eqn{\mu}{m} is the mean of the distribution,
  \eqn{\sigma}{s} the dispersion, \eqn{\nu}{f} is the family
  parameter, and \eqn{K()}{K()} is the fractional Bessel function of
  the third kind.

\eqn{\nu=-1/2}{f=-1/2} yields an inverse Gaussian distribution,
\eqn{\sigma=\infty}{s=infinity}, \eqn{\nu>0}{f>0} a gamma
distribution, and \eqn{\nu=0}{f=0} a hyperbola distribution.
}
\author{J.K. Lindsey}
\seealso{
  \code{\link{dinvgauss}} for the inverse Gaussian distribution.
}
\examples{
dginvgauss(10, 3, 1, 1)
pginvgauss(10, 3, 1, 1)
}
\keyword{distribution}
