\name{Inverse Gaussian}
\title{The Inverse Gaussian Distribution}
\usage{
dinvgauss(y, m, s)
pinvgauss(q, m, s)
hinvgauss(y m, s)
}
\alias{dinvgauss}
\alias{pinvgauss}
\alias{hinvgauss}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of means.}
  \item{s}{vector of dispersion parameters.}
}
\description{
  These functions provide information about the inverse Gaussian
distribution with mean equal to \code{m} and dispersion equal to \code{s}.
  \code{dinvgauss} gives the density, \code{pinvgauss} gives the distribution
  function, \code{hinvgauss} gives the hazard function.

  The inverse Gaussian distribution has density
  \deqn{
    f(y) =
    \frac{1}{\sqrt{2\pi\sigma y^3}} e^{-(y-\mu)^2/(2 y \sigma^2 m^2)}}{
    f(y) = 1/sqrt(2 pi s y^3) e^-((y - m)^2/(2 y s^2 m^2))}
  where \eqn{\mu}{m} is the mean of the distribution and
  \eqn{\sigma}{s} is the dispersion.
}
\author{J.K. Lindsey}
\seealso{
  \code{\link{dnorm}} for the normal distribution and
  \code{\link{dlnorm}} for the \emph{Log}normal distribution.
}
\examples{
dinvgauss(5, 5, 1)
pinvgauss(5, 5, 1)
}
\keyword{distribution}
