\name{Levy}
\title{The Levy Distribution}
\usage{
dlevy(y, m, s)
plevy(q, m, s)
}
\alias{dlevy}
\alias{plevy}
\arguments{
  \item{y}{vector of responses.}
  \item{q}{vector of quantiles.}
  \item{m}{vector of location parameters.}
  \item{s}{vector of dispersion parameters.}
}
\description{
  These functions provide information about the Levy distribution
with location parameter equal to \code{m} and dispersion equal to
\code{s}.  \code{dlevy} gives the density, \code{plevy} gives
the distribution function.

  The Levy distribution has density
  \deqn{
    f(y) = \sqrt{\frac{\sigma}{2 \pi (y-\mu)^3}} \exp(-\sigma/(2 (y-\mu)))}{
    f(y) = sqrt(s/(2 pi (y-m)^3)) exp(-s/(2 (y-m)))
  }
  where \eqn{\mu}{m} is the location parameter of the distribution and
  \eqn{\sigma}{s} is the dispersion, and \eqn{y>\mu}{y>m}.
}
\author{J.K. Lindsey}
\seealso{
  \code{\link{dnorm}} for the normal distribution and
  \code{\link{dcauchy}} for the Cauchy distribution, two other stable
  distributions.
}
\examples{
dlevy(5, 2, 1)
plevy(5, 2, 1)
}
\keyword{distribution}
