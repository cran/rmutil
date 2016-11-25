\name{Levy}
\title{Levy Distribution}
\usage{
dlevy(y, m=0, s=1, log=FALSE)
plevy(q, m=0, s=1)
qlevy(p, m=0, s=1)
rlevy(n, m=0, s=1)
}
\alias{dlevy}
\alias{plevy}
\alias{qlevy}
\alias{rlevy}
\description{
These functions provide information about the Levy distribution
with location parameter equal to \code{m} and dispersion equal to
\code{s}: density, cumulative distribution, quantiles, and
random generation.

The Levy distribution has density
\deqn{
f(y) = \sqrt{\frac{\sigma}{2 \pi (y-\mu)^3}} \exp(-\sigma/(2 (y-\mu)))}{
f(y) = sqrt(s/(2 pi (y-m)^3)) exp(-s/(2 (y-m)))}
where \eqn{\mu}{m} is the location parameter of the distribution and
\eqn{\sigma}{s} is the dispersion, and \eqn{y>\mu}{y>m}.
}
\arguments{
\item{y}{vector of responses.}
\item{q}{vector of quantiles.}
\item{p}{vector of probabilities}
\item{n}{number of values to generate}
\item{m}{vector of location parameters.}
\item{s}{vector of dispersion parameters.}
\item{log}{if TRUE, log probabilities are supplied.}
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
qlevy(0.6, 2, 1)
rlevy(10, 2, 1)
}
\keyword{distribution}
