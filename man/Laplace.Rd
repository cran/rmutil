\name{Laplace}
\title{Laplace Distribution}
\usage{
dlaplace(y, m=0, s=1, log=FALSE)
plaplace(q, m=0, s=1)
qlaplace(p, m=0, s=1)
rlaplace(n=1, m=0, s=1)
}
\alias{dlaplace}
\alias{plaplace}
\alias{qlaplace}
\alias{rlaplace}
\description{
These functions provide information about the Laplace distribution
with location parameter equal to \code{m} and dispersion equal to
\code{s}: density, cumulative distribution, quantiles, log hazard, and
random generation.

The Laplace distribution has density
\deqn{
f(y) = \frac{\exp(-abs(y-\mu)/\sigma)}{(2\sigma)}}{
f(y) = exp(-abs(y-m)/s)/(2*s)}
where \eqn{\mu}{m} is the location parameter of the distribution and
\eqn{\sigma}{s} is the dispersion.
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
\code{\link{dexp}} for the exponential distribution and
\code{\link{dcauchy}} for the Cauchy distribution.
}
\examples{
dlaplace(5, 2, 1)
plaplace(5, 2, 1)
qlaplace(0.95, 2, 1)
rlaplace(10, 2, 1)
}
\keyword{distribution}
