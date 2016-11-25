\name{Pareto}
\title{Pareto Distribution}
\usage{
dpareto(y, m, s, log=FALSE)
ppareto(q, m, s)
qpareto(p, m, s)
rpareto(n, m, s)
}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}
\description{
These functions provide information about the Pareto distribution
with location parameter equal to \code{m} and dispersion equal to
\code{s}: density, cumulative distribution, quantiles, log hazard, and
random generation.

The Pareto distribution has density
\deqn{
f(y) = \frac{\sigma }{\mu (\sigma-1)(1 + y/(\mu (\sigma-1)))^{\sigma+1}}}{
f(y) = s (1 + y/(m (s-1)))^(-s-1)/(m (s-1))}
where \eqn{\mu}{m} is the mean parameter of the distribution and
\eqn{\sigma}{s} is the dispersion.

This distribution can be obtained as a mixture distribution from the
exponential distribution using a gamma mixing distribution.
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
\code{\link{dexp}} for the exponential distribution.
}
\examples{
dpareto(5, 2, 2)
ppareto(5, 2, 2)
qpareto(0.9, 2, 2)
rpareto(10, 2, 2)
}
\keyword{distribution}
