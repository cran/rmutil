\name{Multiplicative Poisson}
\title{The Multiplicative Poisson Distribution}
\usage{
dmultpois(y, m, s)
pmultpois(y, m, s)
}
\alias{dmultpois}
\alias{pmultpois}
\arguments{
\item{y}{vector of counts}
\item{q}{vector of quantiles}
\item{m}{scalar or vector of means}
\item{s}{scalar or vector of overdispersion parameters}
}
\value{
These functions provide information about the multiplicative Poisson
distribution with parameters \code{m} and \code{s}. \code{dmultpois}
gives the density and \code{pmultpois} gives the distribution function.

The multiplicative Poisson distribution with \code{mu} \eqn{= m} has density
\deqn{p(y) = c({m}, {s})\exp({-m}) {m}^{y} {s}^({y}^2) / {y!}%
    }{p(y) = c(m,s) exp(-m) m^y s^(y^2) / y!}
for \eqn{y = 0, \ldots}, where c(.) is a normalizing constant.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dpois}} for the Poisson, \code{\link{ddoublepois}} for
the double Poisson, and \code{\link{dnbinom}}
for the negative binomial distribution.
}
\examples{
dmultpois(5,5,0.9)
pmultpois(5,5,0.9)
}
\keyword{distribution}
