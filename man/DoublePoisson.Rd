\name{Double Poisson}
\title{The Double Poisson Distribution}
\usage{
ddoublepois(y, m, s)
pdoublepois(q, m, s)
}
\alias{ddoublepois}
\alias{pdoublepois}
\arguments{
\item{y}{vector of counts}
\item{q}{vector of quantiles}
\item{m}{vector of means}
\item{s}{vector of overdispersion parameters}
}
\value{
These functions provide information about the double Poisson
distribution with parameters \code{m} and \code{s}. \code{ddoublepois}
gives the density and \code{pdoublepois} gives the distribution function.

The double Poisson distribution with \code{mu} \eqn{= m} has density
\deqn{p(y) = c({m}, {s}) {s}^({y}/{m}) ({m}/{y})^(y\log({s})) {y}^{(y-1)} / {y!}%
    }{p(y) = c(m,s) s^(y/m) (m/y)^(y log(s)) y^(y-1) / y!}
for \eqn{y = 0, \ldots}, where c(.) is a normalizing constant.
}

\seealso{
\code{\link{dpois}} for the Poisson, \code{\link{dmultpois}} for
the multiplicative Poisson, and \code{\link{dnbinom}}
for the negative binomial distribution.
}

\examples{
ddoublepois(5,5,0.9)
pdoublepois(5,5,0.9)
}
\keyword{distribution}
