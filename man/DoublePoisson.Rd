\name{DoublePoisson}
\title{Double Poisson Distribution}
\usage{
ddoublepois(y, m, s, log=FALSE)
pdoublepois(q, m, s)
qdoublepois(p, m, s)
rdoublepois(n, m, s)
}
\alias{ddoublepois}
\alias{pdoublepois}
\alias{qdoublepois}
\alias{rdoublepois}
\description{
These functions provide information about the double Poisson
distribution with parameters \code{m} and \code{s}: density,
cumulative distribution, quantiles, and random generation.

The double Poisson distribution with \code{mu} \eqn{= m} has density
\deqn{p(y) = c({\mu}, {\lambda}) {\lambda}^({y}/{\mu}) ({\mu}/{y})^(y\log({\lambda})) {y}^{(y-1)} / {y!}%
}{p(y) = c(m,s) s^(y/m) (m/y)^(y log(s)) y^(y-1) / y!}
for \eqn{y = 0, \ldots}, where c(.) is a normalizing constant.
}
\arguments{
\item{y}{vector of counts}
\item{q}{vector of quantiles}
\item{p}{vector of probabilities}
\item{n}{number of values to generate}
\item{m}{vector of means}
\item{s}{vector of overdispersion parameters}
\item{log}{if TRUE, log probabilities are supplied.}
}
\seealso{
\code{\link{dpois}} for the Poisson, \code{\link[rmutil]{dconsul}} for
the Consul generalized Poisson, \code{\link[rmutil]{dgammacount}} for
the gamma count, \code{\link[rmutil]{dmultpois}} for the
multiplicative Poisson, \code{\link[rmutil]{dpvfpois}} for the power
variance function Poisson, and \code{\link{dnbinom}} for the negative
binomial distribution.
}

\examples{
ddoublepois(5,10,0.9)
pdoublepois(5,10,0.9)
qdoublepois(0.08,10,0.9)
rdoublepois(10,10,0.9)
}
\keyword{distribution}
