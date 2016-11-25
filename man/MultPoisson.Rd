\name{MultPoisson}
\title{Multiplicative Poisson Distribution}
\usage{
dmultpois(y, m, s, log=FALSE)
pmultpois(q, m, s)
qmultpois(p, m, s)
rmultpois(n, m, s)
}
\alias{dmultpois}
\alias{pmultpois}
\alias{qmultpois}
\alias{rmultpois}
\description{
These functions provide information about the multiplicative Poisson
distribution with parameters \code{m} and \code{s}: density,
cumulative distribution, quantiles, and random generation.

The multiplicative Poisson distribution with \code{mu} \eqn{= m} has density
\deqn{p(y) = c({\mu}, {\lambda})\exp({-\mu}) {\mu}^{y} {\lambda}^({y}^2) / {y!}%
}{p(y) = c(m,s) exp(-m) m^y s^(y^2) / y!}
with \eqn{s <= 1} for \eqn{y = 0, \ldots}, where c(.) is a normalizing
constant.

Note that it only allows for underdispersion, not being defined for
\eqn{s > 1}.
}
\arguments{
\item{y}{vector of counts}
\item{q}{vector of quantiles}
\item{p}{vector of probabilities}
\item{n}{number of values to generate}
\item{m}{scalar or vector of means}
\item{s}{scalar or vector of overdispersion parameters, all of which
must lie in (0,1). }
\item{log}{if TRUE, log probabilities are supplied.}
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dpois}} for the Poisson, \code{\link[rmutil]{ddoublepois}} for
the double Poisson, \code{\link[rmutil]{dpvfpois}} for the power variance
function Poisson, \code{\link[rmutil]{dconsul}} for the Consul
generalized Poisson, \code{\link[rmutil]{dgammacount}} for the gamma count, and
\code{\link{dnbinom}} for the negative binomial distribution.
}
\examples{
dmultpois(5,10,0.9)
pmultpois(5,10,0.9)
qmultpois(0.85,10,0.9)
rmultpois(10,10,0.9)
}
\keyword{distribution}
