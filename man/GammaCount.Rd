\name{Gamma Count}
\title{Gamma Count Distribution}
\usage{
dgammacount(y, m, s, log=FALSE)
pgammacount(q, m, s)
qgammacount(p, m, s)
rgammacount(n, m, s)
}
\alias{dgammacount}
\alias{pgammacount}
\alias{qgammacount}
\alias{rgammacount}
\description{
These functions provide information about the gamma count
distribution with parameters \code{m} and \code{s}: density,
cumulative distribution, quantiles, and random generation.

The gamma count distribution with \code{prob} \eqn{= m} has density
\deqn{p(y) = pgamma(\mu \sigma,y \sigma,1)-pgamma(\mu \sigma,(y+1) \sigma,1)
}{p(y) = pgamma(m s, y s, 1)-pgamma(m s, (y+1) s, 1)}
for \eqn{y = 0, \ldots, n} where
\eqn{pgamma(\mu \sigma,0,1)=1}{pgamma(m s, 0, 1)=1}.
}
\arguments{
\item{y}{vector of frequencies}
\item{q}{vector of quantiles}
\item{p}{vector of probabilities}
\item{n}{number of values to generate}
\item{m}{vector of probabilities}
\item{s}{vector of overdispersion parameters}
\item{log}{if TRUE, log probabilities are supplied.}
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dpois}} for the Poisson, \code{\link[rmutil]{dconsul}} for
the Consul generalized Poisson, \code{\link[rmutil]{ddoublepois}} for
the double Poisson, \code{\link[rmutil]{dmultpois}}
for the multiplicative Poisson distributions, and \code{\link{dnbinom}}
for the negative binomial distribution.
}
\examples{
dgammacount(5,10,0.9)
pgammacount(5,10,0.9)
qgammacount(0.08,10,0.9)
rgammacount(10,10,0.9)
}
\keyword{distribution}
