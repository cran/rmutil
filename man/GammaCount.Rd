\name{Gamma Count}
\title{The Gamma Count Distribution}
\usage{
dgammacount(y, m, s)
pgammacount(q, m, s)
}
\alias{dgammacount}
\alias{pgammacount}
\arguments{
\item{y}{vector of frequencies}
\item{q}{vector of quantiles}
\item{m}{vector of probabilities}
\item{s}{vector of overdispersion parameters}
}
\value{
These functions provide information about the gamma count
distribution with parameters \code{m} and \code{s}. \code{dgammacount}
gives the density and \code{pgammacount} gives the distribution function.

The gamma count distribution with \code{prob} \eqn{= m} has density
\deqn{p(y) = pgamma(\mu \sigma,y \sigma,1)-pgamma(\mu \sigma,(y+1) \sigma,1)
    }{p(y) = pgamma(m s, y s, 1)-pgamma(m s, (y+1) s, 1)}
for \eqn{y = 0, \ldots, n} where
\eqn{pgamma(\mu \sigma,0,1)=1}{pgamma(m s, 0, 1)=1}.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dpois}} for the Poisson, \code{\link{ddoublepois}} for
the double Poisson, and \code{\link{dmultpois}}
for the multiplicative Poisson distributions.
}
\examples{
dgammacount(5,5,0.9)
pgammacount(5,5,0.9)
}
\keyword{distribution}
