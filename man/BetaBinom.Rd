\name{Beta Binomial}
\title{The Beta Binomial Distribution}
\usage{
dbetabinom(y, n, m, s)
pbetabinom(q, n, m, s)
}
\alias{dbetabinom}
\alias{pbetabinom}
\arguments{
\item{y}{vector of frequencies}
\item{q}{vector of quantiles}
\item{n}{vector of totals}
\item{m}{vector of probabilities}
\item{s}{vector of overdispersion parameters}
}
\value{
These functions provide information about the beta binomial
distribution with parameters \code{m} and \code{s}. \code{dbetabinom}
gives the density and \code{pbetabinom} gives the distribution function.

The beta binomial distribution with total \eqn{= n} and
\code{prob} \eqn{= m} has density
\deqn{p(y) = \frac{B(y+\sigma \mu, n-y+\sigma*(1-\mu)) {n \choose y}
}{ B(s m,s (1-m))}%
    }{p(y) = B(y+s m,n-y+s (1-m)) Choose(n,y) / B(s m,s (1-m))}
for \eqn{y = 0, \ldots, n} where \eqn{B()} is the beta function.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dbinom}} for the binomial, \code{\link{ddoublebinom}} for
the double binomial, and \code{\link{dmultbinom}}
for the multiplicative binomial distribution.
}
\examples{
# Compute P(45 < y < 55) for y beta binomial(100,0.5,1.1)
sum(dbetabinom(46:54, 100, 0.5, 1.1))
pbetabinom(54,100,0.5,1.1)-pbetabinom(45,100,0.5,1.1)
}
\keyword{distribution}
