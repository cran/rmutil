\name{Beta Binomial}
\title{Beta Binomial Distribution}
\usage{
dbetabinom(y, size, m, s, log=FALSE)
pbetabinom(q, size, m, s)
qbetabinom(p, size, m, s)
rbetabinom(n, size, m, s)
}
\alias{dbetabinom}
\alias{pbetabinom}
\alias{qbetabinom}
\alias{rbetabinom}
\description{
These functions provide information about the beta binomial
distribution with parameters \code{m} and \code{s}: density,
cumulative distribution, quantiles, and random generation.

The beta binomial distribution with total \eqn{= n} and
\code{prob} \eqn{= m} has density
\deqn{p(y) = \frac{B(y+\sigma \mu, n-y+\sigma*(1-\mu)) {n \choose y}
}{B(s m,s (1-m))}%
}{p(y) = B(y+s m,n-y+s (1-m)) Choose(n,y) / B(s m,s (1-m))}
for \eqn{y = 0, \ldots, n} where \eqn{B()} is the beta function.}
\arguments{
\item{y}{vector of frequencies}
\item{q}{vector of quantiles}
\item{p}{vector of probabilities}
\item{n}{number of values to generate}
\item{size}{vector of totals}
\item{m}{vector of probabilities of success}
\item{s}{vector of overdispersion parameters}
\item{log}{if TRUE, log probabilities are supplied.}
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dbinom}} for the binomial, \code{\link[rmutil]{ddoublebinom}} for
the double binomial, and \code{\link[rmutil]{dmultbinom}}
for the multiplicative binomial distribution.
}
\examples{
# compute P(45 < y < 55) for y beta binomial(100,0.5,1.1)
sum(dbetabinom(46:54, 100, 0.5, 1.1))
pbetabinom(54,100,0.5,1.1)-pbetabinom(45,100,0.5,1.1)
pbetabinom(2,10,0.5,1.1)
qbetabinom(0.33,10,0.5,1.1)
rbetabinom(10,10,0.5,1.1)
}
\keyword{distribution}
