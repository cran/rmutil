\name{Multiplicative Binomial}
\title{The Multiplicative Binomial Distribution}
\usage{
dmultbinom(y, n, m, s)
pmultbinom(q, n, m, s)
}
\alias{dmultbinom}
\alias{pmultbinom}
\arguments{
\item{y}{vector of frequencies}
\item{q}{vector of quantiles}
\item{n}{vector of totals}
\item{m}{vector of probabilities}
\item{s}{vector of overdispersion parameters}
}
\value{
These functions provide information about the multiplicative binomial
distribution with parameters \code{m} and \code{s}. \code{dmultbinom}
gives the density and \code{pmultbinom} gives the distribution function.

The multiplicative binomial distribution with total \eqn{= n} and
\code{prob} \eqn{= m} has density
\deqn{p(y) = c({n}, {m}, {s}){n \choose y} {m}^{y} {(1-m)}^{n-y} {s}^{(y(n-y))}%
    }{p(y) = c(n,m,s) Choose(n,y) m^y (1-m)^(n-y) s^(y(n-y))}
for \eqn{y = 0, \ldots, n}, where c(.) is a normalizing constant.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dbinom}} for the binomial, \code{\link{ddoublebinom}} for
the double binomial, and \code{\link{dbetabinom}}
for the beta binomial distribution.
}
\examples{
# Compute P(45 < y < 55) for y multiplicative binomial(100,0.5,1.1)
sum(dmultbinom(46:54, 100, 0.5, 1.1))
pmultbinom(54, 100, 0.5, 1.1)-pmultbinom(45, 100, 0.5, 1.1)
}
\keyword{distribution}
