\name{Double Binomial}
\title{The Double Binomial Distribution}
\usage{
ddoublebinom(y, n, m, s)
pdoublebinom(q, n, m, s)
}
\alias{ddoublebinom}
\alias{pdoublebinom}
\arguments{
\item{y}{vector of frequencies}
\item{q}{vector of quantiles}
\item{n}{vector of totals}
\item{m}{vector of probabilities}
\item{s}{vector of overdispersion parameters}
}
\value{
These functions provide information about the double binomial
distribution with parameters \code{m} and \code{s}. \code{ddoublebinom}
gives the density and \code{pdoublebinom} gives the distribution function.

The double binomial distribution with total \eqn{= n} and
\code{prob} \eqn{= m} has density
\deqn{p(y) = c({n}, {m}, {s}){n \choose y} ({m}/{y})^({y}\log({s})) {((1-m)/(n-y))}^(({n-y})\log({s})) {y}^{y} {(n-y)}^{(n-y)})%
    }{p(y) = c(n,m,s) Choose(n,y) n^(n log(s)) (m/y)^(y log(s))
    ((1-m)/(n-y))^(s(n-y)) y^y (n-y)^(n-y)}
for \eqn{y = 0, \ldots, n}, where c(.) is a normalizing constant.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dbinom}} for the binomial, \code{\link{dmultbinom}} for
the multiplicative binomial, and \code{\link{dbetabinom}}
for the beta binomial distribution.
}
\examples{
# Compute P(45 < y < 55) for y double binomial(100,0.5,1.1)
sum(ddoublebinom(46:54, 100, 0.5, 1.1))
pdoublebinom(54, 100, 0.5, 1.1)-pdoublebinom(45, 100, 0.5, 1.1)
}
\keyword{distribution}
