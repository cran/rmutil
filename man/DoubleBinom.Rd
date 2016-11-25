\name{Double Binomial}
\title{Double Binomial Distribution}
\usage{
ddoublebinom(y, size, m, s, log=FALSE)
pdoublebinom(q, size, m, s)
qdoublebinom(p, size, m, s)
rdoublebinom(n, size, m, s)
}
\alias{ddoublebinom}
\alias{pdoublebinom}
\alias{qdoublebinom}
\alias{rdoublebinom}
\description{
These functions provide information about the double binomial
distribution with parameters \code{m} and \code{s}: density,
cumulative distribution, quantiles, and random generation.

The double binomial distribution with total \eqn{= n} and
\code{prob} \eqn{= m} has density
\deqn{p(y) = c({n}, {m}, {s}){n \choose y} {n}^{{n}{s}} ({m}/{y})^({y}{s}) {((1-m)/(n-y))}^(({n-y})s y) {y}^{y} {(n-y)}^{(n-y)})%
}{p(y) = c(n,m,s) Choose(n,y) n^(n s) (m/y)^(y s)
((1-m)/(n-y))^(s(n-y)) y^y (n-y)^(n-y)}
for \eqn{y = 0, \ldots, n}, where c(.) is a normalizing constant.
}
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
\code{\link{dbinom}} for the binomial, \code{\link[rmutil]{dmultbinom}} for
the multiplicative binomial, and \code{\link[rmutil]{dbetabinom}}
for the beta binomial distribution.
}
\examples{
# compute P(45 < y < 55) for y double binomial(100,0.5,1.1)
sum(ddoublebinom(46:54, 100, 0.5, 1.1))
pdoublebinom(54, 100, 0.5, 1.1)-pdoublebinom(45, 100, 0.5, 1.1)
pdoublebinom(2,10,0.5,1.1)
qdoublebinom(0.05,10,0.5,1.1)
rdoublebinom(10,10,0.5,1.1)
}
\keyword{distribution}
