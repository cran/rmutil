\name{Multiplicative Binomial}
\title{Multiplicative Binomial Distribution}
\usage{
dmultbinom(y, size, m, s, log=FALSE)
pmultbinom(q, size, m, s)
qmultbinom(p, size, m, s)
rmultbinom(n, size, m, s)
}
\alias{dmultbinom}
\alias{pmultbinom}
\alias{qmultbinom}
\alias{rmultbinom}
\description{
These functions provide information about the multiplicative binomial
distribution with parameters \code{m} and \code{s}: density,
cumulative distribution, quantiles, and random generation.

The multiplicative binomial distribution with total \eqn{= n} and
\code{prob} \eqn{= m} has density
\deqn{p(y) = c({n}, {m}, {s}){n \choose y} {m}^{y} {(1-m)}^{n-y} {s}^{(y(n-y))}
}{p(y) = c(n,m,s) Choose(n,y) m^y (1-m)^(n-y) s^(y(n-y))}
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
\code{\link{dbinom}} for the binomial, \code{\link[rmutil]{ddoublebinom}} for
the double binomial, and \code{\link[rmutil]{dbetabinom}}
for the beta binomial distribution.
}
\examples{
# compute P(45 < y < 55) for y multiplicative binomial(100,0.5,1.1)
sum(dmultbinom(46:54, 100, 0.5, 1.1))
pmultbinom(54, 100, 0.5, 1.1)-pmultbinom(45, 100, 0.5, 1.1)
pmultbinom(2,10,0.5,1.1)
qmultbinom(0.025,10,0.5,1.1)
rmultbinom(10,10,0.5,1.1)
}
\keyword{distribution}
