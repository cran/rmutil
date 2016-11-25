\name{Consul}
\title{Consul Distribution}
\usage{
dconsul(y, m, s, log=FALSE)
pconsul(q, m, s)
qconsul(p, m, s)
rconsul(n, m, s)
}
\alias{dconsul}
\alias{pconsul}
\alias{qconsul}
\alias{rconsul}
\description{
These functions provide information about the Consul
distribution with parameters \code{m} and \code{s}: density,
cumulative distribution, quantiles, and random generation.

The Consul distribution with \code{mu} \eqn{= m} has density
\deqn{p(y) = \mu \exp(-(\mu+y(\lambda-1))/\lambda) (\mu+y(\lambda-1))^(y-1)/(\lambda^y y!)%
}{p(y) = m exp(-(m+y(s-1))/s) (m+y(s-1))^(y-1)/(s^y y!)}
for \eqn{y = 0, \ldots}.
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
\code{\link{dpois}} for the Poisson, \code{\link[rmutil]{ddoublepois}} for
the double Poisson, \code{\link[rmutil]{dmultpois}} for
the multiplicative Poisson, and \code{\link[rmutil]{dpvfpois}} for the power
variance function Poisson.
}

\examples{
dconsul(5,10,0.9)
pconsul(5,10,0.9)
qconsul(0.08,10,0.9)
rconsul(10,10,0.9)
}
\keyword{distribution}
