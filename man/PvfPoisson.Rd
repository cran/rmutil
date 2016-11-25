\name{PvfPoisson}
\title{Power Variance Function Poisson Distribution}
\usage{
dpvfpois(y, m, s, f, log=FALSE)
ppvfpois(q, m, s, f)
qpvfpois(p, m, s, f)
rpvfpois(n, m, s, f)
}
\alias{dpvfpois}
\alias{ppvfpois}
\alias{qpvfpois}
\alias{rpvfpois}
\description{
These functions provide information about the overdispersed power
variance function Poisson distribution with parameters \code{m},
\code{s}, and \code{f}: density, cumulative distribution, quantiles,
and random generation. This function is obtained from a Poisson
distribution as a mixture with a power variance distribution. In the
limit, for \code{f=0}, the mixing distribution is gamma so that it is
a negative binomial distribution.  For \code{f=0.5}, the mixing
distribution is inverse Gaussian. For \code{f<0}, the mixing
distribution is a compound distribution of the sum of a Poisson number
of gamma distributions. For \code{f=1}, it is undefined.

The power variance function Poisson distribution with \code{m} \eqn{=
\mu}, the mean, \code{s} \eqn{= \theta}, and \code{f} \eqn{= \alpha}
has density
\deqn{p(y) =
{\exp(-\mu((\theta+1)^\alpha/\theta^\alpha-\theta)/\alpha)\over y!}
\sum_{i=1}^y c_{yi}(\alpha)\mu^i(\theta+1)^{i\alpha-y}/\theta^{i(\alpha-1)}
}{p(y) = (exp(-m((s+1)^f/s^f-s)/f) / y!)
sum_{i=1}^y c_{yi}(f) m^i (s+1)^{if-y} / s^{i(f-1)}}
for \eqn{y = 0, \ldots}, where \code{c_{yi}(f)} are coefficients
obtained by recursion.
}
\arguments{
\item{y}{vector of counts}
\item{q}{vector of quantiles}
\item{p}{vector of probabilities}
\item{n}{number of values to generate}
\item{m}{scalar or vector of means}
\item{s}{scalar or vector of overdispersion parameters}
\item{f}{scalar or vector of family parameters, all < 1}
\item{log}{if TRUE, log probabilities are supplied.}
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dpois}} for the Poisson, \code{\link[rmutil]{ddoublepois}} for
the double Poisson, \code{\link[rmutil]{dmultpois}} for
the multiplicative Poisson, \code{\link[rmutil]{dconsul}} for the Consul
generalized Poisson, \code{\link[rmutil]{dgammacount}} for the gamma count, and
\code{\link{dnbinom}} for the negative binomial distribution.
}
\examples{
dpvfpois(5,10,0.9,0.5)
ppvfpois(5,10,0.9,0.5)
qpvfpois(0.85,10,0.9,0.5)
rpvfpois(10,10,0.9,0.5)
}
\keyword{distribution}
