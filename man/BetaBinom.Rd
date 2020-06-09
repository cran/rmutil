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
Compared to the parameterization of `VGAM::pbetabinom.ab`,
\code{m = alpha/(alpha+beta)} and \code{s = (alpha+beta)}.
See examples.

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
\item{m}{vector of probabilities of success; Compared to the parameterization of `VGAM::pbetabinom.ab`,
\code{m = alpha/(alpha+beta)} where \code{shape1=alpha} and \code{shape2=beta}. See examples.}
\item{s}{vector of overdispersion parameters; Compared to the parameterization of `VGAM::pbetabinom.ab`, \code{s = (alpha+beta)} where \code{shape1=alpha} and \code{shape2=beta}. See examples.}
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
## compare to VGAM
\dontrun{
# The beta binomial distribution with total = n and prob = m has density
# 
# p(y) = B(y+s m,n-y+s (1-m)) Choose(n,y) / B(s m,s (1-m))
# 
# for y = 0, …, n where B() is the beta function.

## in `rmutil` from the .Rd file (excerpt above), the "alpha" is s*m
## in `rmutil` from the .Rd file (excerpt above), the "beta"  is s*(1-m)

## in `VGAM`, rho is 1/(1+alpha+beta)

qq = 2.2
zz = 100

alpha = 1.1
beta  = 2
VGAM::pbetabinom.ab(q=qq, size=zz, shape1=alpha, shape2=beta)

## for VGAM `rho`
rr = 1/(1+alpha+beta)
VGAM::pbetabinom   (q=qq, size=zz, prob=mm, rho = rr)

## for rmutil `m` and `s`:
mm = alpha / (alpha+beta)
ss = (alpha+beta)
rmutil::pbetabinom(q=qq, size=zz, m=mm, s=ss )
}
}
\keyword{distribution}
