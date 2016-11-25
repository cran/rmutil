\name{Box-Cox}
\title{Box-Cox Distribution}
\usage{
dboxcox(y, m, s=1, f=1, log=FALSE)
pboxcox(q, m, s=1, f=1)
qboxcox(p, m, s=1, f=1)
rboxcox(n, m, s=1, f=1)
}
\alias{dboxcox}
\alias{pboxcox}
\alias{qboxcox}
\alias{rboxcox}
\description{
These functions provide information about the Box-Cox
distribution with location parameter equal to \code{m}, dispersion
equal to \code{s}, and power transformation equal to \code{f}: density,
cumulative distribution, quantiles, log hazard, and random generation.

The Box-Cox distribution has density
\deqn{
f(y) =
\frac{1}{\sqrt{2 \pi \sigma^2}} \exp(-((y^\nu/\nu-\mu)^2/(2 \sigma^2)))/
(1-I(\nu<0)-sign(\nu)*pnorm(0,\mu,sqrt(\sigma)))}{
f(y) = 1/sqrt(2 pi s^2) exp(-((y^f/f - mu)^2/(2 s^2)))/
(1-I(f<0)-sign(f)*pnorm(0,m,sqrt(s)))}
where \eqn{\mu}{m} is the location parameter of the distribution,
\eqn{\sigma}{s} is the dispersion, \eqn{\nu}{f} is the family
parameter, \eqn{I()} is the indicator function, and \eqn{y>0}.

\eqn{\nu=1}{f=1} gives a truncated normal distribution.
}
\arguments{
\item{y}{vector of responses.}
\item{q}{vector of quantiles.}
\item{p}{vector of probabilities}
\item{n}{number of values to generate}
\item{m}{vector of location parameters.}
\item{s}{vector of dispersion parameters.}
\item{f}{vector of power parameters.}
\item{log}{if TRUE, log probabilities are supplied.}
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dnorm}} for the normal or Gaussian distribution.
}
\examples{
dboxcox(2, 5, 5, 2)
pboxcox(2, 5, 5, 2)
qboxcox(0.1, 5, 5, 2)
rboxcox(10, 5, 5, 2)
}
\keyword{distribution}
