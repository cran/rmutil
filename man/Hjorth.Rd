\name{Hjorth}
\title{Hjorth Distribution}
\usage{
dhjorth(y, m, s, f, log=FALSE)
phjorth(q, m, s, f)
qhjorth(p, m, s, f)
rhjorth(n, m, s, f)
}
\alias{dhjorth}
\alias{phjorth}
\alias{qhjorth}
\alias{rhjorth}
\description{
These functions provide information about the Hjorth
distribution with location parameter equal to \code{m}, dispersion equal
to \code{s}, and family parameter equal to \code{f}: density,
cumulative distribution, quantiles, log hazard, and random generation.

The Hjorth distribution has density
\deqn{
f(y) = (1+\sigma y)^{-\nu/\sigma} \exp(-(y/\mu)^2/2)
(\frac{y}{\mu^2}+\frac{\nu}{1+\sigma y})}{
f(y) = (1+s y)^(-f/s) exp(-(y/m)^2/2) (y/m^2+f/(1+s y))}
	
where \eqn{\mu}{m} is the location parameter of the distribution,
\eqn{\sigma}{s} is the dispersion, and \eqn{\nu}{f} is the family
parameter.
}
\arguments{
\item{y}{vector of responses.}
\item{q}{vector of quantiles.}
\item{p}{vector of probabilities}
\item{n}{number of values to generate}
\item{m}{vector of location parameters.}
\item{s}{vector of dispersion parameters.}
\item{f}{vector of family parameters.}
\item{log}{if TRUE, log probabilities are supplied.}
}
\author{J.K. Lindsey}
\examples{
dhjorth(5, 5, 5, 2)
phjorth(5, 5, 5, 2)
qhjorth(0.8, 5, 5, 2)
rhjorth(10, 5, 5, 2)
}
\keyword{distribution}
