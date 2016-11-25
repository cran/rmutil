\name{Generalized Logistic}
\title{Generalized Logistic Distribution}
\usage{
dglogis(y, m=0, s=1, f=1, log=FALSE)
pglogis(q, m=0, s=1, f=1)
qglogis(p, m=0, s=1, f=1)
rglogis(n, m=0, s=1, f=1)
}
\alias{dglogis}
\alias{pglogis}
\alias{qglogis}
\alias{rglogis}
\description{
These functions provide information about the generalized logistic
distribution with location parameter equal to \code{m}, dispersion equal
to \code{s}, and family parameter equal to \code{f}: density,
cumulative distribution, quantiles, log hazard, and random generation.

The generalized logistic distribution has density
\deqn{
f(y) =
\frac{\nu \sqrt{3} \exp(-\sqrt{3} (y-\mu)/(\sigma \pi))}{
\sigma \pi (1+\exp(-\sqrt{3} (y-\mu)/(\sigma \pi)))^{\nu+1}}}{
f(y) = f sqrt(3) exp(-sqrt(3) (y-m)/(s pi))/
(s pi (1+exp(-sqrt(3) (y-m)/(s pi)))^(f+1))}
	
where \eqn{\mu}{m} is the location parameter of the distribution,
\eqn{\sigma}{s} is the dispersion, and \eqn{\nu}{f} is the family
parameter.

\eqn{\nu=1}{f=1} gives a logistic distribution.
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
\seealso{
\code{\link{dlogis}} for the logistic distribution.
}
\examples{
dglogis(5, 5, 1, 2)
pglogis(5, 5, 1, 2)
qglogis(0.25, 5, 1, 2)
rglogis(10, 5, 1, 2)
}
\keyword{distribution}
