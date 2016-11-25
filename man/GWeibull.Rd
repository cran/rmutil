\name{Generalized Weibull}
\title{Generalized Weibull Distribution}
\usage{
dgweibull(y, s, m, f, log=FALSE)
pgweibull(q, s, m, f)
qgweibull(p, s, m, f)
rgweibull(n, s, m, f)
}
\alias{dgweibull}
\alias{pgweibull}
\alias{qgweibull}
\alias{rgweibull}
\description{
These functions provide information about the generalized Weibull
distribution, also called the exponentiated Weibull, with scale
parameter equal to \code{m}, shape equal to \code{s}, and family
parameter equal to \code{f}: density, cumulative distribution,
quantiles, log hazard, and random generation.

The generalized Weibull distribution has density
\deqn{
f(y) = \frac{\sigma \nu y^{\sigma-1} (1-\exp(-(y/\mu)^\sigma))^{\nu-1}
\exp(-(y/\mu)^\sigma)}{\mu^\sigma}}{
f(y) = s f y^(s-1) (1-exp(-(y/m)^s))^(f-1) exp(-(y/m)^s)/m^s}
	
where \eqn{\mu}{m} is the scale parameter of the distribution,
\eqn{\sigma}{s} is the shape, and \eqn{\nu}{f} is the family
parameter.

\eqn{\nu=1}{f=1} gives a Weibull distribution, for
\eqn{\sigma=1}{s=1}, \eqn{\nu<0}{f<0} a generalized F distribution,
and for \eqn{\sigma>0}{s>0}, \eqn{\nu\leq0}{f<=0} a Burr type XII distribution.
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
\code{\link{dweibull}} for the Weibull distribution,
\code{\link{df}} for the F distribution,
\code{\link[rmutil]{dburr}} for the Burr distribution.
}
\examples{
dgweibull(5, 1, 3, 2)
pgweibull(5, 1, 3, 2)
qgweibull(0.65, 1, 3, 2)
rgweibull(10, 1, 3, 2)
}
\keyword{distribution}
