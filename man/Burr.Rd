\name{Burr}
\title{Burr Distribution}
\usage{
dburr(y, m, s, f, log=FALSE)
pburr(q, m, s, f)
qburr(p, m, s, f)
rburr(n, m, s, f)
}
\alias{dburr}
\alias{pburr}
\alias{qburr}
\alias{rburr}
\description{
These functions provide information about the Burr distribution with
location parameter equal to \code{m}, dispersion equal to
\code{s}, and family parameter equal to \code{f}: density,
cumulative distribution, quantiles, log hazard, and random generation.

The Burr distribution has density
\deqn{
f(y) = \frac{\nu \sigma (y / \mu)^{\sigma-1}}
{\mu (1+(y/\mu)^\sigma)^{\nu+1}}}{
f(y) = f s (y/m)^(s-1)/(m (1+(y/m)^s)^(f+1))}
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
dburr(2, 5, 1, 2)
pburr(2, 5, 1, 2)
qburr(0.3, 5, 1, 2)
rburr(10, 5, 1, 2)
}
\keyword{distribution}
