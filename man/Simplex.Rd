\name{Simplex}
\title{Simplex Distribution}
\usage{
dsimplex(y, m, s, log=FALSE)
psimplex(q, m, s)
qsimplex(p, m, s)
rsimplex(n, m, s)
}
\alias{dsimplex}
\alias{psimplex}
\alias{qsimplex}
\alias{rsimplex}
\description{
These functions provide information about the simplex distribution
with location parameter equal to \code{m} and shape equal to
\code{s}: density, cumulative distribution, quantiles, and
random generation.

The simplex distribution has density
\deqn{
f(y) = \frac{1}{\sqrt(2\pi\sigma(y(1-y))^3)}
     \exp(-((y-\mu)/(\mu(1-\mu)))^2/(2y(1-y)\sigma))}{
f(y) = exp(-((y-m)/(m(1-m)))^2/(2y(1-y)s))/sqrt(2 pi s(y(1-y))^3)}
where \eqn{\mu}{m} is the location parameter of the distribution and
\eqn{\sigma}{s} is the shape, and \eqn{0<y<1}{0<y<1}.
}
\arguments{
\item{y}{vector of responses.}
\item{q}{vector of quantiles.}
\item{p}{vector of probabilities}
\item{n}{number of values to generate}
\item{m}{vector of location parameters.}
\item{s}{vector of shape parameters.}
\item{log}{if TRUE, log probabilities are supplied.}
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dbeta}} for the beta distribution and
\code{\link{dtwosidedpower}} for the two-sided power distribution,
other distributions for proportions between zero and one.
}
\examples{
dsimplex(0.3, 0.5, 1)
psimplex(0.3, 0.5, 1)
qsimplex(0.1, 0.5, 1)
rsimplex(10, 0.5, 1)
}
\keyword{distribution}
