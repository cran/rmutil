\name{Inverse Gaussian}
\title{Inverse Gaussian Distribution}
\usage{
dinvgauss(y, m, s, log=FALSE)
pinvgauss(q, m, s)
qinvgauss(p, m, s)
rinvgauss(n, m, s)
}
\alias{dinvgauss}
\alias{pinvgauss}
\alias{qinvgauss}
\alias{rinvgauss}
\description{
These functions provide information about the inverse Gaussian
distribution with mean equal to \code{m} and dispersion equal to
\code{s}: density, cumulative distribution, quantiles, log hazard, and
random generation.

The inverse Gaussian distribution has density
\deqn{
f(y) =
\frac{1}{\sqrt{2\pi\sigma y^3}} e^{-(y-\mu)^2/(2 y \sigma m^2)}}{
f(y) = 1/sqrt(2 pi s y^3) e^-((y - m)^2/(2 y s m^2))}
where \eqn{\mu}{m} is the mean of the distribution and
\eqn{\sigma}{s} is the dispersion.
}
\arguments{
\item{y}{vector of responses.}
\item{q}{vector of quantiles.}
\item{p}{vector of probabilities}
\item{n}{number of values to generate}
\item{m}{vector of means.}
\item{s}{vector of dispersion parameters.}
\item{log}{if TRUE, log probabilities are supplied.}
}
\author{J.K. Lindsey}
\seealso{
\code{\link{dnorm}} for the normal distribution and
\code{\link{dlnorm}} for the \emph{Log}normal distribution.
}
\examples{
dinvgauss(5, 5, 1)
pinvgauss(5, 5, 1)
qinvgauss(0.8, 5, 1)
rinvgauss(10, 5, 1)
}
\keyword{distribution}
