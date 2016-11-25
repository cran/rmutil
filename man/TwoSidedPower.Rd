\name{Two-Sided Power}
\title{Two-Sided Power Distribution}
\usage{
dtwosidedpower(y, m, s=2, log=FALSE)
ptwosidedpower(q, m, s=2)
qtwosidedpower(p, m, s=2)
rtwosidedpower(n, m, s=2)
}
\alias{dtwosidedpower}
\alias{ptwosidedpower}
\alias{qtwosidedpower}
\alias{rtwosidedpower}
\description{
These functions provide information about the two-sided power distribution
with location parameter equal to \code{m} and shape equal to
\code{s}: density, cumulative distribution, quantiles, and
random generation.

The two-sided power distribution has density
\deqn{
f(y) = s(\frac{y}{m})^{s-1}, y<=m}{
f(y) = s(y/m)^(s-1), y<=m}
\deqn{
f(y) =s(\frac{1-y}{1-m})^{s-1}, y>=m}{
f(y) = s((1-y)/(1-m))^(s-1), y>=m}
where \eqn{\mu}{m} is the location parameter of the distribution and
\eqn{\sigma}{s} is the shape, and \eqn{0<y<1}{0<y<1}.

For \eqn{\sigma=1}{s=1}, this is the uniform distribution and for
\eqn{\sigma=2}{s=2}, it is the triangular distribution.
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
\references{
van Dorp, J.R. and Kotz, S. (2002) A novel extension of the triangular
distribution and its parameter estimation. The Statistician 51, 63-79.
}
\seealso{
\code{\link{dbeta}} for the beta distribution and
\code{\link{dsimplex}} for the simplex distribution, other
distributions for proportions between zero and one.
}
\examples{
dtwosidedpower(0.3, 0.5, 3)
ptwosidedpower(0.3, 0.5, 3)
qtwosidedpower(0.1, 0.5, 3)
rtwosidedpower(10, 0.5, 3)
}
\keyword{distribution}
