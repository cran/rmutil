\name{gauss.hermite}
\title{Calculate Gauss-Hermite Quadrature Points}
\usage{
gauss.hermite(points, iterlim=10)
}
\arguments{
\item{points}{The number of points.}
\item{iterlim}{Maximum number of iterations in Newton-Raphson.}
}
\description{
\code{gauss.hermite} calculates the Gauss-Hermite quadrature points.
}
\value{
\code{gauss.hermite} returns a two-column matrix containing the points
and their corresponding weights.
}
\author{J.K. Lindsey}
\examples{
gauss.hermite(10)
}
