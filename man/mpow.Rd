\name{mpower}
\title{Power of a Matrix}
\alias{\%^\%}
\usage{
x\%^\%p
}
\description{
\code{\%^\%} calculates \code{x^p} for the square matrix, \code{x}, by
spectral decomposition.
}
\arguments{
\item{x}{A square matrix.}
\item{p}{The power to which the matrix is to be raised.}
}
\value{
\code{\%^\%} returns the power of a matrix.
}
\author{J.K. Lindsey}
\examples{
\dontrun{
x <- matrix(c(0.4,0.6,0.6,0.4),nrow=2)
x\%^\%2
x\%^\%10
x\%^\%20
}
}
\keyword{array}
