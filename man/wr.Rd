\name{wr}
\title{Find the response vector and design matrix for a model formula}
\usage{
wr(formula)
}
\arguments{
\item{formula}{A model formula.}
}
\description{
\code{wr} gives the response vector and design matrix for a formula in
Wilkinson and Rogers notation.
}
\value{
\code{wr} returns a list containing the response vector (z$response)
and the design matrix (z$design) for the formula supplied.
e}
\author{J.K. Lindsey}
\examples{
y <- rnorm(20)
x <- gl(4,5)
z <- rpois(20,2)
wr(y~x+z)
}
