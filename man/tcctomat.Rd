\name{tcctomat}
\title{Create or add to a time-constant covariate (tccov) object}
\usage{
tcctomat(ccov, names=NULL, oldccov=NULL)
}
\alias{tcctomat}
\alias{print.tccov}
\alias{covind}
\arguments{
\item{ccov}{A vector or matrix containing time-constant baseline
covariates with one row per individual, a model formula using
vectors of the same size, or an object of class, tccov.}
\item{names}{The names of the covariates (if the matrix does not have
column names).}
\item{oldccov}{An object of class, tccov, to which ccov is to be
added.}
}
\description{
\code{tcctomat} transforms a vector or matrix containing time-constant
baseline covariates or a model formula into an object of class, tccov,
or combines two such objects.

Such objects can be printed. Methods are available for extracting the
covariates, their names, and the formula.

To obtain the indexing to expand time-constant covariates to the size
of a repeated measurements response, use \code{covind}.
}
\value{
Returns an object of class, tccov, containing one matrix for the
covariates (z$ccov) with one row per individual and possibly the model
formula (z$linear).
}
\author{J.K. Lindsey}
\seealso{
\code{\link{read.list}}, \code{\link{restovec}}, \code{\link{rmna}},
\code{\link{tvctomat}}.
}
\examples{
x1 <- gl(4,1)
tcctomat(~x1)
x2 <- c(0,0,1,1)
tcctomat(~x2)
}
