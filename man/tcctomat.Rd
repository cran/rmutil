\name{tcctomat}
\title{Create or add to a time-constant covariate (tccov) object}
\usage{
tcctomat(ccov, names=NULL, oldccov=NULL, dataframe=TRUE)
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
\item{dataframe}{If TRUE and factor variables are present, the
covariates are stored as a dataframe; if FALSE, they are expanded to
indicator variables. If no factor variables are present, covariates
are always stored as a matrix.}
}
\description{
\code{tcctomat} creates an object of class, tccov, from a vector or
matrix containing time-constant baseline covariates or a model formula,
or combines two such objects.

Such objects can be printed. Methods are available for extracting the
covariates, their names, and the formula: \code{covariates},
\code{names}, and \code{formula}. The method, \code{link{transform}},
can transform variables in place or by adding new variables to the
object.

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
\code{\link{transform}}, \code{\link{tvctomat}}.
}
\examples{
x1 <- gl(4,1)
print(tcc1 <- tcctomat(~x1))
covariates(tcc1)
covariates(tcc1, name="x12")
tcctomat(x1)
tcctomat(x1, dataframe=T)
x2 <- c(0,0,1,1)
print(tcc2 <- tcctomat(~x2))
covariates(tcc2)
print(tcc3 <- tcctomat(~x1+x2))
covariates(tcc3)
covariates(tcc3, names=c("x12","x2"))
formula(tcc3)
names(tcc3)
print(tcc4 <- tcctomat(data.frame(x1,x2)))
covariates(tcc4)
print(tcc5 <- tcctomat(data.frame(x1,x2), dataframe=T))
covariates(tcc5)
}
