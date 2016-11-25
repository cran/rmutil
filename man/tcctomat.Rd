\name{tcctomat}
\title{Create a Time-constant, Inter-individual Covariate (tccov) Object}
\alias{tcctomat}
\usage{
tcctomat(ccov, names=NULL, units=NULL, oldccov=NULL, dataframe=TRUE,
	description=NULL)
}
\description{
\code{tcctomat} creates an object of class, \code{tccov}, from a
vector or matrix containing time-constant or inter-individual baseline
covariates or a model formula. It can also combine two such objects.

Such objects can be printed. Methods are available for extracting the
covariates, their names, and the formula: \code{\link[rmutil]{covariates}},
\code{\link[rmutil]{names}}, and \code{\link[rmutil]{formula}}. The method,
\code{\link[rmutil]{transform}}, can transform variables in place or by adding
new variables to the object.

To obtain the indexing to expand time-constant or inter-individual
covariates to the size of a repeated measurements response, use
\code{\link[rmutil]{covind}}.
}
\arguments{
\item{ccov}{A vector, matrix, or dataframe containing time-constant or
inter-individual baseline covariates with one row per individual, a
model formula using vectors of the same size, or an object of class,
\code{tccov}. In the first two cases, the variables may be factors; if
\code{dataframe=FALSE}, these are transformed to indicator variables.}
\item{units}{Optional character vector specifying units of
measurements of covariates.}
\item{names}{The names of the covariates (if the matrix does not have
column names).}
\item{oldccov}{An object of class, \code{tccov}, to which \code{ccov}
is to be added.}
\item{dataframe}{If TRUE and factor variables are present, the
covariates are stored as a dataframe; if FALSE, they are expanded to
indicator variables. If no factor variables are present, covariates
are always stored as a matrix.}
\item{description}{An optional named list of character vectors with
names of some or all covariates containing their descriptions.}
}
\value{
Returns an object of class, \code{tccov}, containing one matrix or
dataframe for the covariates (\code{z$ccov}) with one row per
individual and possibly the model formula (\code{z$linear}).
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{DataMethods}}, \code{\link[rmutil]{covariates}},
\code{\link[rmutil]{description}}, \code{\link[rmutil]{formula}},
\code{\link[rmutil]{lvna}}, \code{\link[rmutil]{names}},
\code{\link[rmutil]{restovec}}, \code{\link[rmutil]{rmna}},
\code{\link[rmutil]{transform}}, \code{\link[rmutil]{tvctomat}},
\code{\link[rmutil]{units}}
}
\examples{
x1 <- gl(4,1)
print(tcc1 <- tcctomat(~x1))
covariates(tcc1)
covariates(tcc1, name="x12")
tcctomat(x1)
tcctomat(x1, dataframe=FALSE)
x2 <- c(0,0,1,1)
print(tcc2 <- tcctomat(~x2, units="days"))
covariates(tcc2)
print(tcc3 <- tcctomat(~x1+x2))
covariates(tcc3)
covariates(tcc3, names=c("x12","x2"))
formula(tcc3)
names(tcc3)
print(tcc4 <- tcctomat(data.frame(x1,x2), units=c(NA,"days")))
covariates(tcc4)
print(tcc5 <- tcctomat(data.frame(x1,x2), dataframe=FALSE, units=c(NA,"days")))
covariates(tcc5)
}
\keyword{manip}
