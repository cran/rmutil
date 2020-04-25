\name{tvctomat}
\title{Create a Time-varying, Intra-individual Covariate (tvcov) Object}
\alias{tvctomat}
\usage{
tvctomat(tvcov, names=NULL, units=NULL, interaction=NULL, ccov=NULL,
	oldtvcov=NULL, dataframe=TRUE, description=NULL)
}
\description{
\code{tvctovmat} creates an object of class, \code{tvcov}, from a list of
matrices with time-varying or intra-individual covariates for each
individual or one matrix or dataframe of such covariate values. It can
also combine two such objects or add interactions among covariates.

Such objects can be printed. Methods are available for extracting the
covariates and their names: \code{\link[rmutil]{covariates}} and
\code{\link[rmutil]{names}}. The method,
\code{\link[rmutil]{transform}}, can transform variables in place or
by adding new variables to the object.
}
\arguments{
\item{tvcov}{Either (1) if unbalanced, a list of matrices or
dataframes with time-varying or intra-individual covariate values for
each individual (one column per variable), (2) if balanced, one matrix
or dataframe of such covariate values (when there is only one such
covariate) with dimensions: number of individuals by number of
observations/individual, or (3) an object of class, \code{tvcov}.
In the first two cases, the variables may be factors; if
\code{dataframe=FALSE}, these are transformed to indicator variables.}
\item{names}{The names of the time-varying or intra-individual
covariates in \code{tvcov} (if the matrices do not have column names)
or the names of the time-constant covariates for which interactions
are to be created.}
\item{units}{Optional character vector specifying units of
measurements of covariates.}
\item{interaction}{A pair of index numbers or names of variables in
\code{tvcov}, with that class, for which an interaction is to be added
or, if \code{ccov} is provided, a set of such names of time-varying or
intra-individual covariates for creating interactions with the
time-constant covariates.}
\item{ccov}{Time-constant or inter-individual covariates for which an
interaction is to be introduced with time-varying or intra-individual
covariates in \code{tvcov}.}
\item{oldtvcov}{An object of class, \code{tvcov}, to which
\code{tvcov} is to be added.}
\item{dataframe}{If TRUE and factor variables are present, the
covariates are stored as a dataframe; if FALSE, they are expanded to
indicator variables. If no factor variables are present, covariates
are always stored as a matrix.}
\item{description}{An optional named list of character vectors with
names of some or all covariates containing their descriptions.}
}
\value{
Returns an object of class, \code{tvcov}, containing a matrix or
dataframe for the covariates (\code{z$tvcov}) with one row per
response per individual and a vector giving the number of observations
per individual (\code{z$nobs}).
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{DataMethods}}, \code{\link[rmutil]{covariates}},
\code{\link[rmutil]{description}}, \code{\link[rmutil]{formula}},
\code{\link[rmutil]{gettvc}}, \code{\link[rmutil]{lvna}},
\code{\link[rmutil]{names}}, \code{\link[rmutil]{restovec}},
\code{\link[rmutil]{rmna}}, \code{\link[rmutil]{tcctomat}},
\code{\link[rmutil]{transform}}, \code{\link[rmutil]{units}}
}
\examples{
z <- matrix(rpois(20,5),ncol=5)
print(tvc <- tvctomat(z, units="days"))
covariates(tvc)
names(tvc)
v <- data.frame(matrix(rep(c("a","b","c","d","e"),4),ncol=5),stringsAsFactors=TRUE)
print(tvc2 <- tvctomat(v, oldtvc=tvc, units=NA))
covariates(tvc2)
print(tvc3 <- tvctomat(v, oldtvc=tvc, dataframe=FALSE, units=NA))
covariates(tvc3)
print(tvc4 <- tvctomat(tvc2, interaction=c("z","v")))
covariates(tvc4)
x1 <- 1:4
x2 <- gl(4,1)
xx <- tcctomat(data.frame(x1,x2), dataframe=FALSE)
tvctomat(tvc3, interaction="z", ccov=xx)
tvctomat(tvc3, interaction="z", ccov=xx, names="x1")
tvctomat(tvc3, interaction="z", ccov=xx, names=c("x22","x23","x24"))
xx <- tcctomat(data.frame(x1,x2), dataframe=TRUE)
tvctomat(tvc2, interaction="z", ccov=xx)
tvctomat(tvc2, interaction="z", ccov=xx, names="x1")
tvctomat(tvc2, interaction="z", ccov=xx, names="x2")
}
\keyword{manip}
