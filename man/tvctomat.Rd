\name{tvctomat}
\title{Create or add to a time-varying covariate (tvcov) object}
\usage{
tvctomat(tvcov, names=NULL, interaction=NULL, ccov=NULL,
	oldtvcov=NULL, dataframe=TRUE)
}
\alias{tvctomat}
\alias{print.tvccov}
\arguments{
\item{tvcov}{A list of matrices with time-varying covariate values
for each individual (one column per variable), one matrix or
dataframe of such covariate values (when there is only one such
covariate), or an object of class, tvcov. In the first two cases, the
variables may be factors.}
\item{names}{The names of the time-varying covariates in tvcov (if the
matrices do not have column names) or the names of the time-constant
covariates for interactions.}
\item{interaction}{A pair of index numbers or names of variables in tvcov,
with that class, for which an interaction is to be added or, if
\code{ccov} is provided, a set of such names of time-varying
covariates for creating interactions with the time-constant covariates.}
\item{ccov}{Time-constant covariates for which an interaction is to
be introduced with time-varying covariates in tvcov.}
\item{oldtvcov}{An object of class, tvcov, to which tvcov is to be
added.}
\item{dataframe}{If TRUE and factor variables are present, the
covariates are stored as a dataframe; if FALSE, they are expanded to
indicator variables. If no factor variables are present, covariates
are always stored as a matrix.}
}
\description{
\code{tvctovmat} creates an object of class, tvcov, from a list of
matrices with time-varying covariates for each individual or one
matrix or dataframe of such covariate values or combines two such
objects. It can also add interactions among covariates.

Such objects can be printed. Methods are available for extracting the
covariates and their names: \code{covariates} and \code{names}. The
method, \code{link{transform}}, can transform variables in place or by
adding new variables to the object.
}
\value{
Returns an object of class, tvcov, containing a matrix for the
covariates (z$tvcov) with one row per response per individual and a
vector giving the number of observations per individual (z$nobs).
}
\author{J.K. Lindsey}
\seealso{
\code{\link{gettvc}}, \code{\link{read.list}}, \code{\link{restovec}},
\code{\link{rmna}}, \code{\link{tcctomat}}, \code{\link{transform}}.
}

\author{J.K. Lindsey}
\examples{
z <- matrix(rpois(20,5),ncol=5)
print(tvc <- tvctomat(z))
covariates(tvc)
names(tvc)
v <- data.frame(matrix(rep(c("a","b","c","d","e"),4),ncol=5))
print(tvc2 <- tvctomat(v, oldtvc=tvc))
covariates(tvc2)
print(tvc3 <- tvctomat(v, oldtvc=tvc, dataframe=T))
covariates(tvc3)
print(tvc4 <- tvctomat(tvc3, interaction=c("z","v")))
covariates(tvc4)
x1 <- 1:4
x2 <- gl(4,1)
xx <- tcctomat(data.frame(x1,x2))
tvctomat(tvc3, interaction="z", ccov=xx)
tvctomat(tvc3, interaction="z", ccov=xx, names="x1")
tvctomat(tvc3, interaction="z", ccov=xx, names=c("x22","x23","x24"))
xx <- tcctomat(data.frame(x1,x2), dataframe=T)
tvctomat(tvc3, interaction="z", ccov=xx)
tvctomat(tvc3, interaction="z", ccov=xx, names="x1")
tvctomat(tvc3, interaction="z", ccov=xx, names="x2")
}
