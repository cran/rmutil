\name{tvctomat}
\title{Create or add to a time-varying covariate (tvcov) object}
\usage{
tvctomat(tvcov, names=NULL, index=NULL, tccov=NULL, oldtvcov=NULL)
}
\alias{tvctomat}
\alias{print.tvccov}
\arguments{
\item{tvcov}{A list of matrices with time-varying covariate values
for each individual (one column per variable), one matrix or
dataframe of such covariate values (when there is only one such
covariate), or an object of class, tvcov.}
\item{names}{The names of the time-varying covariates in tvcov (if the
matrices do not have column names) or the name of the time-constant
covariate.}
\item{index}{The index numbers or names of variables in tvcov, with
that class, for which an interaction is to be added or one such number
or name for an interaction with a time-constant covariate.}
\index{tccov}{A time-constant covariate for which an interaction is to
be introduced with a time-varying covariate in tvcov.}
\item{oldtvcov}{An object of class, tvcov, to which tvcov is to be
added.}
}
\description{
\code{tvctovmat} transforms a list of matrices with time-varying
covariates for each individual or one matrix or dataframe of such
covariate values into an object of class, tvcov, or combines two such
objects. It can also add interactions among covariates.

Such objects can be printed. Methods are available for extracting the
covariates and their names.
}
\value{
Returns an object of class, tvcov, containing a matrix for the
covariates (z$tvcov) with one row per response per individual and a
vector giving the number of observations per individual (z$nobs).
}
\author{J.K. Lindsey}
\seealso{
\code{\link{gettvc}}, \code{\link{read.list}}, \code{\link{restovec}},
\code{\link{rmna}}, \code{\link{tcctomat}}.
}

\author{J.K. Lindsey}
\examples{
z <- matrix(rpois(20,5),ncol=5)
tvctomat(z)
}
