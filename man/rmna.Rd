\name{rmna}
\title{Create a repeated object, removing NAs}
\usage{
rmna(response, tvcov=NULL, ccov=NULL)
}
\alias{rmna}
\alias{print.repeated}
\alias{plot.repeated}
\arguments{
\item{response}{An object of class, response (created by
\code{\link{restovec}}), containing the response variable information.}
\item{tvcov}{An object of class, tvcov (created by \code{\link{tvctomat}}),
containing the time-varying covariate information.}
\item{tccov}{An object of class, tccov (created by \code{\link{tcctomat}}),
containing the time-constant covariate information.}
}
\description{
\code{rmna} forms an object of class, repeated, from a response object
and possibly time-varying covariate (tvcov), and time-constant
covariate (tccov) objects, removing any response and covariate values
that have NAs.

Such objects can be printed and plotted.
}
\value{
Returns an object of class, repeated, containing a list of the
response object (z$response, so that, for example, the response vector
is z$response$y; see \code{\link{restovec}}), and possibly the two
classes of covariate objects (z$ccov and z$tvcov).

Methods are available for extracting the response, the numbers of
observations per individual, the times, the weights, the nesting
variable, and the covariates or their names: \code{response},
\code{nobs}, \code{times}, \code{weights}, \code{nesting},
\code{covariates}, and \code{names}.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{carma}}, \code{\link{elliptic}}, \code{\link{gettvc}},
\code{\link{kalcount}}, \code{\link{kalseries}}, \code{\link{nbkal}},
\code{\link{read.list}}, \code{\link{restovec}}, \code{\link{tcctomat}},
\code{\link{tvctomat}}.
}
\examples{
y <- matrix(rnorm(20),ncol=5)
tt <- c(1,3,6,10,15)
print(resp <- restovec(y,times=tt))
x <- c(0,0,1,1)
tcc <- tcctomat(x)
z <- matrix(rpois(20,5),ncol=5)
tvc <- tvctomat(z)
print(reps <- rmna(resp, tvcov=tvc, ccov=tcc))
response(reps)
response(reps, nind=2:3)
times(reps)
nobs(reps)
weights(reps)
covariates(reps)
covariates(reps,names="x")
covariates(reps,names="z")
names(reps)
nesting(reps)
# because individuals are the only nesting, this is the same as
covind(reps)
# binomial
y <- matrix(rpois(20,5),ncol=5)
print(respb <- restovec(y,totals=y+matrix(rpois(20,5),ncol=5),times=tt))
print(repsb <- rmna(respb, tvcov=tvc, ccov=tcc))
response(repsb)
# censored data
y <- matrix(rweibull(20,2,5),ncol=5)
print(respc <- restovec(y,censor=matrix(rbinom(20,1,0.9),ncol=5),times=tt))
print(repsc <- rmna(respc, tvcov=tvc, ccov=tcc))
# if there is no censoring, censor indicator is not printed
response(repsc)
# nesting clustered within individuals
nest <- c(1,1,2,2,2)
print(respn <- restovec(y,censor=matrix(rbinom(20,1,0.9),ncol=5),
	times=tt,nest=nest))
print(repsn <- rmna(respn, tvcov=tvc, ccov=tcc))
response(respn)
times(respn)
nesting(respn)
}
