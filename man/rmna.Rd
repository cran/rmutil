\name{rmna}
\title{Create a repeated Object, Removing NAs}
\alias{rmna}
\usage{
rmna(response, ccov=NULL, tvcov=NULL)}
\description{
\code{rmna} forms an object of class, \code{repeated}, from a
\code{response} object and possibly time-varying or
intra-individual covariate (\code{tvcov}), and time-constant or
inter-individual covariate (\code{tccov}) objects, removing any
observations where response and covariate values have NAs. Subjects
must be in the same order in all (three) objects to be combined.

Such objects can be printed and plotted. Methods are available for
extracting the response, the numbers of observations per individual,
the times, the weights, the units of measurement/Jacobian, the nesting
variable, the covariates, and their names: \code{\link[rmutil]{response}},
\code{\link[rmutil]{nobs}}, \code{\link[rmutil]{times}},
\code{\link[rmutil]{weights}}, \code{\link[rmutil]{delta}},
\code{\link[rmutil]{nesting}}, \code{\link[rmutil]{covariates}}, and
\code{\link[rmutil]{names}}.
}
\arguments{
\item{response}{An object of class, \code{response} (created by
\code{\link[rmutil]{restovec}}), containing the response variable information.}
\item{ccov}{An object of class, \code{tccov} (created by
\code{\link[rmutil]{tcctomat}}), containing the time-constant or
inter-individual covariate information.}
\item{tvcov}{An object of class, \code{tvcov} (created by
\code{\link[rmutil]{tvctomat}}), containing the time-varying or
intra-individual covariate information.}
}
\value{
Returns an object of class, \code{repeated}, containing a list of the
response object (\code{z$response}, so that, for example, the response vector
is \code{z$response$y}; see \code{\link[rmutil]{restovec}}), and
possibly the two classes of covariate objects (\code{z$ccov} and
\code{z$tvcov}; see \code{\link[rmutil]{tcctomat}} and
\code{\link[rmutil]{tvctomat}}).
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{DataMethods}}, \code{\link[rmutil]{covariates}},
\code{\link[rmutil]{covind}}, \code{\link[rmutil]{delta}},
\code{\link[rmutil]{dftorep}}, \code{\link[rmutil]{lvna}},
\code{\link[rmutil]{names}}, \code{\link[rmutil]{nesting}},
\code{\link[rmutil]{nobs}}, \code{\link[rmutil]{read.list}},
\code{\link[rmutil]{read.surv}}, \code{\link[rmutil]{response}},
\code{\link[rmutil]{resptype}}, \code{\link[rmutil]{restovec}},
\code{\link[rmutil]{tcctomat}}, \code{\link[rmutil]{times}},
\code{\link[rmutil]{transform}}, \code{\link[rmutil]{tvctomat}},
\code{\link[rmutil]{units}}, \code{\link[rmutil]{weights}}
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
#
# use in glm
rm(y,x,z)
glm(y~x+z,data=as.data.frame(reps))
#
# binomial
y <- matrix(rpois(20,5),ncol=5)
print(respb <- restovec(y,totals=y+matrix(rpois(20,5),ncol=5),times=tt))
print(repsb <- rmna(respb, tvcov=tvc, ccov=tcc))
response(repsb)
#
# censored data
y <- matrix(rweibull(20,2,5),ncol=5)
print(respc <- restovec(y,censor=matrix(rbinom(20,1,0.9),ncol=5),times=tt))
print(repsc <- rmna(respc, tvcov=tvc, ccov=tcc))
# if there is no censoring, censor indicator is not printed
response(repsc)
#
# nesting clustered within individuals
nest <- c(1,1,2,2,2)
print(respn <- restovec(y,censor=matrix(rbinom(20,1,0.9),ncol=5),
	times=tt,nest=nest))
print(repsn <- rmna(respn, tvcov=tvc, ccov=tcc))
response(respn)
times(respn)
nesting(respn)
}
\keyword{manip}
