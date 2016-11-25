\name{gettvc}
\title{Find the Most Recent Value of a Time-varying Covariate
Before Each Observed Response}
\alias{gettvc}
\usage{
gettvc(response, times=NULL, tvcov=NULL, tvctimes=NULL,
	oldtvcov=NULL, ties=TRUE)
}
\description{
\code{gettvc} finds the most recent value of a time-varying covariate
before each observed response and possibly adds them to a list of other
time-varying covariates.
It compares the times of response observations with those of
time-varying covariates to find the most recent observed time-varying
covariate for each response. These are either placed in a new object of
class, \code{tvcov}, added to an already existing list of matrices containing
other time-varying covariates and a new object of class, \code{tvcov},
created, or added to an existing object of class, \code{tvcov}.

If there are response observation times before the first covariate
time, the covariate for these times is set to zero.
}
\arguments{
\item{response}{A list of two column matrices with response values and
times for each individual, one matrix or dataframe of response
values, or an object of class, \code{response} (created by
\code{\link[rmutil]{restovec}}).}
\item{times}{When \code{response} is a matrix, a vector of possibly
unequally spaced times for the response, when they are the same for
all individuals or a matrix of times. Not necessary if equally spaced.}
\item{tvcov}{A list of two column matrices with time-varying covariate
values and corresponding times for each individual or one matrix or
dataframe of such covariate values. Times need not be the same as for
responses.}
\item{tvctimes}{When the time-varying covariate is a matrix, a vector
of possibly unequally spaced times for the covariate, when they are
the same for all individuals or a matrix of times. Not necessary if
equally spaced.}
\item{oldtvcov}{A list of matrices with time-varying covariate values,
observed at the event times in \code{response}, for each individual,
or an object of class, \code{tvcov}. If not provided, a new object is
created.}
\item{ties}{If TRUE, when the response and covariate times are
identical, the response depends on that new value (as in observational
studies); if FALSE, only the next response depends on that value (for
example, if the covariate is a new treatment just applied at that time).}
}
\value{
An object of class, \code{tvcov}, is returned containing the new time-varying
covariate and, possibly, those in \code{oldtvcov}.
}
\author{J.K. Lindsey and D.F. Heitjan}
\seealso{
\code{\link[rmutil]{read.list}}, \code{\link[rmutil]{restovec}},
\code{\link[rmutil]{tvctomat}}.
}
\examples{
\dontrun{
y <- matrix(rnorm(20), ncol=5)
resp <- restovec(y, times=c(1,3,6,10,15))
z <- matrix(rpois(20,5),ncol=5)
z
# create a new time-varying covariate object for the response
newtvc <- gettvc(resp, tvcov=z, tvctimes=c(1,2,5,12,14))
covariates(newtvc)
# add another time-varying covariate to the object
z2 <- matrix(rpois(20,5),ncol=5)
z2
newtvc2 <- gettvc(resp, tvcov=z2, tvctimes=c(0,4,5,12,16), oldtvc=newtvc)
covariates(newtvc2)
}
}
\keyword{manip}
