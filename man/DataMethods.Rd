\name{DataMethods}
\title{Methods for response, tccov, tvcov, and repeated Data Objects}
\usage{
as.data.frame(x, ...)
as.matrix(x, ...)
covariates(z, ...)
covind(z, ...)
delta(z, ...)
\method{formula}{tccov}(x, ...)
\method{formula}{repeated}(x, ...)
\method{names}{tccov}(x, ...)
\method{names}{repeated}(x, ...)
nesting(z, ...)
nobs(z, ...)
\method{plot}{response}(x, name=NULL, nind=NULL, nest=1, ccov=NULL, add=FALSE, lty=NULL, pch=NULL,
	main=NULL, ylim=NULL, xlim=NULL, xlab=NULL, ylab=NULL, ...)
\method{plot}{repeated}(x, name=NULL, nind=NULL, nest=1, ccov=NULL, add=FALSE, lty=NULL, pch=NULL,
	main=NULL, ylim=NULL, xlim=NULL, xlab=NULL, ylab=NULL, ...)
\method{print}{tccov}(x, ...)
\method{print}{repeated}(x, nindmax=50, ...)
response(z, ...)
resptype(z, ...)
times(z, ...)
\method{transform}{response}(`_data`, times=NULL, units=NULL, ...)
\method{transform}{repeated}(`_data`, times=NULL, ...)
units(x, ...)
weights(object, ...)
}
\alias{DataMethods}
\alias{as.data.frame}
\alias{as.data.frame.repeated}
\alias{as.data.frame.response}
\alias{as.data.frame.tccov}
\alias{as.data.frame.tvcov}
\alias{as.matrix}
\alias{as.matrix.repeated}
\alias{as.matrix.response}
\alias{as.matrix.tccov}
\alias{as.matrix.tvcov}
\alias{coef.gnlm}
\alias{covariates}
\alias{covariates.tccov}
\alias{covariates.tvcov}
\alias{covariates.repeated}
\alias{covind}
\alias{covind.default}
\alias{delta}
\alias{delta.response}
\alias{delta.repeated}
\alias{description}
\alias{description.default}
\alias{description.repeated}
\alias{deviance.gnlm}
\alias{df.residual.gnlm}
\alias{formula}
\alias{formula.tccov}
\alias{formula.repeated}
\alias{names}
\alias{names.response}
\alias{names.tccov}
\alias{names.tvcov}
\alias{names.repeated}
\alias{nesting}
\alias{nesting.response}
\alias{nesting.repeated}
\alias{nobs}
\alias{nobs.default}
\alias{nobs.response}
\alias{nobs.tvcov}
\alias{nobs.data.frame}
\alias{plot.response}
\alias{plot.repeated}
\alias{print.response}
\alias{print.tccov}
\alias{print.tvcov}
\alias{print.repeated}
\alias{print.fmobj}
\alias{print.gnlm}
\alias{response}
\alias{response.response}
\alias{response.repeated}
\alias{resptype}
\alias{resptype.response}
\alias{resptype.repeated}
\alias{times}
\alias{times.default}
\alias{times.response}
\alias{transform}
\alias{transform.response}
\alias{transform.tccov}
\alias{transform.tvcov}
\alias{transform.repeated}
\alias{units}
\alias{units.default}
\alias{units.repeated}
\alias{vcov.gnlm}
\alias{weights}
\alias{weights.response}
\alias{weights.repeated}
\alias{weights.gnlm}
\description{
Objects of class, \code{response}, contain response values, and possibly
the corresponding times, binomial totals, nesting categories, censor
indicators, and/or units of precision/Jacobian. Objects of class,
\code{tccov}, contain time-constant or inter-individual, baseline
covariates. Objects of class, \code{tvcov}, contain time-varying or
intra-individual covariates. Objects of class, \code{repeated},
contain a \code{response} object and possibly \code{tccov} and
\code{tvcov} objects.

In formula and functions, the key words, \code{times} can be used to
refer to the response times from the data object as a covariate,
\code{individuals} to the index for individuals as a factor covariate,
and \code{nesting} the index for nesting as a factor covariate. The
latter two only work for W&R notation.

The following methods are available for accessing the contents of such
data objects.

\code{as.data.frame}: places all of the variables in the data object
in one dataframe, extending time-constant covariates to the length of
the others unless the object has class, \code{tccov}. Binomial and
censored response variables have two columns, respectively `yes' and
`no' and response and censoring indicator, with the name given to the
response.

\code{as.matrix}: places all of the variables in the data object
in one matrix, extending time-constant covariates to the length of
the others unless the object has class, \code{tccov}. If any
covariates are factor variables (instead of the corresponding sets of
indicator variables), the matrix will be character instead of numeric.

\code{covariates}: extracts covariate matrices from a data object (for
formulae and functions, possibly for selected individuals. See
\code{\link[rmutil]{covariates.formulafn}}).

\code{covind}: gives the indexing of the response by individual (that
is, the nesting indicator for observations within individuals). It can
be used to expand time-constant covariates to the size of the repeated
measurements response.

\code{delta}: extracts the units of measurement vector and Jacobian of
any transformation of the response, possibly for selected individuals.
Note that, if the unit of measurement/Jacobian is available in the
\code{response} object, this is automatically included in the
calculation of the likelihood function in all library model functions.

\code{units}: prints the variable names and their description
and returns the latter.

\code{formula}: gives the formula used to create the time-constant
covariate matrix of a data object (for formulae and functions, see
\code{\link[rmutil]{formula.formulafn}}).

\code{names}: extracts the names of the response and/or covariates.

\code{nesting}: gives the coding variable(s) for individuals (same as
\code{covind}) and also for nesting within individuals if available,
possibly for selected individuals.

\code{nobs}: gives the number of observations per individual.

\code{plot}: plots the variables in the data object in various ways.
For \code{repeated} objects, \code{name} can be a response or a
time-varying covariate.

\code{print}: prints summary information about the variables in a data object.

\code{response}: extracts the response vector, possibly for selected
individuals. If there are censored observations, this is a two-column
matrix, with the censor indicator in the second column. For binomial
data, it is a two-column matrix with "positive" (y) and "negative"
(totals-y) frequencies.

\code{resptype}: extracts the type of each response.

\code{times}: extracts the times vector, possibly for selected
individuals.

\code{transform}: transforms variables. For example,
\code{transform(z, y=fcn1(y), times=fcn2(times))} where \code{fcn1}
and \code{fcn2} are transformation functions. When the response is
transformed, the Jacobian is automatically calculated. New response
variables and covariates can be created in this way, if the left hand
side is a new name (\code{ynew=fcn3(y)}), as well as replacing an old
variable with the transformed one. If the transformation reverses the
order of the responses, use its negative to keep the ordering and have
a positive Jacobian; for example, \code{ry=-1/y}. For \code{repeated}
objects, only the response and the times can be transformed.

\code{units}: prints the variable names and their units of measurement
and returns the latter.

\code{weights}: extracts the weight vector, possibly for selected
individuals.
}
\arguments{
\item{x,z}{A \code{response}, \code{tccov}, \code{tvcov}, or
\code{repeated} data object. For \code{covind} and \code{nobs}, this
may also be a model.}
%\item{y}{The function, when the response is to be transformed.}
\item{times}{The function, when the times are to be transformed.}
\item{names}{The names of the response variable(s) or covariate(s).}
\item{nind}{The numbers of individuals to be used. (For plotting,
cannot be used simultaneously with \code{ccov}.)}
%\item{expand}{For intra-class (time-constant) covariates, if TRUE,
%expand them to give one value per observation rather than one per
%individual. Only works with \code{repeated} objects when all
%individuals are requested (\code{nind=NULL}).}
\item{ccov}{For plotting: If a vector of values for the time-constant
covariates is supplied, only individuals having that set of values
will have profiles plotted. These values must be in the order in which
the covariates appear when the data object is printed. For factor
variables, the codes must be given. If the name of a covariate is
supplied, a set of graphs is plotted, one for each covariate value,
showing profiles of all individuals having that value. (The covariate
can have a maximum of six values.) Cannot be used simultaneously with
\code{nind}.} 
\item{nest}{For plotting: nesting category to plot.}
\item{add}{For plotting: add to previous plot.}
\item{nindmax}{For printing a \code{response}, \code{tvcov}, or
\code{repeated} object, if the number of individuals is greater than
\code{nindmax}, the range of numbers of observations per individual is
printed instead of the vector of numbers.}
\item{name,lty,pch,main,ylim,xlim,xlab,ylab}{ See base plot.}
\item{_data,units,object}{ TBD.}
\item{...}{ Arguments to other methods}

}
\value{
These methods extract information stored in \code{response},
\code{tccov}, \code{tvcov}, and \code{repeated} data objects created
respectively by \code{\link[rmutil]{restovec}}, \code{\link[rmutil]{tcctomat}},
\code{\link[rmutil]{tvctomat}}, and \code{\link[rmutil]{rmna}}.

Note that if a vector of binomial totals or a censoring indicator is
present, this is extract as the second column of the matrix produced
by the \code{response} method.
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{restovec}}, \code{\link[rmutil]{rmna}},
\code{\link[rmutil]{tcctomat}}, \code{\link[rmutil]{tvctomat}}.
}
\examples{
# set up some data and create the objects
#
y <- matrix(rnorm(20),ncol=5)
tt <- c(1,3,6,10,15)
print(resp <- restovec(y, times=tt, units="m", type="duration"))
x <- c(0,0,1,1)
x2 <- as.factor(c("a","b","a","b"))
tcc <- tcctomat(data.frame(x=x,x2=x2))
z <- matrix(rpois(20,5),ncol=5)
tvc <- tvctomat(z)
print(reps <- rmna(resp, tvcov=tvc, ccov=tcc))
#
plot(resp)
plot(reps)
plot(reps, nind=1:2)
plot(reps, ccov=c(0,1))
plot(reps, ccov="x2")
plot(reps, name="z", nind=3:4, pch=1:2)
plot(reps, name="z", ccov="x2")
#
response(resp)
response(transform(resp, y=1/y))
response(reps)
response(reps, nind=2:3)
response(transform(reps,y=1/y))
#
times(resp)
times(transform(resp,times=times-6))
times(reps)
#
delta(resp)
delta(reps)
delta(transform(reps,y=1/y))
delta(transform(reps,y=1/y), nind=3)
#
nobs(resp)
nobs(tvc)
nobs(reps)
#
units(resp)
units(reps)
#
resptype(resp)
resptype(reps)
#
weights(resp)
weights(reps)
#
covariates(tcc)
covariates(tcc, nind=2:3)
covariates(tvc)
covariates(tvc, nind=3)
covariates(reps)
covariates(reps, nind=3)
covariates(reps,names="x")
covariates(reps,names="z")
#
names(tcc)
names(tvc)
names(reps)
#
nesting(resp)
nesting(reps)
#
# because individuals are the only nesting, this is the same as
covind(resp)
covind(reps)
#
as.data.frame(resp)
as.data.frame(tcc)
as.data.frame(tvc)
as.data.frame(reps)
#
# use in glm
rm(y,x,z)
glm(y~x+z, data=as.data.frame(reps))
}
\keyword{manip}
