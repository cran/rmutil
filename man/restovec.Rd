\name{restovec}
\title{Create a response object}
\usage{
restovec(response, times=NULL, nest=NULL, coordinates=NULL,
	censor=NULL, totals=NULL, weights=NULL, delta=NULL, type=NULL)
}
\arguments{
\item{response}{A matrix or dataframe of response values, a list of
vectors of event histories, or a list of one or more column matrices,
for each individual, with response values in the first column and
times in the second, possibly followed by columns with nesting
categories, binomial totals, censoring indicators, and/or units of
measurement. For independent responses with one observation per
individual or for a single time series, one vector may be supplied (in
the latter case, the times must be given even if equally spaced).}
\item{times}{When response is a matrix, a vector of possibly unequally
spaced times when they are the same for all individuals or a matrix of
times. Not necessary if equally spaced, except if a vector containing
a single time series is supplied (if not given in this case, it
considers the responses to be independent, not a time series). For
clustered data with no time ordering, set to FALSE.}
\item{nest}{When \code{response} is a matrix, a vector of length equal
to the number of responses per individual indicating which responses
belong to which nesting category. Values must be consecutive increasing
integers. This is the second level of nesting, with the individual
being the first level.}
\item{coordinates}{When \code{response} is a vector, a two-column
matrix giving the coordinates for spatial data.}
\item{censor}{A vector of the same length as the number of individuals
containing a binary indicator, with a one indicating that the last
time period in the series terminated with an event and zero that it
was censored, or, when response is matrix, a matrix of the same size.
For event history data, even with no censoring, an appropriate vector
of ones must be supplied.}
\item{totals}{If the response is a matrix of binomial counts, a
corresponding vector (one total or one per individual) or matrix of
totals.}
\item{weights}{A vector, matrix, or list of vectors of frequencies or
weights, of the same dimensions as the response.}
\item{delta}{The unit of measurement (if not equal to unity), a scalar
or a vector. If the response has been transformed, this should contain
the numerical values of the Jacobian. When the \code{transform} method
is applied to the response, this is automatically updated.}
\item{type}{The type of measurement: nominal, ordinal, discrete,
duration, continuous, or unknown. This becomes an attribute of the
\code{response} object.}
}
\description{
\code{restovec} produces an object of class, response, from a list of
vectors of event histories or of two or more column matrices with
times, response values, and possibly binomial totals, nesting
categories, censor indicators, and/or units of measurement, for each
individual, or a matrix or dataframe of response values.

Such objects can be printed and plotted. NAs are removed with
\code{rmna} (where necessary, incoordination with the appropriate
covariates).

Methods are available for extracting the response, the numbers of
observations per individual, the times, the weights, and the nesting
variable: \code{response}, \code{nobs}, \code{times}, \code{weights},
and \code{nesting}. The response and or the times may be transformed
using \code{transform(z, y=fcn1(y), times=fcn2(times))} where fcn1 and
fcn2 are transformations. When the response is transformed, the
Jacobian is automatically calculated.
}
\value{
Returns an object of class, response, containing a vector with the
responses (z$y), a corresponding vector of times (z$times) if
applicable, a vector giving the number of observations per individual
(z$nobs, set to a scalar 1 if observations are independent), and
possibly binomial totals (z$n), nesting (clustering, z$nest),
censoring (z$censor), weights (z$wt), and unit of measurement
information (z$delta).
}
\seealso{
\code{\link{read.list}}, \code{\link{rmna}}, \code{\link{tcctomat}},
\code{\link{transform}}, \code{\link{tvctomat}}.
}
\author{J.K. Lindsey}
\examples{
y <- matrix(rnorm(20),ncol=5)
# times assumed to be 1:5
restovec(y)
#unequally-spaced times
tt <- c(1,3,6,10,15)
print(resp <- restovec(y,times=tt))
response(resp)
response(resp, nind=2:3)
times(resp)
nobs(resp)
weights(resp)
nesting(resp)
# because individuals are the only nesting, this is the same as
covind(resp)
# binomial
y <- matrix(rpois(20,5),ncol=5)
# responses summarized as relative frequencies
print(respb <- restovec(y,totals=y+matrix(rpois(20,5),ncol=5),times=tt))
response(respb)
# censored data
y <- matrix(rweibull(20,2,5),ncol=5)
print(respc <- restovec(y,censor=matrix(rbinom(20,1,0.9),ncol=5),times=tt))
# if there is no censoring, censor indicator is not printed
response(respc)
# nesting clustered within individuals
nest <- c(1,1,2,2,2)
print(respn <- restovec(y,censor=matrix(rbinom(20,1,0.9),ncol=5),
	times=tt,nest=nest))
response(respn)
times(respn)
nesting(respn)
}
