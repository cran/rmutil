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
\item{weights}{A vector of frequencies or weights. Only applicable
when \code{response} is a vector of independent observations.}
\item{delta}{The unit of measurement (if not equal to unity), a scalar
or a vector.}
\item{type}{The type of measurement: nominal, ordinal, discrete,
duration, continuous, or unknown. This becomes an attribute of the
\code{response} object.}
}
\description{
\code{restovec} transforms a list of vectors of event histories or of
two or more column matrices with times, response values, and possibly
binomial totals, nesting categories, censor indicators, and/or units of
measurement, for each individual, or a matrix or dataframe of response
values into an object of class, response.

Such objects can be printed and plotted. NAs are removed with
\code{rmna} (where necessary, incoordination with the appropriate
covariates).

Methods are available for extracting the response variable, the times,
the nesting variable, and the weights.
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
\code{\link{tvctomat}}.
}
\author{J.K. Lindsey}
\examples{
y <- matrix(rnorm(20),ncol=5)
restovec(y)
times <- c(1,3,6,10,15)
restovec(y,times=times)
}
