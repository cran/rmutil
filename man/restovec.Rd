\name{restovec}
\title{Create a response Object}
\alias{restovec}
\usage{
restovec(response=NULL, times=NULL, nest=NULL, coordinates=NULL,
	censor=NULL, totals=NULL, weights=NULL, delta=NULL,
	type=NULL, names=NULL, units=NULL, oldresponse=NULL,
	description=NULL)
}
\description{
\code{restovec} can produce an object of class, \code{response}, from
a vector of (1) independent univariate responses or (2) a single time
series.

It can produce such an object from repeated measurements in the form
of (1) a list of vectors of event histories, (2) a list of two or more
column matrices with times, response values, and and other information
or (3) a matrix or dataframe of response values. The first two are for
unbalanced data and the third for balanced data.

Multivariate responses can be supplied as (1) a three-dimensional
array of balanced repeated measurements, (2) lists of matrices for
unbalanced repeated measurements, or (3) a matrix with either (a)
several time series or (b) single observations per individual on
several variables.

In formula and functions, the key words, \code{times} can be used to
refer to the response times from the data object as a covariate,
\code{individuals} to the index for individuals as a factor covariate,
and \code{nesting} the index for nesting as a factor covariate. The
latter two only work for W&R notation.

NAs can be detected with \code{\link[rmutil]{lvna}} or removed with
\code{\link[rmutil]{rmna}} (where necessary, in coordination with the
appropriate covariates) to create a \code{repeated} object.

\code{response} objects can be printed and plotted. Methods are
available for extracting the response, the numbers of observations per
individual, the times, the weights, the units of measurement/Jacobian,
and the nesting variable: \code{\link[rmutil]{response}},
\code{\link[rmutil]{nobs}}, \code{\link[rmutil]{times}},
\code{\link[rmutil]{weights}}, \code{\link[rmutil]{delta}}, and
\code{\link[rmutil]{nesting}}.

The response and or the times may be \link[rmutil]{transform}ed using
\code{transform(z, newy=fcn1(y), times=fcn2(times))} where
\code{fcn1} and \code{fcn2} are transformations and \code{y} is the
name of a response variable. When the response is
transformed, the Jacobian is automatically calculated. Note that, if
the unit of precision/Jacobian (\code{\link[rmutil]{delta}}) is available in
the \code{response} object, this is automatically included in the
calculation of the likelihood function in all library model functions.
}
\arguments{
\item{response}{For (1) independent univariate responses with one
observation per individual or (2) a single time series, one vector may
be supplied (in the latter case, the times must be given even if
equally spaced).

Univariate repeated measurements responses can be given (1) if
balanced, as a matrix or dataframe of response values with dimensions:
number of individuals by number of responses/individual, (2) a list of
vectors of event histories, or (3) a list of one or more column
matrices, for each individual, with response values in the first
column and times in the second (if there are no times, set
\code{times} to FALSE), possibly followed by columns with nesting
categories, binomial totals, censoring indicators, and/or units of
measurement.

Multivariate responses can be supplied as (1) a three-dimensional
array of balanced repeated measurements with dimensions: number of
individuals by number of responses/individual by number of variables,
(2) a list of matrices for unbalanced repeated measurements each with
dimensions: number of responses on that individual by number of
variables, plus a column for times if available (otherwise set
\code{times} to FALSE), or (3) a matrix with either (a) several time
series, having dimensions: length of time series by by number of times
series, or (b) single observations per individual on several variables
with dimensions: number of individuals by number of variables. In all
but case (1), \code{type} must be a character vector with length equal
to the number of responses. In case (2), where applicable, \code{censor},
\code{totals}, and \code{delta} must be supplied as lists of matrices
of the same size as for \code{response}, and \code{nest} and
\code{weights} as lists of vectors of lengths equal to the number of
observations on each individual.}

\item{times}{When \code{response} is a matrix or multivariate array,
these can be (1) a vector when the times are the same for all
individuals, possibly unequally-spaced, or (2) a matrix with
dimensions: number of individuals by number of responses/individual.
Not necessary if times are equally spaced, except if a vector
containing a single time series is supplied (if not given in this
case, it takes the responses to be independent, not a time
series). For clustered data with no time ordering, set to FALSE.}

\item{nest}{This is the second level of nesting, with the individual
being the first level. Values for an individual must be consecutive
increasing integers with all responses in the same cluster grouped
together. For example, with three clusters of four observations each,
the code would be 1,1,1,1,2,2,2,2,3,3,3,3.

When \code{response} is a matrix or multivariate array, this can be a
vector of length equal to the number of responses/individual
indicating which responses belong to which nesting category.

If \code{response} is a multivariate list, this must also be a list.

When \code{response} is a univariate list of unbalanced repeated
measurements, the nesting indicator may instead be included in that
list but must respect the same ordering as described above.}

\item{coordinates}{When \code{response} is a vector, a two-column
matrix giving the coordinates for spatial data.}

\item{censor}{When \code{response} is a matrix, this can be (1) a vector
of the same length as the number of individuals, containing a binary
indicator, with a one indicating that the last time period in the
series terminated with an event and zero that it was censored, or
(2) a matrix of the same size as \code{response}.

When \code{response} is a multivariate array, this can be (1) a matrix
with dimensions: number of individuals by number of responses, or (2)
an array of the same size as \code{response}. In the first case, for
each column corresponding to a duration response, it should contain a
binary indicator, with a one indicating that the last time period in
the series terminated with an event and zero that it was censored, and
NAs in columns not containing durations. In the second case, layers
not corresponding to duration responses should contain NAs.

If \code{response} is a multivariate list, this must also be a list.

For event history data, even with no censoring, an appropriate vector
of ones must be supplied.

When \code{response} is a univariate list of unbalanced repeated
measurements, the censoring indicator may instead be included in that
list.}

\item{totals}{If the response is a matrix of binomial counts, this can
be (1) a corresponding vector (one total per individual) or (2) a
matrix of totals.

When \code{response} is a multivariate array, this can be (1) a matrix
with dimensions: number of individuals by number of responses if all
binomial responses for an individual have the same total, or (2) an
array of the same size as \code{response}. In the first case, for each
column corresponding to a binomial response, it should contain the
corresponding totals, with NAs in columns not containing binomial. In
the second case, layers not corresponding to binomial responses should
contain NAs.

If \code{response} is a multivariate list, this must also be a list.

When \code{response} is a univariate list of unbalanced repeated
measurements, the totals may instead be included in that list.}

\item{weights}{A vector, matrix, array, or list of vectors of
frequencies or weights, with one value per \code{response}. In other
words, a multivariate response has only one corresponding weight value.}

\item{delta}{For continuous measurements, the unit of precision (if
not equal to unity) for each response: a scalar, vector, matrix,
array, or list of the same dimensions as \code{response}. For example,
if responses have two decimal places (\code{12.34}), \code{delta=0.01}.
If the response has been transformed, this should be multiplied by the
numerical values of the Jacobian. When the \code{transform} method is
applied to the \code{response} object, this is automatically updated.}

\item{type}{The type(s) of observations: nominal, ordinal, discrete,
duration, continuous, or unknown. If not specified otherwise, those
responses with \code{delta} and no \code{censor} are assumed to be
\code{continuous}, those with \code{censor} indicator are assumed to
be \code{duration}, those with \code{totals} are assumed to be
\code{nominal}, and all others \code{unknown}.}

\item{names}{Optional name(s) of the response variable(s).}
\item{units}{Optional character vector giving units of measurement of
response(s).}
\item{oldresponse}{An existing \code{response} object to which the new
data are to be added.}
\item{description}{An optional named list of character vectors with
names of some or all response variables containing their descriptions.}
}
\value{
Returns an object of class, \code{response}, containing a vector with the
responses (\code{z$y}), a corresponding vector of times (\code{z$times}) if
applicable, a vector giving the number of observations per individual
(\code{z$nobs}, set to a scalar 1 if observations are independent),
type  (\code{z$delta}), and possibly binomial totals (\code{z$n}),
nesting (clustering, \code{z$nest}), censoring (\code{z$censor}),
weights (\code{z$wt}),  unit of precision/Jacobian (\code{z$delta}),
units of measurement (\code{z$units}), and description
(\code{z$description}) information.
}
\seealso{
\code{\link[rmutil]{DataMethods}}, \code{\link[rmutil]{covind}},
\code{\link[rmutil]{delta}}, \code{\link[rmutil]{description}},
\code{\link[rmutil]{lvna}}, \code{\link[rmutil]{names}},
\code{\link[rmutil]{nesting}}, \code{\link[rmutil]{nobs}},
\code{\link[rmutil]{read.list}}, \code{\link[rmutil]{read.surv}},
\code{\link[rmutil]{response}}, \code{\link[rmutil]{resptype}},
\code{\link[rmutil]{rmna}}, \code{\link[rmutil]{tcctomat}},
\code{\link[rmutil]{times}}, \code{\link[rmutil]{transform}},
\code{\link[rmutil]{tvctomat}}, \code{\link[rmutil]{units}},
\code{\link[rmutil]{weights}}
}
\author{J.K. Lindsey}
\examples{
#
#continuous response
y <- matrix(rnorm(20),ncol=5)
# times assumed to be 1:5
restovec(y, units="m")
#unequally-spaced times
tt <- c(1,3,6,10,15)
print(resp <- restovec(y, times=tt, units="m",
	description=list(y="Response measured in metres")))
response(resp)
response(resp, nind=2:3)
response(transform(resp, y=1/y))
transform(resp, y=1/y, units="1/m")
units(resp)
description(resp)
times(resp)
times(transform(resp, times=times-6))
nobs(resp)
weights(resp)
nesting(resp)
# because individuals are the only nesting, this is the same as
covind(resp)
#
# binomial response
y <- matrix(rpois(20,5),ncol=5)
# responses summarized as relative frequencies
print(respb <- restovec(y, totals=y+matrix(rpois(20,5),ncol=5), times=tt))
response(respb)
#
# censored data
y <- matrix(rweibull(20,2,5),ncol=5)
print(respc <- restovec(y, censor=matrix(rbinom(20,1,0.9),ncol=5), times=tt))
# if there is no censoring, censor indicator is not printed
response(respc)
# nesting clustered within individuals
nest <- c(1,1,2,2,2)
print(respn <- restovec(y, censor=matrix(rbinom(20,1,0.9),ncol=5),
	times=tt,nest=nest))
response(respn)
times(respn)
nesting(respn)
#
# multivariate response
restovec(y, censor=matrix(rbinom(20,1,0.9),ncol=5),
	units=c("m","days","l","cm","mon"),
	type=c("continuous","duration","continuous","continuous","duration"),
	description=list(y1="First continuous variable",
		y2="First duration variable",y3="Second continuous variable",
		y4="Third continuous variable",y5="Second duration variable"))
restovec(y, censor=matrix(rbinom(20,1,0.9),ncol=5), 
	names=c("a","b","c","d","e"), units=c("m","days","l","cm","mon"),
	type=c("continuous","duration","continuous","continuous","duration"),
	description=list(a="First continuous variable",
		b="First duration variable",c="Second continuous variable",
		d="Third continuous variable",e="Second duration variable"))
}
\keyword{manip}
