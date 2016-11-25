\name{read.rep}
\title{Read a Rectangular Data Set from a File to Create a repeated Object}
\alias{read.rep}
\usage{
read.rep(file, header=TRUE, skip=0, sep = "", na.strings="NA",
	response, id=NULL, times=NULL, censor=NULL, totals=NULL,
	weights=NULL, nest=NULL, delta=NULL, coordinates=NULL,
	type=NULL, ccov=NULL, tvcov=NULL, na.rm=TRUE)
}
\description{
\code{dftorep} forms an object of class, \code{repeated}, from data
read from a file with the option of removing any observations where
response and covariate values have NAs. For repeated measurements,
observations on the same individual must be together in the file. A
number of validity checks are performed on the data.

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
\item{file}{A file name from which to read the data with variables as
columns and observations as rows.}

\item{header}{A logical value indicating whether the file contains the
names of the variables as the line before the first row of data.}

\item{skip}{The number of lines of the file to skip before beginning
to read data.}

\item{sep}{The field separator character. Values on each line of the
file are separated by this character.}

\item{na.strings}{A vector of strings defining what values are to be
assigned NA.}

\item{response}{A character vector giving the column name(s) of the
dataframe for the response variable(s).}

\item{id}{A character vector giving the column name of the
dataframe for the identification numbers of the individuals. If the
numbers are not consecutive integers, a warning is given.

If NULL, one observation per individual is assumed if \code{times} is
also NULL, other time series is assumed.}

\item{times}{An optional character vector giving the column name of the
dataframe for the times vector.}

\item{censor}{An optional character vector giving the column name(s)
of the dataframe for the censor indicator(s). This must be the same
length as \code{response}. Responses without censor indicator can have
a column either of all NAs or all 1s.}

\item{totals}{An optional character vector giving the column name(s)
of the dataframe for the totals for binomial data. This must be the same
length as \code{response}. Responses without censor indicator can have
a column all NAs.}

\item{weights}{An optional character vector giving the column name of
the dataframe for the weights vector.}

\item{nest}{An optional character vector giving the column name of
the dataframe for the nesting vector within individuals.

This is the second level of nesting for repeated measurements, with
the individual being the first level. Values for an individual must be
consecutive increasing integers.}

\item{delta}{An optional character vector giving the column name(s)
of the dataframe for the units of measurement/Jacobian(s) of the
response(s). This must be the same length as \code{response}.
Responses without units of measurement/Jacobian can have a column all
NAs.

If all response variables have the same unit of measurement, this can
be that one number. If each response variable has the same unit of
measurement for all its values, this can be a numeric vector of length
the number of response variables.}

\item{coordinates}{An optional character vector giving the two or
three column name(s) of the dataframe for the spatial coordinates.}

\item{type}{An optional character vector giving the types of response
variables: nominal, ordinal, discrete, duration, continuous,
multivariate, or unknown.}

\item{ccov}{An optional character vector giving the column names of
the dataframe for the time-constant or inter-individual covariates.
For repeated measurements, if the value is not constant for all
observations on an individual, an error is produced.}

\item{tvcov}{An optional character vector giving the column names of
the dataframe for the time-varying or intra-individual covariates.}

\item{na.rm}{If TRUE, observations with NAs in any variables selected
are removed in the object returned. Otherwise, the corresponding
indicator variable is returned in a slot in the object.}
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
\code{\link[rmutil]{dftorep}}, \code{\link[rmutil]{lvna}},
\code{\link[rmutil]{read.list}}, \code{\link[rmutil]{restovec}},
\code{\link[rmutil]{rmna}}, \code{\link[rmutil]{tcctomat}},
\code{\link[rmutil]{tvctomat}}
}
\examples{
\dontrun{read.rep("test.dat", resp=c("y1","y2"), times="tt", id="id",}
\dontrun{	totals=c("tot1","tot2"), tvcov="x",ccov="x2")}
}
\keyword{file}
