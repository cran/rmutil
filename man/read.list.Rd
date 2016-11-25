\name{read.list}
\title{Read a List of Matrices from a File for Unbalanced Repeated
Measurements}
\alias{read.list}
\usage{
read.list(file="", skip=0, nlines=2, order=NULL)
}
\description{
\code{read.list} reads sets of lines of data from a file and creates a
list of matrices. Different sets of lines may be have different lengths.
}
\arguments{
\item{file}{Name of the file to read}
\item{skip}{Number of lines to skip at the beginning of the file}
\item{nlines}{Number of lines per matrix}
\item{order}{Order in which the lines are to be used as columns of the
matrix. If NULL, they are placed in the order read.}
}
\value{
The list of matrices, each with \code{nlines} columns, is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{lvna}}, \code{\link[rmutil]{read.rep}},
\code{\link[rmutil]{read.surv}}, \code{\link[rmutil]{restovec}},
\code{\link[rmutil]{rmna}}, \code{\link[rmutil]{tvctomat}}
}
\examples{
\dontrun{y <- read.list("test.dat")}
}
\keyword{file}

