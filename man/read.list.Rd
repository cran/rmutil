\name{read.list}
\title{Read a list of matrices from a file for unbalanced repeated
measurements.}
\usage{
read.list(file="", skip=0, nlines=2, order=NULL)
}
\arguments{
\item{file}{Name of the file to read}
\item{skip}{Number of lines to skip at the beginning of the file}
\item{nlines}{Number of lines per matrix}
\item{order}{Order in which the lines are to be used as columns of the
matrix. If NULL, they are placed in the order read.}
}
\description{
\code{read.list} reads sets of lines of data from a file and creates a
list of matrices. Different sets of lines may be have different lengths.
}
\value{
The list of matrices, each with nlines columns, is returned.
}
\author{J.K. Lindsey}
\seealso{
\link{read.surv}}, \code{\link{restovec}}, \code{\link{rmna}}.
}
\examples{
y <- read.list("test.dat")
}
