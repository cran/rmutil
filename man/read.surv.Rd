\name{read.surv}
\title{Read a list of matrices from a file for use with kalsurv.}
\usage{
read.list(file="", skip=0, nlines=1, cumulative=T, all=T)
}
\arguments{
\item{file}{Name of the file to read}
\item{skip}{Number of lines to skip at the beginning of the file}
\item{nlines}{Number of lines in each series of duration times}
\item{cumulative}{If TRUE, the times are cumulative and differences
are taken to obtain times between events. Otherwise, the times are
used unchanged.}
\item{all}{If TRUE, all times have accompanying censor indicators;
otherwise, only the last one does.}
}
\description{
\code{read.surv} reads sets of lines of data from a file. Each set may
contain a series of duration times followed by a censor indicator for
the last value (all=FALSE) or a series of pairs of times followed by
their censor indicators (all=TRUE).
}
\value{
A list containing a list of vectors containing the series of times
and a vector of censor indicators for the last time of each series is
returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{read.list}}, \code{\link{restovec}}, \code{\link{rmna}}.
}
\examples{
y <- read.surv("test.dat")
}
