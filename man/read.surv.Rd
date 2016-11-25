\name{read.surv}
\title{Read a List of Matrices from a File for Repeated Times to Events}
\alias{read.surv}
\usage{
read.surv(file="", skip=0, nlines=1, cumulative=TRUE, all=TRUE)
}
\description{
\code{read.surv} reads sets of lines of data from a file. Each set may
contain a series of duration times followed by a censor indicator for
the last value (\code{all=FALSE}) or a series of pairs of times followed by
their censor indicators (\code{all=TRUE}).
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
\value{
A list containing a list of vectors with the series of times and a
vector of censor indicators for the last time of each series is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{lvna}}, \code{\link[rmutil]{read.list}},
\code{\link[rmutil]{read.rep}}, \code{\link[rmutil]{restovec}},
\code{\link[rmutil]{rmna}}
}
\examples{
\dontrun{y <- read.surv("test.dat")}
}
\keyword{file}

