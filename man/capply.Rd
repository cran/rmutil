\name{capply}
\alias{capply}
\title{A Fast Simplified Version of \code{tapply}}
\description{
a fast simplified version of \code{tapply}
}
\usage{
capply(x, index, fcn=sum)
}
\arguments{
\item{x}{x}
\item{index}{index}
\item{fcn}{default sum}
}
\details{
a fast simplified version of tapply
}
\value{
Returns \code{ans} where \code{for(i in split(x,index))ans <- c(ans,fcn(i))}.
}
