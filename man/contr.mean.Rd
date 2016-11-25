\name{contrast}
\alias{contr.mean}
\title{Contrast Matrix for Constraints about the Mean}
\description{
Return a matrix of contrasts for constraints about the mean.
}
\usage{
contr.mean(n, contrasts = TRUE)
}
\arguments{
\item{n}{A vector of levels for a factor or the number of levels.}
\item{contrasts}{A logical value indicating whether or not contrasts
should be computed.}
}
\details{
This function corrects \code{\link{contr.sum}} to display labels properly.
}
\value{
A matrix of computed contrasts with \code{n} rows and \code{k}
columns, with \code{k=n-1} if \code{contrasts} is \code{TRUE} and
\code{k=n} if \code{contrasts} is \code{FALSE}. The columns of the
resulting matrices contain contrasts which can be used for coding a
factor with \code{n} levels.
}
\seealso{
\code{\link{contrasts}}, \code{\link{C}},
and \code{\link{contr.sum}}.
}
\examples{
oldop <- options(contrasts=c("contr.sum","contra.poly"))
y <- rnorm(30)
x <- gl(3,10,labels=c("First","Second","Third"))
glm(y~x)
options(contrasts=c("contr.mean","contra.poly"))
x <- gl(3,10,labels=c("First","Second","Third"))
glm(y~x)
options(oldop)
}
\keyword{design}
\keyword{regression}
\keyword{array}
