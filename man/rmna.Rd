\name{rmna}
\title{Create a repeated object, removing NAs}
\usage{
rmna(response, tvcov=NULL, ccov=NULL)
}
\alias{rmna}
\alias{print.repeated}
\alias{plot.repeated}
\arguments{
\item{response}{An object of class, response (created by
\code{\link{restovec}}), containing the response variable information.}
\item{tvcov}{An object of class, tvcov (created by \code{\link{tvctomat}}),
containing the time-varying covariate information.}
\item{tccov}{An object of class, tccov (created by \code{\link{tcctomat}}),
containing the time-constant covariate information.}
}
\description{
\code{rmna} forms an object of class, repeated, from a response object
and possibly time-varying covariate (tvcov), and time-constant
covariate (tccov) objects, removing any response and covariate values
that have NAs.

Such objects can be printed and plotted.
}
\value{
Returns an object of class, repeated, containing a list of the
response object (z$response, so that, for example, the response vector
is z$response$y; see \code{\link{restovec}}), and possibly the two
classes of covariate objects (z$ccov and z$tvcov).
}
\author{J.K. Lindsey}
\seealso{
\code{\link{carma}}, \code{\link{elliptic}}, \code{\link{gettvc}},
\code{\link{kalcount}}, \code{\link{kalseries}}, \code{\link{nbkal}},
\code{\link{read.list}}, \code{\link{restovec}}, \code{\link{tcctomat}},
\code{\link{tvctomat}}.
}
\examples{
y <- matrix(rnorm(20),ncol=5)
times <- c(1,3,6,10,15)
resp <- restovec(y,times=times)
x <- c(0,0,1,1)
tcc <- tcctomat(x)
z <- matrix(rpois(20,5),ncol=5)
tvc <- tvctomat(z)
rmna(resp, tvcov=tvc, ccov=tcc)
}
