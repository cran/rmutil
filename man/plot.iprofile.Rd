\name{plot.iprofile}
\title{Plot Individual Time Profiles}
\usage{
plot.iprofile(z, nind=1, obs=TRUE, add=FALSE, plotsd=FALSE, lty=NULL,
	pch=NULL, ylab="Recursive fitted value", xlab="Time",
	main=NULL, ylim=NULL, xlim=NULL, ...)
}
\alias{plot.iprofile}
\arguments{
\item{z}{An object of class recursive, from \code{\link{carma}},
\code{\link{gar}}, \code{\link{kalcount}}, \code{\link{kalseries}},
\code{\link{kalsurv}}, or \code{\link{nbkal}}.}
\item{nind}{Observation number(s) of individual(s) to be plotted.}
\item{obs}{If TRUE, plots observed responses.}
\item{add}{If TRUE, the graph is added to an existing plot.}
\item{plotsd}{If TRUE, plots standard deviations around profile
(\code{\link{carma}} only).}
\item{others}{Plotting control options.}
}
\value{
\code{plot.iprofile} is used for plotting individual profiles over time
for models obtained from Kalman fitting. See
\code{\link{plot.profile}} for plotting marginal profiles.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{carma}}, \code{\link{gar}}, \code{\link{kalcount}},
\code{\link{kalseries}}, \code{\link{kalsurv}}, \code{\link{nbkal}}
\code{\link{plot.profile}} \code{\link{plot.residuals}}.
}
\examples{
times <- rep(1:20,2)
dose <- c(rep(2,20),rep(5,20))
mu <- function(p) exp(p[1]-p[3])*(dose/(exp(p[1])-exp(p[2]))*
	(exp(-exp(p[2])*times)-exp(-exp(p[1])*times)))
shape <- function(p) exp(p[1]-p[2])*times*dose*exp(-exp(p[1])*times)
conc <- matrix(rgamma(40,1,mu(log(c(1,0.3,0.2)))),ncol=20,byrow=T)
conc[,2:20] <- conc[,2:20]+0.5*(conc[,1:19]-matrix(mu(log(c(1,0.3,0.2))),
	ncol=20,byrow=T)[,1:19])
conc <- ifelse(conc>0,conc,0.01)
z <- gar(conc, dist="gamma", times=1:20, mu=mu, shape=shape,
	preg=log(c(1,0.4,0.1)), pdepend=0.5, pshape=log(c(1,0.2)))
# plot individual profiles and the average profile
plot.iprofile(z, nind=1:2, pch=c(1,20), lty=3:4)
plot.profile(z, nind=1:2, lty=1:2, add=T)
}
