\name{iprofile}
\title{Produce Individual Time Profiles for Plotting}
\usage{
plot(iprofile(z, plotsd=FALSE), nind=1, observed=TRUE, intensity=F,
	add=FALSE, lty=NULL, pch=NULL, ylab=NULL, xlab=NULL,
	main=NULL, ylim=NULL, xlim=NULL, ...)
}
\arguments{
\item{z}{An object of class recursive, from \code{\link{carma}},
\code{\link{elliptic}}, \code{\link{gar}}, \code{\link{kalcount}},
\code{\link{kalseries}}, \code{\link{kalsurv}}, or \code{\link{nbkal}}.}
\item{plotsd}{If TRUE, plots standard deviations around profile
(\code{\link{carma}} and \code{\link{elliptic}} only).}
\item{nind}{Observation number(s) of individual(s) to be plotted.}
\item{observed}{If TRUE, plots observed responses.}
\item{intensity}{If z has class, \code{kalsurv}, and this is TRUE, the
intensity is plotted instead of the time between events.}
\item{add}{If TRUE, the graph is added to an existing plot.}
\item{others}{Plotting control options.}
}
\value{
\code{iprofile} is used for plotting individual profiles over time
for models obtained from Kalman fitting. See \code{\link{profile}} for
plotting marginal profiles.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{carma}}, \code{\link{elliptic}}, \code{\link{gar}},
\code{\link{kalcount}}, \code{\link{kalseries}},
\code{\link{kalsurv}}, \code{\link{nbkal}} \code{\link{profile}}
\code{\link{plot.residuals}}.
}
\examples{
library(repeated)
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
plot(iprofile(z), nind=1:2, pch=c(1,20), lty=3:4)
plot(profile(z), nind=1:2, lty=1:2, add=T)
}
