\name{plot.residuals}
\title{Plot Residuals}
\usage{
plot.residuals(z, x=NULL, subset=NULL, ccov=NULL, nind=NULL,
	recursive=TRUE, pch=20, ylab="Residual", xlab=NULL,
	main=NULL, ...)
}
\alias{plot.residuals}
\arguments{
\item{z}{An object of class recursive, from \code{\link{carma}},
\code{\link{gar}}, \code{\link{kalcount}}, \code{\link{kalseries}},
\code{\link{kalsurv}}, or \code{\link{nbkal}}.}
\item{x}{Vector of of values for the x-axis. If missing, time is
used. It can also be specified by the strings "response" or "fitted".}
\item{subset}{A logical vector defining which observations are to be
used.}
\item{ccov}{If the index number of a time-constant covariate is
supplied, separate plots are made for each distinct value of that
covariate.}
\item{nind}{Observation number(s) of individual(s) to be plotted.}
\item{recursive}{If TRUE, plot recursive residuals, otherwise ordinary
residuals.}
\item{others}{Plotting control options.}
}
\value{
\code{plot.residuals} is used for plotting residuals from models
obtained from Kalman fitting for given subsets of the data.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{carma}}, \code{\link{gar}}, \code{\link{kalcount}},
\code{\link{kalseries}}, \code{\link{kalsurv}}, \code{\link{nbkal}}
\code{\link{plot.iprofile}}, \code{\link{plot.profile}}.
}
\examples{
times <- rep(1:20,2)
dose <- c(rep(2,20),rep(5,20))
mu <- function(p) exp(p[1]-p[3])*(dose/(exp(p[1])-exp(p[2]))*
	(exp(-exp(p[2])*times)-exp(-exp(p[1])*times)))
shape <- function(p) exp(p[1]-p[2])*times*dose*exp(-exp(p[1])*times)
conc <- matrix(rgamma(40,2,mu(log(c(1,0.3,0.2)))),ncol=20,byrow=T)
conc[,2:20] <- conc[,2:20]+0.5*(conc[,1:19]-matrix(mu(log(c(1,0.3,0.2))),
	ncol=20,byrow=T)[,1:19])
conc <- ifelse(conc>0,conc,0.01)
z <- gar(conc, dist="gamma", times=1:20, mu=mu, shape=shape,
	preg=log(c(1,0.4,0.1)), pdepend=0.1, pshape=log(c(1,0.2)))
plot.residuals(z, subset=1:20, main="Dose 1")
plot.residuals(z, x="fitted", subset=1:20, main="Dose 1")
plot.residuals(z, x="response", subset=1:20, main="Dose 1")
}
