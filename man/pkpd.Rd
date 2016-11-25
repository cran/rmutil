\name{pkpd}
\title{Pharmacokinetic Compartment Models}
\usage{
mu1.0o1c(p, times, dose=1, end=0.5)
mu1.1o1c(p, times, dose=1)
mu1.1o2c(p, times, dose=1)
mu1.1o2cl(p, times, dose=1)
mu1.1o2cc(p, times, dose=1)
mu2.0o1c(p, times, dose=1, ind, end=0.5)
mu2.0o2c1(p, times, dose=1, ind, end=0.5)
mu2.0o2c2(p, times, dose=1, ind, end=0.5)
mu2.1o1c(p, times, dose=1, ind)
mu2.0o1cfp(p, times, dose=1, ind, end=0.5)
mu2.0o2c1fp(p, times, dose=1, ind, end=0.5)
mu2.0o2c2fp(p, times, dose=1, ind, end=0.5)
mu2.1o1cfp(p, times, dose=1, ind)
}
\alias{pkpd}
\alias{mu1.0o1c}
\alias{mu1.1o1c}
\alias{mu1.1o2c}
\alias{mu1.1o2cl}
\alias{mu1.1o2cc}
\alias{mu2.0o1c}
\alias{mu2.0o2c1}
\alias{mu2.0o2c2}
\alias{mu2.1o1c}
\alias{mu2.0o1cfp}
\alias{mu2.0o2c1fp}
\alias{mu2.0o2c2fp}
\alias{mu2.1o1cfp}
\description{
Mean functions for use in fitting pharmacokineticcompartment models models.

\code{mu1.0o1c}: open zero-order one-compartment model

\code{mu1.1o1c}: open first-order one-compartment model

\code{mu1.1o2c}: open first-order two-compartment model (ordered)

\code{mu1.1o2cl}: open first-order two-compartment model (ordered,
	absorption and transfer equal)

\code{mu1.1o2cc}: open first-order two-compartment model (circular)

Simultaneous models for parent drug and metabolite:

\code{mu2.0o1c}: zero-order one-compartment model

\code{mu2.0o2c1}: zero-order two-compartment for parent,
one-compartment for metabolite, model

\code{mu2.0o2c2}: zero-order two-compartment model for both parent and
metabolite

\code{mu2.1o1c}: first-order one-compartment model

\code{mu2.0o1cfp}: zero-order one-compartment first-pass model

\code{mu2.0o2c1fp}: zero-order two-compartment for parent,
one-compartment for metabolite, model with first-pass

\code{mu2.0o2c2fp}: zero-order two-compartment model for both parent and
metabolite with first-pass

\code{mu2.1o1cfp}: first-order one-compartment first-pass model
}
\arguments{
\item{p}{Vector of parameters. See the source file for details.}
\item{times}{Vector of times.}
\item{dose}{Vector of dose levels.}
\item{ind}{Indicator whether parent drug or metabolite.}
\item{end}{Time infusion ends.}
}
\value{
The profile of mean concentrations for the given times and doses is
returned.
}
\author{J.K. Lindsey}
\examples{
\dontrun{
library(repeated)
times <- rep(1:20,2)
dose <- c(rep(2,20),rep(5,20))
# set up a mean function for gar based on mu1.1o1c:
mu <- function(p) {
	ka <- exp(p[2])
	ke <- exp(p[3])
	exp(p[2]-p[1])/(ka-ke)*(exp(-ke*times)-exp(-ka*times))}
conc <- matrix(rgamma(40,2,scale=mu(log(c(1,0.3,0.2)))/2),ncol=20,byrow=TRUE)
conc[,2:20] <- conc[,2:20]+0.5*(conc[,1:19]-matrix(mu(log(c(1,0.3,0.2))),
	ncol=20,byrow=TRUE)[,1:19])
conc <- ifelse(conc>0,conc,0.01)
gar(conc, dist="gamma", times=1:20, mu=mu, preg=log(c(1,0.4,0.1)),
	pdepend=0.1, pshape=1)
# changing variance
shape <- mu
gar(conc, dist="gamma", times=1:20, mu=mu, preg=log(c(0.5,0.4,0.1)),
	pdep=0.1, shape=shape, pshape=log(c(0.5,0.4,0.1)))
}
}
\keyword{models}
