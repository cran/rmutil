\name{fnenvir}
\title{Check Covariates and Parameters of a Function}
\alias{fnenvir}
\alias{fnenvir.data.frame}
\alias{fnenvir.default}
\alias{fnenvir.repeated}
\alias{fnenvir.tccov}
\alias{fnenvir.tvcov}
\usage{
fnenvir(.z, ...)
\method{fnenvir}{default}(.z, .envir=parent.frame(), .name=NULL, .expand=TRUE,
	.response=FALSE, ...)
}
\description{
\code{fnenvir} finds the covariates and parameters in a function and
can modify it so that the covariates used in it are found in the data
object specified by \code{.envir}.

If the data object has class, \code{repeated}, the key word
\code{times} in a function will use the response times from the data
object as a covariate.
}
\arguments{
\item{.z}{A function.}
\item{.envir}{The environment or data object of class,
\code{repeated}, \code{tccov}, or \code{tvcov}, in which the function
is to be interpreted.}
\item{.name}{Character string giving the name of the data object
specified by \code{.envir}. Ignored unless the latter is such an
object and only necessary when \code{fnenvir} is called within other
functions.}
\item{.expand}{If TRUE, expand functions with only time-constant
covariates to return one value per observation instead of one value
per individual. Ignored unless \code{.envir} is an object of class,
\code{repeated}.}
\item{.response}{If TRUE, any response variable can be used in the
function. If FALSE, checks are made that the response is not also used
as a covariate.}
\item{...}{Arguments passed to other functions.}
}
\value{
The (modified) function, of class \code{formulafn}, is returned with its
attributes giving the (new) model function, the covariate names, and
the parameter names.}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{FormulaMethods}},\code{\link[rmutil]{covariates}},
\code{\link[rmutil]{finterp}}, \code{\link[rmutil]{model}},
\code{\link[rmutil]{parameters}}
}
\examples{
fn <- function(p) a+b*x
fnenvir(fn)
fn <- function(p) a+p*x
fnenvir(fn)
x <- 1:4
fnenvir(fn)
fn <- function(p) p[1]+exp(p[2]*x)
fnenvir(fn)
#
y <- matrix(rnorm(20),ncol=5)
y[3,3] <- y[2,2] <- NA
resp <- restovec(y)
xx <- tcctomat(x)
z1 <- matrix(rnorm(20),ncol=5)
z2 <- matrix(rnorm(20),ncol=5)
z3 <- matrix(rnorm(20),ncol=5)
zz <- tvctomat(z1)
zz <- tvctomat(z2,old=zz)
reps <- rmna(resp, ccov=xx, tvcov=zz)
rm(y, x, z1, z2)
#
# repeated objects
func1 <- function(p) p[1]+p[2]*x+p[3]*z2
print(fn1 <- fnenvir(func1, .envir=reps))
fn1(2:4)
#
# time-constant covariates
func2 <- function(p) p[1]+p[2]*x
print(fn2 <- fnenvir(func2, .envir=reps))
fn2(2:3)
print(fn2a <- fnenvir(func2, .envir=xx))
fn2a(2:3)
#
# time-varying covariates
func3 <- function(p) p[1]+p[2]*z1+p[3]*z2
print(fn3 <- fnenvir(func3, .envir=reps))
fn3(2:4)
print(fn3a <- fnenvir(func3, .envir=zz))
fn3a(2:4)
# including times
func3b <- function(p) p[1]+p[2]*z1+p[3]*z2+p[4]*times
print(fn3b <- fnenvir(func3b, .envir=reps))
fn3b(2:5)
#
# with typing error and a variable not in the data object
func4 <- function(p) p[1]+p2[2]*z1+p[3]*z2+p[4]*z3
print(fn4 <- fnenvir(func4, .envir=reps))
#
# first-order one-compartment model
# data objects for formulae
dose <- c(2,5)
dd <- tcctomat(dose)
times <- matrix(rep(1:20,2), nrow=2, byrow=TRUE)
tt <- tvctomat(times)
# vector covariates for functions
dose <- c(rep(2,20),rep(5,20))
times <- rep(1:20,2)
# functions
mu <- function(p) {
	absorption <- exp(p[1])
	elimination <- exp(p[2])
	absorption*exp(-p[3])*dose/(absorption-elimination)*
		(exp(-elimination*times)-exp(-absorption*times))}
shape <- function(p) exp(p[1]-p[2])*times*dose*exp(-exp(p[1])*times)
# response
conc <- matrix(rgamma(40,shape(log(c(0.1,0.4))),
	scale=mu(log(c(1,0.3,0.2))))/shape(log(c(0.1,0.4))),ncol=20,byrow=TRUE)
conc[,2:20] <- conc[,2:20]+0.5*(conc[,1:19]-matrix(mu(log(c(1,0.3,0.2))),
	ncol=20,byrow=TRUE)[,1:19])
conc <- restovec(ifelse(conc>0,conc,0.01))
reps <- rmna(conc, ccov=dd, tvcov=tt)
#
print(fn5 <- fnenvir(mu,.envir=reps))
fn5(c(0,-1.2,-1.6))
}
\keyword{programming}
