\name{finterp}
\title{Formula interpreter}
\usage{
finterp(z, envir=sys.frame(sys.parent()), formula=FALSE,
	vector=TRUE, start=1, name=NULL, expand=TRUE)
}
\arguments{
\item{z}{A model formula beginning with ~, either in Wilkinson and
Rogers notation or containing unknown parameters. If it also has
class, repeated, tccov, or tvcov, \code{envir} must be an object of
one of these classes.}
\item{envir}{The environment in which the formula is to be
interpreted or a data object of class, repeated, tccov, or tvcov.}
\item{formula}{If TRUE and the formula is in Wilkinson and Rogers
notation, just returns the formula.}
\item{vector}{If FALSE and the formula contains unknown parameters,
the function returned has them as separate arguments, if TRUE, it has
one argument, the unknowns as a vector. Always true if \code{envir} is
a data object.}
\item{start}{The starting index value of the parameter vector in the
function returned.}
\item{name}{Character string giving the name of the data object
specified by \code{envir}. Ignored unless the latter is such an
object and only necessary when \code{finterp} is called within other
functions.}
\item{expand}{If TRUE, expand functions with only time-constant
covariates to return one value per observation instead of one value
per individual. Ignored unless \code{envir} is an object of class,
repeated.}
}
\description{
\code{finterp} translates a model formula into a function of the
unknown parameters or of a vector of them. Such language formulae can
either be in Wilkinson and Rogers notation or be expressions
containing both known (existing) covariates and unknown (not existing)
parameters. In the latter, factor variables cannot be used and
parameters must be scalars. The covariates in the formula are sought
in the environment or in the data object provided.

Note that, in parameter displays, formulae in Wilkinson and Rogers
notation use variable names whereas those with unknowns use the names
of these parameters, as given in the formulae, and that the meaning of
operators (*, /, :, etc.) is different in the two cases.
}
\value{
A function, of class formulafn, of the unknown parameters or of a
vector of them is returned. Its attributes give the formula supplied,
the model function produced, the covariate names, the parameter names,
and the range of values of the index of the parameter vector. If
\code{formula} is TRUE and a Wilkinson and Rogers formula was
supplied, it is returned instead of a function.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{fnenvir}}
}
\examples{
x1 <- rpois(20,2)
x2 <- rnorm(20)
#
# Wilkinson and Rogers formula with three parameters
fn1 <- finterp(~x1+x2)
fn1
fn1(rep(2,3))
# the same formula with unknowns
fn2 <- finterp(~b0+b1*x1+b2*x2)
fn2
fn2(rep(2,3))
#
# nonlinear formulae with unknowns
# log link
fn2a <- finterp(~exp(b0+b1*x1+b2*x2))
fn2a
fn2a(rep(0.2,3))
# compartment model
times <- 1:20
# exp() parameters to ensure that they are positive
fn3 <- finterp(~exp(volume)*exp(absorption)/(exp(absorption)-
	exp(elimination))*(exp(-exp(elimination)*times)-
	exp(-exp(absorption)*times)))
fn3
fn3(log(c(3,0.3,0.2)))
#
# Poisson density
y <- rpois(20,5)
fn4 <- finterp(~mu^y*exp(-mu)/gamma(y+1))
fn4
fn4(5)
dpois(y,5)
#
# Poisson likelihood
# mean parameter
fn5 <- finterp(~-y*log(mu)+mu+lgamma(y+1),vector=F)
fn5
likefn1 <- function(p) sum(fn5(mu=p))
nlm(likefn1,p=1)
mean(y)
# canonical parameter
fn5a <- finterp(~-y*theta+exp(theta)+lgamma(y+1),vector=F)
fn5a
likefn1a <- function(p) sum(fn5a(theta=p))
nlm(likefn1a,p=1)
#
# likelihood for Poisson log linear regression
y <- rpois(20,fn2a(c(0.2,1,0.4)))
nlm(likefn1,p=1)
mean(y)
likefn2 <- function(p) sum(fn5(mu=fn2a(p)))
nlm(likefn2,p=c(1,0,0))
# or
likefn2a <- function(p) sum(fn5a(theta=fn2(p)))
nlm(likefn2a,p=c(1,0,0))
#
# likelihood for Poisson nonlinear regression
y <- rpois(20,fn3(log(c(3,0.3,0.2))))
nlm(likefn1,p=1)
mean(y)
likefn3 <- function(p) sum(fn5(mu=fn3(p)))
nlm(likefn3,p=log(c(1,0.4,0.1)))
#
# envir as data objects
y <- matrix(rnorm(20),ncol=5)
y[3,3] <- y[2,2] <- NA
x <- 1:4
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
# Wilkinson and Rogers notation
form1 <- ~x+z2
class(form1) <- c(class(form1),"repeated")
print(fn1 <- finterp(form1, envir=reps))
fn1(2:4)
# with unknown parameters
form2 <- ~a+b*x+d*z2
class(form2) <- c(class(form2),"repeated")
print(fn2 <- finterp(form2, envir=reps))
fn2(2:4)
#
# time-constant covariates
# Wilkinson and Rogers notation
form3 <- ~x
class(form3) <- c(class(form3),"tccov")
print(fn3 <- finterp(form3, envir=reps))
fn3(2:3)
print(fn3a <- finterp(form3, envir=xx))
fn3a(2:3)
# with unknown parameters
form4 <- ~a+b*x
class(form4) <- c(class(form4),"tccov")
print(fn4 <- finterp(form4, envir=reps))
fn4(2:3)
print(fn4a <- finterp(form4, envir=xx))
fn4a(2:3)
#
# time-varying covariates
# Wilkinson and Rogers notation
form5 <- ~z1+z2
class(form5) <- c(class(form5),"tvcov")
print(fn5 <- finterp(form5, envir=reps))
fn5(2:4)
print(fn5a <- finterp(form5, envir=zz))
fn5a(2:4)
# with unknown parameters
form6 <- ~a+b*z1+d*z2
class(form6) <- c(class(form6),"tvcov")
print(fn6 <- finterp(form6, envir=reps))
fn6(2:4)
print(fn6a <- finterp(form6, envir=zz))
fn6a(2:4)
# with a variable not in the data object
form7 <- ~a+b*z1+d*z2+e*z3
class(form7) <- c(class(form7),"tvcov")
print(fn7 <- finterp(form7, envir=reps))
fn7(2:5)
}
