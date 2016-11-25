\name{finterp}
\title{Formula Interpreter}
\alias{finterp}
\alias{finterp.data.frame}
\alias{finterp.default}
\alias{finterp.repeated}
\alias{finterp.tccov}
\alias{finterp.tvcov}
\usage{
finterp(.z, ...)
\method{finterp}{default}(.z, .envir=parent.frame(), .formula=FALSE, .vector=TRUE,
	.args=NULL, .start=1, .name=NULL, .expand=TRUE, .intercept=TRUE,
	.old=NULL, .response=FALSE, ...)
}
\description{
\code{finterp} translates a model formula into a function of the
unknown parameters or of a vector of them. Such language formulae can
either be in Wilkinson and Rogers notation or be expressions
containing both known (existing) covariates and unknown (not existing)
parameters. In the latter, factor variables cannot be used and
parameters must be scalars.

The covariates in the formula are sought in the environment or in the
data object provided. If the data object has class, \code{repeated} or
\code{response}, then the key words, \code{times} will use the
response times from the data object as a covariate, \code{individuals}
will use the index for individuals as a factor covariate, and
\code{nesting} the index for nesting as a factor covariate. The latter
two only work for W&R notation.

Note that, in parameter displays, formulae in Wilkinson and Rogers
notation use variable names whereas those with unknowns use the names
of these parameters, as given in the formulae, and that the meaning of
operators (*, /, :, etc.) is different in the two cases.
}
\arguments{
\item{.z}{A model formula beginning with ~, either in Wilkinson and
Rogers notation or containing unknown parameters. If it contains
unknown parameters, it can have several lines so that, for example,
local variables can be assigned temporary values. In this case,
enclose the formula in curly brackets.}
\item{.envir}{The environment in which the formula is to be
interpreted or a data object of class, \code{repeated}, \code{tccov},
or \code{tvcov}.}
\item{.formula}{If TRUE and the formula is in Wilkinson and Rogers
notation, just returns the formula.}
\item{.vector}{If FALSE and the formula contains unknown parameters,
the function returned has them as separate arguments. If TRUE, it has
one argument, the unknowns as a vector, unless certain parameter names
are specified in \code{.args}. Always TRUE if \code{.envir} is
a data object.}
\item{.args}{If \code{.vector} is TRUE, names of parameters that are to
be function arguments and not included in the vector.}
\item{.start}{The starting index value of the parameter vector in the
function returned when \code{.vector} is TRUE.}
\item{.name}{Character string giving the name of the data object
specified by \code{.envir}. Ignored unless the latter is such an
object and only necessary when \code{finterp} is called within other
functions.}
\item{.expand}{If TRUE, expand functions with only time-constant
covariates to return one value per observation instead of one value
per individual. Ignored unless \code{.envir} is an object of class,
\code{repeated}.}
\item{.intercept}{If W&R notation is supplied and \code{.intercept=F},
a model function without intercept is returned.}
\item{.old}{The name of an existing object of class \code{formulafn}
which has common parameters with the one being created, or a list of
such objects. Only used if \code{.vector}=TRUE. The value of
\code{.start} should ensure that there is no conflict in indexing the
vector.}
\item{.response}{If TRUE, any response variable can be used in the
function. If FALSE, checks are made that the response is not also used
as a covariate.}
\item{...}{Arguments passed to other functions.}
}
\value{
A function, of class \code{formulafn}, of the unknown parameters or of a
vector of them is returned. Its attributes give the formula supplied,
the model function produced, the covariate names, the parameter names,
and the range of values of the index of the parameter vector. If
\code{formula} is TRUE and a Wilkinson and Rogers formula was
supplied, it is simply returned instead of creating a function.
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{FormulaMethods}}, \code{\link[rmutil]{covariates}},
\code{\link[rmutil]{fnenvir}}, \code{\link[rmutil]{formula}},
\code{\link[rmutil]{model}}, \code{\link[rmutil]{parameters}}
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
# parameters common to two functions
fn2b <- finterp(~c0+c1*exp(b0+b1*x1+b2*x2), .old=fn2a, .start=4)
fn2b
# function returned also depends on values of another function
fn2c <- finterp(~fn2+c1*exp(b0+b1*x1+b2*x2), .old=fn2a,
	.start=4, .args="fn2")
fn2c
args(fn2c)
fn2c(rep(0.2,4),fn2(rep(2,3)))
#
# compartment model
times <- 1:20
# exp() parameters to ensure that they are positive
fn3 <- finterp(~exp(absorption-volume)/(exp(absorption)-
	exp(elimination))*(exp(-exp(elimination)*times)-
	exp(-exp(absorption)*times)))
fn3
fn3(log(c(0.3,3,0.2)))
# a more efficient way
# (note that parameters do not appear in the same order)
form <- ~{
	ka <- exp(absorption)
	ke <- exp(elimination)
	ka*exp(-volume)/(ka-ke)*(exp(-ke*times)-exp(-ka*times))}
fn3a <- finterp(form)
fn3a(log(c(0.3,0.2,3)))
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
fn5 <- finterp(~-y*log(mu)+mu+lgamma(y+1),.vector=FALSE)
fn5
likefn1 <- function(p) sum(fn5(mu=p))
nlm(likefn1,p=1)
mean(y)
# canonical parameter
fn5a <- finterp(~-y*theta+exp(theta)+lgamma(y+1),.vector=FALSE)
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
x1 <- 1:4
x2 <- c("a","b","c","d")
resp <- restovec(y)
xx <- tcctomat(x1)
xx2 <- tcctomat(data.frame(x1,x2))
z1 <- matrix(rnorm(20),ncol=5)
z2 <- matrix(rnorm(20),ncol=5)
z3 <- matrix(rnorm(20),ncol=5)
zz <- tvctomat(z1)
zz <- tvctomat(z2,old=zz)
reps <- rmna(resp, ccov=xx, tvcov=zz)
reps2 <- rmna(resp, ccov=xx2, tvcov=zz)
rm(y, x1, x2 , z1, z2)
#
# repeated objects
#
# time-constant covariates
# Wilkinson and Rogers notation
form1 <- ~x1
print(fn1 <- finterp(form1, .envir=reps))
fn1(2:3)
print(fn1a <- finterp(form1, .envir=xx))
fn1a(2:3)
form1b <- ~x1+x2
print(fn1b <- finterp(form1b, .envir=reps2))
fn1b(2:6)
print(fn1c <- finterp(form1b, .envir=xx2))
fn1c(2:6)
# with unknown parameters
form2 <- ~a+b*x1
print(fn2 <- finterp(form2, .envir=reps))
fn2(2:3)
print(fn2a <- finterp(form2, .envir=xx))
fn2a(2:3)
#
# time-varying covariates
# Wilkinson and Rogers notation
form3 <- ~z1+z2
print(fn3 <- finterp(form3, .envir=reps))
fn3(2:4)
print(fn3a <- finterp(form3, .envir=zz))
fn3a(2:4)
# with unknown parameters
form4 <- ~a+b*z1+c*z2
print(fn4 <- finterp(form4, .envir=reps))
fn4(2:4)
print(fn4a <- finterp(form4, .envir=zz))
fn4a(2:4)
#
# note: lengths of x1 and z2 differ
# Wilkinson and Rogers notation
form5 <- ~x1+z2
print(fn5 <- finterp(form5, .envir=reps))
fn5(2:4)
# with unknown parameters
form6 <- ~a+b*x1+c*z2
print(fn6 <- finterp(form6, .envir=reps))
fn6(2:4)
#
# with times
# Wilkinson and Rogers notation
form7 <- ~x1+z2+times
print(fn7 <- finterp(form7, .envir=reps))
fn7(2:5)
form7a <- ~x1+x2+z2+times
print(fn7a <- finterp(form7a, .envir=reps2))
fn7a(2:8)
# with unknown parameters
form8 <- ~a+b*x1+c*z2+e*times
print(fn8 <- finterp(form8, .envir=reps))
fn8(2:5)
#
# with a variable not in the data object
form9 <- ~a+b*z1+c*z2+e*z3
print(fn9 <- finterp(form9, .envir=reps))
fn9(2:5)
# z3 assumed to be an unknown parameter:
fn9(2:6)
#
# multiline formula
form10 <- ~{
	tmp <- exp(b)
	a+tmp*z1+c*z2+d*times}
print(fn10 <- finterp(form10, .envir=reps))
fn10(2:5)
}
\keyword{programming}
