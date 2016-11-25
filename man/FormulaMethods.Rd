\name{FormulaMethods}
\title{Methods for formulafn Functions}
\usage{
\method{covariates}{formulafn}(z, ...)
\method{formula}{formulafn}(x, ...)
model(z, ...)
parameters(z, ...)
\method{print}{formulafn}(x, ...)
}
\alias{covariates.formulafn}
\alias{formula.formulafn}
\alias{model}
\alias{model.formulafn}
\alias{parameters}
\alias{parameters.formulafn}
\alias{print.formulafn}
\description{
Methods for accessing the contents of a function created from formula
produced by \code{\link[rmutil]{finterp}} or a function modified by
\code{\link[rmutil]{fnenvir}}.

\code{covariates}: extract the names of the covariates.

\code{formula}: extract the formula used to produce the function
(\code{\link[rmutil]{finterp}} only).

\code{model}: extract the model function or model matrix if W&R
notation was used.

\code{parameters}: extract the names of the parameters.
}
\arguments{
\item{x,z}{A function of class, \code{formulafn}.}
\item{...}{Arguments to other functions.}

}
\value{
These methods extract information about functions of class, \code{formulafn},
created by \code{\link[rmutil]{finterp}} or \code{\link[rmutil]{fnenvir}}.
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{finterp}}, \code{\link[rmutil]{fnenvir}}.
}
\examples{
x1 <- rpois(20,2)
x2 <- rnorm(20)
#
# Wilkinson and Rogers formula with three parameters
fn1 <- finterp(~x1+x2)
fn1
covariates(fn1)
formula(fn1)
model(fn1)
parameters(fn1)
#
# nonlinear formula with unknowns
fn2 <- finterp(~exp(b0+b1*x1+b2*x2))
fn2
covariates(fn2)
formula(fn2)
model(fn2)
parameters(fn2)
#
# function transformed by fnenvir
fn3 <- fnenvir(function(p) p[1]+p[2]*x1)
covariates(fn3)
formula(fn3)
model(fn3)
parameters(fn3)
}
\keyword{manip}
