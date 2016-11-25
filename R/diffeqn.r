#
#  rmutil : A Library of Special Functions for Repeated Measurements
#  Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public Licence as published by
#  the Free Software Foundation; either version 2 of the Licence, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public Licence for more details.
#
#  You should have received a copy of the GNU General Public Licence
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     runge.kutta(f, initial, x)
#     lin.diff.eqn(A, initial, t=1)
#
#  DESCRIPTION
#
#	Functions for solving differential equations

###    Runge-Kutta method for solving a differential equation
###
runge.kutta <- function(f, initial, x){
if(!is.function(f))stop("f must be a function")
if(!is.numeric(initial)||length(initial)!=1)stop("initial must be a scalar")
if(!is.vector(x,mode="numeric"))stop("x must be a numeric vector")
y <- initial
for(i in 1:(length(x)-1)){
	stepsize <- x[i+1]-x[i]
	f1 <- stepsize*f(y[i],x[i])
	f2 <- stepsize*f(y[i]+f1/2,x[i]+stepsize/2)
	f3 <- stepsize*f(y[i]+f2/2,x[i]+stepsize/2)
	f4 <- stepsize*f(y[i]+f3,x[i]+stepsize)
	y <- c(y,y[i]+(f1+2*f2+2*f3+f4)/6)}
y
}

### autonomous linear differential equations
###
lin.diff.eqn <- function(A, initial, t=1){
if(!is.matrix(A)||dim(A)[1]!=dim(A)[2])stop("A must be a square matrix")
if(!is.vector(initial,mode="numeric")||length(initial)!=dim(A)[1])
	stop("initial must be a numeric vector with length equal to the dimensions of A")
if(!is.vector(t,mode="numeric"))stop("t must be a numeric scalar or vector")
z <- NULL
for(i in 1:length(t))z <- cbind(z,mexp(A,t=t[i])%*%initial)
t(z)
}
