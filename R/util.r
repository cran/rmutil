#
#  rmutil : A Library of Special Functions for Repeated Measurements
#  Copyright (C) 1998 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     wr(formula)
#     det(x)
#     collapse(x,index,fcn=sum)
#     mexp(x, t=1, n=20, k=3)
#
#  DESCRIPTION
#
#    Utility functions for repeated measurements

wr <- function(formula){
	mt <- terms(formula)
	mf <- model.frame(mt, sys.frame(sys.parent()), na.action=na.fail)
	list(response=model.response(mf, "numeric"),
		design=model.matrix(mt, mf))}

# det <- function(x) abs(Re(prod(eigen(x,only.values=T)$values)))
det <- function(x) abs(prod(diag(qr(x)$qr)))

collapse <- function(x,index,fcn=sum){
	ans <- NULL
	for(i in split(x,index))ans <- c(ans,fcn(i))
	ans}

mexp <- function(x, type="spectral decomposition", t=1, n=20, k=3){
	if(!is.matrix(x))stop("x must be a matrix")
	if(length(dim(x))!=2)stop("x must be a two dimensional matrix")
	if(dim(x)[1]!=dim(x)[2])stop("x must be a square matrix")
	type <- match.arg(type,c("spectral decomposition","series approximation"))
	d <- ncol(x)
	if(type=="spectral decomposition"){
		z <- eigen(t*x,sym=F)
		p <- z$vectors%*%diag(exp(z$values))%*%solve(z$vectors)}
	else {
		xx <- x*t/2^k
		p <- diag(d)
		q <- p
		for(r in 1:n){
			q <- xx%*%q/r
			p <- p+q}
		for(i in 1:k) p <- p%*%p}
	p}
