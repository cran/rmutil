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
#     wr(formula, data=NULL, expand=T))
#     det(x)
#     capply(x, index, fcn=sum)
#     mexp(x, t=1, n=20, k=3)
#     contr.mean(n, contrasts=TRUE)
#
#  DESCRIPTION
#
#    Utility functions for repeated measurements

### function to find the response vector and design matrix from a W&R formula
###
wr <- function(formula, data=NULL, expand=TRUE){
if(is.null(data))data <- parent.frame()
else if(!is.data.frame(data)&&!is.environment(data))
	data <- if(expand||inherits(data,"tccov"))as.data.frame(data)
		else as.data.frame(data$ccov)
mt <- terms(formula)
mf <- model.frame(mt,data=data,na.action=NULL)
list(response=model.response(mf,"numeric"), design=model.matrix(mt,mf))}

### a fast simplified version of tapply
###
capply <- function(x, index, fcn=sum){
ans <- NULL
for(i in split(x,index))ans <- c(ans,fcn(i))
ans}

### matrix exponentiation
###
mexp <- function(x, t=1, type="spectral decomposition", n=20, k=3){
if(!is.matrix(x))stop("x must be a matrix")
if(dim(x)[1]!=dim(x)[2])stop("x must be a square matrix")
type <- match.arg(type,c("spectral decomposition","series approximation"))
if(type=="spectral decomposition"){
	z <- eigen(t*x,symmetric=FALSE)
	p <- z$vectors%*%diag(exp(z$values))%*%solve(z$vectors)}
else {
	xx <- x*t/2^k
	p <- diag(dim(x)[2])
	q <- p
	for(r in 1:n){
		q <- xx%*%q/r
		p <- p+q}
	for(i in 1:k) p <- p%*%p}
p}

### matrix power
###
"%^%" <- function(x, p){
if(!is.matrix(x))stop("x must be a matrix")
if(dim(x)[1]!=dim(x)[2])stop("x must be a square matrix")
z <- eigen(x,symmetric=FALSE)
z$vectors%*%diag(z$values^p)%*%solve(z$vectors)}

#    A function to provide correct constraints about the mean
#  (correcting contr.sum)

contr.mean <- function(n, contrasts=TRUE){
if(length(n) <= 1){
	if(is.numeric(n)&&length(n)==1&&n>1)levels <- 1:n
	else stop("Not enough degrees of freedom to define contrasts")}
else levels <- n
lenglev <- length(levels)
if(contrasts){
	cont <- array(0,c(lenglev,lenglev-1),
		list(levels,levels[1:(lenglev-1)]))
	cont[col(cont)==row(cont)] <- 1
	cont[lenglev,] <- -1}
else {
	cont <- array(0,c(lenglev,lenglev),list(levels,levels))
	cont[col(cont) == row(cont)] <- 1}
cont}
