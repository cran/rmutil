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
#     int(f, a="-infty", b="infty", type="Romberg", eps=1.0e-6,
#	max, d, p=0)
#
#  DESCRIPTION
#
#    A function to perform vectorized Romberg integration

.First.lib <- function(lib, pkg) {
	library.dynam("rmutil", pkg, lib)
	provide(rmutil)
}

int <- function(f, a="-infty", b="infty", type="Romberg", eps=1.0e-6,
	max, d, p=0)
{
int1 <- function(ff, aa, bb){
	z <- .C("romberg",
		ff,
		as.double(aa),
		as.double(bb),
		len=as.integer(len),
		eps=as.double(eps),
		pts=as.integer(d),
		max=as.integer(max),
		err=integer(1),
		res=double(len))
	if(z$err==1)warning("Unable to allocate memory for int")
	if(z$err==2)warning("Division by zero in int")
	else if(z$err==3)warning("No convergence in int")
	z$res}
if(!missing(type))type <- match.arg(type,c("Romberg","TOMS614"))
if(missing(max)) max <- if(type=="Romberg") 16 else 100
if(missing(d)) d <- if(type=="Romberg") 5 else 1
if(is.numeric(a)){
	if(is.numeric(b)){
		if(any(a>=b))stop("Some a>=b")
		len <- length(f((a+b)/2))}
	else len <- length(f(a+1))}
else if(is.numeric(b))len <- length(f(b-1))
else len <- length(f(0))
if(is.numeric(a)&length(a)!=len){
	if(length(a)!=1)stop("a has incorrect length")
	else a <- rep(a,len)}
if(is.numeric(b)&length(b)!=len){
	if(length(b)!=1)stop("b has incorrect length")
	else b <- rep(b,len)}
if(len>1&type!="Romberg")stop("vector functions only allowed with Romberg")
if(type=="Romberg"){
	ff <- function(x) f(1/x)/(x*x)
	if(!is.numeric(b)){
		if(!is.numeric(a)) z <- int1(ff,rep(-1,len),rep(0,len)) +
			int1(f,rep(-1,len),rep(1,len)) +
			int1(ff,rep(0,len),rep(1,len))
		else {
			if(any(a>0)){
				if(any(a<=0))a1 <- ifelse(a>0,a,1)
				else a1 <- a
				z1 <- int1(ff,rep(0,len), 1/a1)}
			else z1 <- rep(0,len)
			if(any(a<=0)){
				if(any(a>0))a1 <- ifelse(a<=0,a,0)
				else a1 <- a
				z2 <- int1(f,a1,rep(1,len)) +
				   int1(ff,rep(0,len),rep(1,len))}
			else z2 <- rep(0,len)
			z <- z1*(a>0)+z2*(a<=0)}}
	else if(!is.numeric(a)){
		if(any(b<0)){
			if(any(b>=0))b1 <- ifelse(b<0,b,1)
			else b1 <- b
			z1 <- int1(ff, 1/b1,rep(0,len))}
		else z1 <- rep(0,len)
		if(any(b>=0)){
			if(any(b<0))b1 <- ifelse(b>=0,b,0)
			else b1 <- b
			z2 <- int1(f,rep(-1,len), b1) +
				int1(ff,rep(-1,len),rep(0,len))}
		else z2 <- rep(0,len)
		z <- z1*(b<0)+z2*(b>=0)}
	else z <- int1(f, a, b)
	z}
else {
	left <- !is.numeric(a)&&is.numeric(b)
	if(!is.numeric(b)){
		if(!is.numeric(a)){
			inf <- 1
			a <- b <- 1}
		else {
			inf <- 2
			b <- 1}}
	else if(!is.numeric(a)){
		a <- 1
		inf <- 1}
	else inf <- 4
	if(left){
		z2 <- .Fortran("inthp",
		a=as.double(b),
		b=as.double(b),
		d=as.double(d),
		f=f,
		m=as.integer(max),
		p=as.double(p),
		eps=as.double(eps),
		inf=as.integer(2),
		quadr=as.double(1),
		DUP=F)
		if(z2$inf==3||z2$inf==4)warning(paste("error",z2$inf,"- integration incomplete - try larger max"))
		else if(z2$inf>4)stop(paste("error",z2$inf,"- incorrect arguments"))}
	z1 <- .Fortran("inthp",
		a=as.double(a),
		b=as.double(b),
		d=as.double(d),
		f=f,
		m=as.integer(max),
		p=as.double(p),
		eps=as.double(eps),
		inf=as.integer(inf),
		quadr=as.double(1),
		DUP=F)
	if(z1$inf==3||z1$inf==4)warning(paste("error",z1$inf,"- integration incomplete - try larger max"))
	else if(z1$inf>4)stop(paste("error",z1$inf,"- incorrect arguments"))
	if(left)z1$quadr <- z1$quadr-z2$quadr
	z1$quadr}}

	
	