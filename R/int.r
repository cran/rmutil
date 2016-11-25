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
#     int(f, a=-Inf, b=Inf, type="Romberg", eps=0.0001,
#	max=NULL, d=NULL, p=0)
#
#  DESCRIPTION
#
#    A function to perform vectorized Romberg integration
# Now using LazyLoad: true in DESCRIPTION
# http://stackoverflow.com/a/4369551/2727349
#.First.lib <- function(lib, pkg)
#	library.dynam("rmutil", pkg, lib)
###
### vectorized one-dimensional integration
###
int <- function(f, a=-Inf, b=Inf, type="Romberg", eps=0.0001,
	max=NULL, d=NULL, p=0){
#
# function to call the C code
#
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
		res=double(len),
		PACKAGE="rmutil")
	if(z$err==1)warning("Unable to allocate memory for int")
	if(z$err==2)warning("Division by zero in int")
	else if(z$err==3)warning("No convergence in int")
	z$res}
#
# check algorithm to be used and initialize parameters
#
type <- match.arg(type,c("Romberg","TOMS614"))
if(is.null(max))max <- if(type=="Romberg") 16 else 100
if(is.null(d))d <- if(type=="Romberg") 5 else 1
#
# check function and integration limits
#
if(length(formals(f))!=1)stop("f must have one argument")
if(any(a==Inf))stop("lower bound cannot be Inf")
if(any(b==-Inf))stop("upper bound cannot be -Inf")
if((length(a)>1&&any(a==-Inf)&&all(a!=-Inf))||(length(b)>1&&any(b==Inf)&&all(b!=Inf)))stop("int cannot be vectorized with some limits infinite")
#
# determine length of vector to be integrated
#
if(all(a!=-Inf)){
	if(all(b!=Inf)){
		if(any(a>=b))stop("some a>=b")
		len <- length(f((a+b)/2))}
	else len <- length(f(a+1))}
else if(all(b!=Inf))len <- length(f(b-1))
else len <- length(f(0))
if(len>1&&type!="Romberg")stop("vector functions only allowed with Romberg")
#
# if a vector and there are infinite limits, check that all limits are infinite
#
if(all(a!=-Inf)&&length(a)!=len){
	if(length(a)!=1)stop("a has incorrect length")
	else a <- rep(a,len)}
if(all(b!=Inf)&&length(b)!=len){
	if(length(b)!=1)stop("b has incorrect length")
	else b <- rep(b,len)}
if(type=="Romberg"){
	# invert function for infinite limits
	ff <- function(x) f(1/x)/(x*x)
	if(b==Inf){
		if(all(a==-Inf))
		# both limits infinite
			z <- int1(ff,rep(-1,len),rep(0,len))+
				int1(f,rep(-1,len),rep(1,len))+
				int1(ff,rep(0,len),rep(1,len))
		else {
		# only upper limit infinite, cut in 2 pieces about 0
			if(any(a>0)){
				if(any(a<=0))a1 <- ifelse(a>0,a,1)
				else a1 <- a
				z1 <- int1(ff,rep(0,len), 1/a1)}
			else z1 <- rep(0,len)
			if(any(a<=0)){
				if(any(a>0))a1 <- ifelse(a<=0,a,0)
				else a1 <- a
				z2 <- int1(f,a1,rep(1,len))+
					int1(ff,rep(0,len),rep(1,len))}
			else z2 <- rep(0,len)
			z <- z1*(a>0)+z2*(a<=0)}}
	else if(all(a==-Inf)){
	# only lower limit infinite, cut in 2 pieces about 0
		if(any(b<0)){
			if(any(b>=0))b1 <- ifelse(b<0,b,1)
			else b1 <- b
			z1 <- int1(ff, 1/b1,rep(0,len))}
		else z1 <- rep(0,len)
		if(any(b>=0)){
			if(any(b<0))b1 <- ifelse(b>=0,b,0)
			else b1 <- b
			z2 <- int1(f,rep(-1,len), b1)+
				int1(ff,rep(-1,len),rep(0,len))}
		else z2 <- rep(0,len)
		z <- z1*(b<0)+z2*(b>=0)}
	else z <- int1(f, a, b)
	z}
else {
#
# TOMS614
#
	left <- a==-Inf&&b!=Inf
	if(b==Inf){
		if(all(a==-Inf)){
		# both limits infinite
			inf <- 1
			a <- b <- 1}
		else {
		# only upper limit is infinite
			inf <- 2
			b <- 1}}
	else if(all(a==-Inf)){
	# only lower limit is infinite
		a <- 1
		inf <- 1}
	else inf <- 4
	if(left){
	# lower limit infinite, upper limit numeric: calculate for
	# whole real line first
		z2 <- .C("inthp",
		a=as.double(b),
		b=as.double(b),
		d=as.double(d),
		f=f,
		m=as.integer(max),
		p=as.double(p),
		eps=as.double(eps),
		inf=as.integer(2),
		quadr=as.double(1),
		## DUP=FALSE,
		PACKAGE="rmutil")
		if(z2$inf==3||z2$inf==4)warning(paste("error",z2$inf,"- integration incomplete - try larger max"))
		else if(z2$inf>4)stop(paste("error",z2$inf,"- incorrect arguments"))}
	# integrate either for both limits finite or with upper limit infinite
	z1 <- .C("inthp",
		a=as.double(a),
		b=as.double(b),
		d=as.double(d),
		f=f,
		m=as.integer(max),
		p=as.double(p),
		eps=as.double(eps),
		inf=as.integer(inf),
		quadr=as.double(1),
		## DUP=FALSE,
		PACKAGE="rmutil")
	if(z1$inf==3||z1$inf==4)warning(paste("error",z1$inf,"- integration incomplete - try larger max"))
	else if(z1$inf>4)stop(paste("error",z1$inf,"- incorrect arguments"))
	# if lower limit infinite, upper limit numeric, subtract upper
	# part from that for whole real line
	if(left)z1$quadr <- z1$quadr-z2$quadr
	z1$quadr}}
###
### vectorized two-dimensional integration
###
int2 <- function(f, a=c(-Inf,-Inf), b=c(Inf,Inf), eps=1.0e-6, max=16, d=5){
#
# function adapted from Gentleman and Ihaka (2000)
#    Jr Comp Graph Stat 9, 491-508
#
g <- function(y){
	fx <- function(x) f(x,y)
	romberg(fx,a[,2],b[,2])}
#
# function to call the C code
#
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
		res=double(len),
		PACKAGE="rmutil")
	if(z$err==1)warning("Unable to allocate memory for int2")
	if(z$err==2)warning("Division by zero in int2")
	else if(z$err==3)warning("No convergence in int2")
	z$res}
#
# function for Romberg integration
#
romberg <- function(f, a=-Inf, b=Inf){
	# invert function for infinite limits
        ff <- function(x) f(1/x)/(x*x)
        if(b==Inf){
        	if(all(a==-Inf))
        	# both limits infinite
        		z <- int1(ff,rep(-1,len),rep(0,len))+
				int1(f,rep(-1,len),rep(1,len))+
				int1(ff,rep(0,len),rep(1,len))
        	else {
        	# only upper limit infinite, cut in 2 pieces about 0
        		if(any(a>0)){
        			if(any(a<=0))a1 <- ifelse(a>0,a,1)
        			else a1 <- a
        			z1 <- int1(ff,rep(0,len), 1/a1)}
        		else z1 <- rep(0,len)
        		if(any(a<=0)){
        			if(any(a>0))a1 <- ifelse(a<=0,a,0)
        			else a1 <- a
        			z2 <- int1(f,a1,rep(1,len))+
					int1(ff,rep(0,len),rep(1,len))}
        		else z2 <- rep(0,len)
        		z <- z1*(a>0)+z2*(a<=0)}}
        else if(all(a==-Inf)){
        # only lower limit infinite, cut in 2 pieces about 0
        	if(any(b<0)){
        		if(any(b>=0))b1 <- ifelse(b<0,b,1)
        		else b1 <- b
        		z1 <- int1(ff, 1/b1,rep(0,len))}
        	else z1 <- rep(0,len)
        	if(any(b>=0)){
        		if(any(b<0))b1 <- ifelse(b>=0,b,0)
        		else b1 <- b
        		z2 <- int1(f,rep(-1,len), b1)+
        			int1(ff,rep(-1,len),rep(0,len))}
        	else z2 <- rep(0,len)
        	z <- z1*(b<0)+z2*(b>=0)}
        else z <- int1(f, a, b)
        z}
#
# check dimensions of limits
#
if(is.vector(a,mode="numeric")&&length(a)==2)a <- matrix(a,ncol=2)
else if(is.matrix(a)){
	if(dim(a)[2]!=2)stop("a must be 2-column matrix")}
else stop("a must be a 2-element vector or a 2-column matrix")
if(is.vector(b,mode="numeric")&&length(b)==2)b <- matrix(b,ncol=2)
else if(is.matrix(b)){
	if(dim(b)[2]!=2)stop("b must be 2-column matrix")}
else stop("b must be a 2-element vector or a 2-column matrix")
if((dim(a)[1]>1&&((any(a[,1]==-Inf)&&all(a[,1]!=-Inf))||(any(a[,2]==-Inf)&&all(a[,1]!=-Inf))))||(dim(b)[1]>1&&((any(b[,1]==Inf)&&all(b[,1]!=Inf))||(any(b[,2]==Inf)&&all(b[,1]!=Inf)))))stop("int2 cannot have only some limits infinite")
if(length(formals(f))!=2)stop("f must have two arguments")
#
# determine length of vectors to be integrated
#
if(all(a!=-Inf)){
	if(all(b!=Inf)){
		if(any(a[,1]>=b[,1])||any(a[,2]>=b[,2]))stop("some a>=b")
		len <- length(f((a[,1]+b[,1])/2,(a[,2]+b[,2])/2))}
	else len <- length(f(a[,1]+1,a[,2]+1))}
else if(all(b!=Inf))len <- length(f(b[,1]-1,b[,2]-1))
else len <- length(f(0,0))
#
# if a matrix and there are infinite limits, check that all limits are infinite
#
if(any(a!=-Inf)&&dim(a)[1]!=len){
	if(dim(a)[1]!=1)stop("a has incorrect size")
	else a <- matrix(rep(a,len),ncol=2,byrow=TRUE)}
if(any(b!=Inf)&&dim(b)[1]!=len){
	if(dim(b)[1]!=1)stop("b has incorrect size")
	else b <- matrix(rep(b,len),ncol=2,byrow=TRUE)}
#
# integrate
#
romberg(g,a[,1],b[,1])}
