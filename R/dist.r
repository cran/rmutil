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
#     pinvgauss(q, m, s)
#     dinvgauss(y, m, s)
#     plaplace(q, m=0, s=1)
#     dlaplace(y, m=0, s=1)
#     plevy(y, m=0, s=1)
#     dlevy(y, m=0, s=1)
#     ppareto(q, m, s)
#     dpareto(y, m, s)
#
#     pboxcox(y, m, s, f)
#     dboxcox(y, m, s, f)
#     pburr(y, m, s, f)
#     dburr(y, m, s, f)
#     pgextval(y, s, m, f)
#     dgextval(y, s, m, f)
#     pggamma(y, s, m, f)
#     dggamma(y, s, m, f)
#     pginvgauss(y, m, s, f)
#     dginvgauss(y, m, s, f)
#     pglogis(y, m, s, f)
#     dglogis(y, m, s, f)
#     pgweibull(y, s, m, f)
#     dgweibull(y, s, m, f)
#     phjorth(y, m, s, f)
#     dhjorth(y, m, s, f)
#     ppowexp(y, m, s, f)
#     dpowexp(y, m, s, f)
#
#     pdoublepois(q, m, s)
#     ddoublepois(y, m, s)
#     pmultpois(q, m, s)
#     dmultpois(y, m, s)
#     pgammacount(q, m, s)
#     dgammacount(y, m, s)
#     pdoublebinom(q, n, m, s)
#     ddoublebinom(y, n, m, s)
#     pmultbinom(q, n, m, s)
#     dmultbinom(y, n, m, s)
#     pbetabinom(q, n, m, s)
#     dbetabinom(y, n, m, s)
#
#  DESCRIPTION
#
#    Functions to compute the probability and cumulative probability
# functions for
# continuous two parameter distributions:
#  inverse Gaussian, Laplace, Levy, Pareto
# continuous three parameter distributions:
#  Box-Cox, Burr, generalized extreme value, generalized gamma, generalized
#  inverse Gaussian, generalized logistic, generalized Weibull, Hjorth
# discrete two parameter distributions:
#  double Poisson, multiplicative Poisson, double binomial,
#  multiplicative binomial, beta binomial 

# continuous two parameter distributions

pinvgauss <- function(q, m, s){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	t <- q/m
	v <- sqrt(q*s)
	pnorm((t-1)/v)+exp(2/(m*s))*pnorm(-(t+1)/v)}

dinvgauss <- function(y, m, s){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	exp(-(y-m)^2/(2*y*s*m^2))/sqrt(2*pi*s*y^3)}

plaplace <- function(q, m=0, s=1){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(s<=0))stop("s must be positive")
	u <- (q-m)/s
	t <- exp(-abs(u))/2
	ifelse(u<0,t,1-t)}

dlaplace <- function(y, m=0, s=1){
	if(any(s<=0))stop("s must be positive")
	exp(-abs(y-m)/s)/(2*s)}

plevy <- function(q, m=0, s=1){
	if(any(q<m))stop("some y <= m")
	if(any(s<=0))stop("s must be positive")
	len <- length(q)
	if(length(m)!=len){
		if(length(m)!=1)stop("m has incorrect length")
		else m <- rep(m,len)}
	if(length(s)!=len){
		if(length(s)!=1)stop("s has incorrect length")
		else s <- rep(s,len)}
	z <- .C("plevy",
		as.double(q),
		as.double(m),
		as.double(s),
		as.double(1),
		len=as.integer(len),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(len),
		DUP=F)
	if(z$err==1)warning("Unable to allocate memory for integration")
	if(z$err==2)warning("Division by zero in integration")
	else if(z$err==3)warning("No convergence in integration")
	z$res}

dlevy <- function(y, m=0, s=1){
	if(any(y<=m))stop("some y <= m")
	if(any(s<=0))stop("s must be positive")
	sqrt(s/(2*pi*(y-m)^3))*exp(-s/(2*(y-m)))}

ppareto <- function(q, m, s){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=1))stop("s must be > 1")
	1-(1+q/(m*(s-1)))^-s}

dpareto <- function(y, m, s){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=1))stop("s must be > 1")
	m <- m*(s-1)
	s*(1+y/m)^(-s-1)/m}

# continuous three parameter distributions

# normed to make it a real distribution with y > 0
pboxcox <- function(q, m, s, f){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	norm <- sign(f)*pnorm(0,m,sqrt(s))
	ind <- f<0
	(pnorm(q^f/f,m,sqrt(s))-(f>0)*norm)/(1-ind-norm)}

dboxcox <- function(y, m, s, f){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	norm <- sign(f)*pnorm(0,m,sqrt(s))
	ind <- f<0
	y^(f-1)*dnorm(y^f/f,m,sqrt(s))/(1-ind-norm)}

pburr <- function(q, m, s, f){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	if(any(f<=0))stop("f must be positive")
	1-(1+(q/m)^s/f)^-f}

dburr <- function(y, m, s, f){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	if(any(f<=0))stop("f must be positive")
	y1 <- y/m
	s*y1^(s-1)/(m^s*(1+y1^s/f)^(f+1))}

# normed to make it a real distribution with y > 0
pgextval <- function(q, s, m, f){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	norm <- sign(f)*exp(-m^-s)
	ind <- f>0
	(pweibull(exp(q^f/f),s,m)-ind+(f>0)*norm)/(1-ind+norm)}

dgextval <- function(y, s, m, f){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	norm <- sign(f)*exp(-m^-s)
	ind <- f>0
	y1 <- exp(y^f/f)
	y^(f-1)*y1*dweibull(y1,s,m)/(1-ind+norm)}

pggamma <- function(q, s, m, f){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	if(any(f<=0))stop("f must be positive")
	pgamma(q^f,s,(m/s)^f)}

dggamma <- function(y, s, m, f){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	if(any(f<=0))stop("f must be positive")
	f*y^(f-1)*dgamma(y^f,s,(m/s)^f)}

pginvgauss <- function(q, m, s, f){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	len <- length(q)
	if(length(m)!=len){
		if(length(m)!=1)stop("m has incorrect length")
		else m <- rep(m,len)}
	if(length(s)!=len){
		if(length(s)!=1)stop("s has incorrect length")
		else s <- rep(s,len)}
	if(length(f)!=len){
		if(length(f)!=1)stop("f has incorrect length")
		else f <- rep(f,len)}
	z <- .C("pginvgauss",
		as.double(q),
		as.double(m),
		as.double(s),
		as.double(f),
		len=as.integer(len),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(len),
		DUP=F)
	if(z$err==1)warning("Unable to allocate memory for integration")
	if(z$err==2)warning("Division by zero in integration")
	else if(z$err==3)warning("No convergence in integration")
	z$res}

dginvgauss <- function(y, m, s, f){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	y^(f-1)*exp(-(1/y+y/m^2)/(2*s))/(m^f*(2*besselK(1/(s*m),abs(f))))}

pglogis <- function(q, m, s, f){
	if(any(s<=0))stop("s must be positive")
	if(any(f<=0))stop("f must be positive")
	(1+exp(-sqrt(3)*(q-m)/(s*pi)))^-f}

dglogis <- function(y, m, s, f) {
	if(any(s<=0))stop("s must be positive")
	if(any(f<=0))stop("f must be positive")
	y1 <- exp(-sqrt(3)*(y-m)/(s*pi))
	sqrt(3)*f*y1/(pi*s*(1+y1)^(f+1))}

pgweibull <- function(q, s, m, f){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	if(any(f<=0))stop("f must be positive")
	(1-exp(-(q/m)^s))^f}

dgweibull <- function(y, s, m, f){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	if(any(f<=0))stop("f must be positive")
	y1 <- exp(-(y/m)^s)
	s*f*y^(s-1)*(1-y1)^(f-1)*y1/m^s}

phjorth <- function(q, m, s, f){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	1-(1+s*q)^(-f/s)*exp(-(q/m)^2/2)}

dhjorth <- function(y, m, s, f){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	(1+s*y)^(-f/s)*exp(-(y/m)^2/2)*(y/m^2+f/(1+s*y))}

ppowexp <- function(q, m=0, s=1, f=1){
	if(any(s<=0))stop("s must be positive")
	if(any(f<=0))stop("f must be positive")
	len <- length(q)
	if(length(m)!=len){
		if(length(m)!=1)stop("m has incorrect length")
		else m <- rep(m,len)}
	if(length(s)!=len){
		if(length(s)!=1)stop("s has incorrect length")
		else s <- rep(s,len)}
	if(length(f)!=len){
		if(length(f)!=1)stop("f has incorrect length")
		else f <- rep(f,len)}
	z <- .C("ppowexp",
		as.double(q),
		as.double(m),
		as.double(s),
		as.double(f),
		len=as.integer(len),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(len),
		DUP=F)
	if(z$err==1)warning("Unable to allocate memory for integration")
	if(z$err==2)warning("Division by zero in integration")
	else if(z$err==3)warning("No convergence in integration")
	ifelse(q-m>0,0.5+z$res,0.5-z$res)}

dpowexp <- function(y, m=0, s=1, f=1){
	if(any(s<=0))stop("s must be positive")
	if(any(f<=0))stop("f must be positive")
	s <- sqrt(s)
	b <- 1+1/(2*f)
	exp(-(abs(y-m)/s)^(2*f)/2)/(s*gamma(b)*2^b)}

# discrete (overdispersed) two parameter distributions

pdoublepois <- function(q, m, s){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	if(length(m)!=length(q)){
		if(length(m)==1)m <- rep(m,length(q))
		else stop("m and q must have the same length")}
	if(length(s)!=length(q)){
		if(length(s)==1)s <- rep(s,length(q))
		else stop("s and q must have the same length")}
	.C("pdp",
		as.integer(q),
		as.integer(3*max(c(q,100))),
		as.double(m),
		as.double(s),
		as.integer(length(q)),
		res=double(length(q)),
		DUP=F)$res}

ddoublepois <- function(y, m, s){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	if(length(m)!=length(y)){
		if(length(m)==1)m <- rep(m,length(y))
		else stop("m and y must have the same length")}
	if(length(s)!=length(y)){
		if(length(s)==1)s <- rep(s,length(y))
		else stop("s and y must have the same length")}
	exp(.C("ddp",
		as.integer(y),
		as.integer(3*max(c(y,100))),
		as.double(m),
		as.double(s),
		as.integer(length(y)),
		as.double(rep(1,length(y))),
		res=double(length(y)),
		DUP=F)$res)}

pmultpois <- function(q, m, s){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	if(length(m)!=length(q)){
		if(length(m)==1)m <- rep(m,length(q))
		else stop("m and q must have the same length")}
	if(length(s)!=length(q)){
		if(length(s)==1)s <- rep(s,length(q))
		else stop("s and q must have the same length")}
	.C("pmp",
		as.integer(q),
		as.integer(3*max(c(q,100))),
		as.double(m),
		as.double(s),
		as.integer(length(q)),
		res=double(length(q)),
		DUP=F)$res}

dmultpois <- function(y, m, s){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	if(length(m)!=length(y)){
		if(length(m)==1)m <- rep(m,length(y))
		else stop("m and y must have the same length")}
	if(length(s)!=length(y)){
		if(length(s)==1)s <- rep(s,length(y))
		else stop("s and y must have the same length")}
	exp(.C("dmp",
		as.integer(y),
		as.integer(3*max(c(y,100))),
		as.double(m),
		as.double(s),
		as.integer(length(y)),
		as.double(rep(1,length(y))),
		res=double(length(y)),
		DUP=F)$res)}

pgammacount <- function(q, m, s){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	1-pgamma(m*s,(q+1)*s,1)}

dgammacount <- function(y, m, s){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(m<=0))stop("m must be positive")
	if(any(s<=0))stop("s must be positive")
	ifelse(y==0,1-pgamma(m*s,(y+1)*s,1),
		pgamma(m*s,y*s+(y==0),1)-pgamma(m*s,(y+1)*s,1))}

pdoublebinom <- function(q, n, m, s){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(n<0))stop("n must contain non-negative values")
	if(length(n)!=length(q)){
		if(length(n)==1)n <- rep(n,length(q))
		else stop("n must be the same length as q")}
	if(any(q>n))stop("q must be <= n")
	if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
	if(any(s<=0))stop("s must be positive")
	if(length(m)!=length(q)){
		if(length(m)==1)m <- rep(m,length(q))
		else stop("m and q must have the same length")}
	if(length(s)!=length(q)){
		if(length(s)==1)s <- rep(s,length(q))
		else stop("s and q must have the same length")}
	.C("pdb",
		as.integer(q),
		as.integer(n),
		as.double(m),
		as.double(s),
		as.integer(length(q)),
		res=double(length(q)),
		DUP=F)$res}

ddoublebinom <- function(y, n, m, s){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(n<0))stop("n must contain non-negative values")
	if(length(n)!=length(y)){
		if(length(n)==1)n <- rep(n,length(y))
		else stop("n must be the same length as y")}
	if(any(y>n))stop("y must be <= n")
	if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
	if(any(s<=0))stop("s must be positive")
	if(length(m)!=length(y)){
		if(length(m)==1)m <- rep(m,length(y))
		else stop("m and y must have the same length")}
	if(length(s)!=length(y)){
		if(length(s)==1)s <- rep(s,length(y))
		else stop("s and y must have the same length")}
	exp(.C("ddb",
		as.integer(y),
		as.integer(n),
		as.double(m),
		as.double(s),
		as.integer(length(y)),
		as.double(rep(1,length(y))),
		res=double(length(y)),
		DUP=F)$res)}

pmultbinom <- function(q, n, m, s){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(n<0))stop("n must contain non-negative values")
	if(length(n)!=length(q)){
		if(length(n)==1)n <- rep(n,length(q))
		else stop("n must be the same length as q")}
	if(any(q>n))stop("q must be <= n")
	if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
	if(any(s<=0))stop("s must be positive")
	if(length(m)!=length(q)){
		if(length(m)==1)m <- rep(m,length(q))
		else stop("m and q must have the same length")}
	if(length(s)!=length(q)){
		if(length(s)==1)s <- rep(s,length(q))
		else stop("s and q must have the same length")}
	.C("pmb",
		as.integer(q),
		as.integer(n),
		as.double(m),
		as.double(s),
		as.integer(length(q)),
		res=double(length(q)),
		DUP=F)$res}

dmultbinom <- function(y, n, m, s){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(n<0))stop("n must contain non-negative values")
	if(length(n)!=length(y)){
		if(length(n)==1)n <- rep(n,length(y))
		else stop("n must be the same length as y")}
	if(any(y>n))stop("y must be <= n")
	if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
	if(any(s<=0))stop("s must be positive")
	if(length(m)!=length(y)){
		if(length(m)==1)m <- rep(m,length(y))
		else stop("m and y must have the same length")}
	if(length(s)!=length(y)){
		if(length(s)==1)s <- rep(s,length(y))
		else stop("s and y must have the same length")}
	exp(.C("dmb",
		as.integer(y),
		as.integer(n),
		as.double(m),
		as.double(s),
		as.integer(length(y)),
		as.double(rep(1,length(y))),
		res=double(length(y)),
		DUP=F)$res)}

pbetabinom <- function(q, n, m, s){
	if(any(q<0))stop("q must contain non-negative values")
	if(any(n<0))stop("n must contain non-negative values")
	if(length(n)!=length(q)){
		if(length(n)==1)n <- rep(n,length(q))
		else stop("n must be the same length as q")}
	if(any(q>n))stop("q must be <= n")
	if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
	if(any(s<=0))stop("s must be positive")
	if(length(m)!=length(q)){
		if(length(m)==1)m <- rep(m,length(q))
		else stop("m and q must have the same length")}
	if(length(s)!=length(q)){
		if(length(s)==1)s <- rep(s,length(q))
		else stop("s and q must have the same length")}
	t <- s*m
	u <- s*(1-m)
	res <- NULL
	for(i in 1:length(q)){
		qq <- 0:q[i]
		res <- c(res,sum(exp(lbeta(qq+t[i],n[i]-qq+u[i])-
			lbeta(t[i],u[i])+lchoose(n[i],qq))))}
	res}

dbetabinom <- function(y, n, m, s){
	if(any(y<0))stop("y must contain non-negative values")
	if(any(n<0))stop("n must contain non-negative values")
	if(length(n)!=length(y)){
		if(length(n)==1)n <- rep(n,length(y))
		else stop("n must be the same length as y")}
	if(any(y>n))stop("y must be <= n")
	if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
	if(any(s<=0))stop("s must be positive")
	if(length(m)!=length(y)&&length(m)!=1)
		stop("m and y must have the same length")
	if(length(s)!=length(y)&&length(s)!=1)
		stop("s and y must have the same length")
	t <- s*m
	u <- s*(1-m)
	exp(lbeta(y+t,n-y+u)-lbeta(t,u)+lchoose(n,y))}
