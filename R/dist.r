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
#     pinvgauss(q, m, s)
#     dinvgauss(y, m, s, log=FALSE)
#     qinvgauss(p, m, s)
#     rinvgauss(n, m, s)
#     plaplace(q, m=0, s=1)
#     dlaplace(y, m=0, s=1)
#     qlaplace(p, m=0, s=1)
#     rlaplace(n, m=0, s=1)
#     plevy(q, m=0, s=1)
#     dlevy(y, m=0, s=1)
#     qlevy(p, m=0, s=1)
#     rlevy(n, m=0, s=1)
#     ppareto(q, m, s)
#     dpareto(y, m, s)
#     qpareto(p, m, s)
#     rpareto(n, m, s)
#     psimplex(q, m, s)
#     dsimplex(y, m, s)
#     qsimplex(p, m, s)
#     rsimplex(n, m, s)
#     ptwosidedpower(q, m, s=2)
#     dtwosidedpower(y, m, s=2)
#     qtwosidedpower(p, m, s=2)
#     rtwosidedpower(n, m, s=2)
#
#     pboxcox(q, m=0, s=1, f=1)
#     dboxcox(y, m=0, s=1, f=1)
#     qboxcox(p, m=0, s=1, f=1)
#     rboxcox(n, m=0, s=1, f=1)
#     pburr(q, m, s, f)
#     dburr(y, m, s, f)
#     qburr(p, m, s, f)
#     rburr(n, m, s, f)
#     pgextval(q, s, m, f)
#     dgextval(y, s, m, f)
#     qgextval(p, s, m, f)
#     rgextval(n, s, m, f)
#     pggamma(q, s, m, f)
#     dggamma(y, s, m, f)
#     qggamma(p, s, m, f)
#     rggamma(n, s, m, f)
#     pginvgauss(q, m, s, f)
#     dginvgauss(y, m, s, f)
#     qginvgauss(p, m, s, f)
#     rginvgauss(n, m, s, f)
#     pglogis(q, m=0, s=1, f=1)
#     dglogis(y, m=0, s=1, f=1)
#     qglogis(p, m=0, s=1, f=1)
#     rglogis(n, m=0, s=1, f=1)
#     pgweibull(q, s, m, f)
#     dgweibull(y, s, m, f)
#     qgweibull(p, s, m, f)
#     rgweibull(n, s, m, f)
#     phjorth(q, m, s, f)
#     dhjorth(y, m, s, f)
#     qhjorth(p, m, s, f)
#     rhjorth(n, m, s, f)
#     ppowexp(q, m=0, s=1, f=1)
#     dpowexp(y, m=0, s=1, f=1)
#     qpowexp(p, m=0, s=1, f=1)
#     rpowexp(n, m=0, s=1, f=1)
#     pskewlaplace(q, m=0, s=1, f=1)
#     dskewlaplace(y, m=0, s=1, f=1)
#     qskewlaplace(p, m=0, s=1, f=1)
#     rskewlaplace(n, m=0, s=1, f=1)
#
#     pbetabinom(q, size, m, s)
#     dbetabinom(y, size, m, s)
#     qbetabinom(p, size, m, s)
#     rbetabinom(n, size, m, s)
#     pdoublebinom(q, size, m, s)
#     ddoublebinom(y, size, m, s)
#     qdoublebinom(p, size, m, s)
#     rdoublebinom(n, size, m, s)
#     pmultbinom(q, size, m, s)
#     dmultbinom(y, size, m, s)
#     qmultbinom(p, size, m, s)
#     rmultbinom(n, size, m, s)
#     pdoublepois(q, m, s)
#     ddoublepois(y, m, s)
#     qdoublepois(p, m, s)
#     rdoublepois(n, m, s)
#     pmultpois(q, m, s)
#     dmultpois(y, m, s)
#     qmultpois(p, m, s)
#     rmultpois(n, m, s)
#     ppvfpois(q, m, s, f)
#     dpvfpois(y, m, s, f)
#     qpvfpois(p, m, s, f)
#     rpvfpois(n, m, s, f)
#     pgammacount(q, m, s)
#     dgammacount(y, m, s)
#     qgammacount(p, m, s)
#     rgammacount(n, m, s)
#     pconsul(q, m, s)
#     dconsul(y, m, s)
#     qconsul(p, m, s)
#     rconsul(n, m, s)
#
#  DESCRIPTION
#
#    Functions to compute the probability and cumulative probability
# functions for
# continuous two parameter distributions:
#  inverse Gaussian, Laplace, Levy, Pareto, simplex
# continuous three parameter distributions:
#  Box-Cox, Burr, generalized extreme value, generalized gamma, generalized
#  inverse Gaussian, generalized logistic, generalized Weibull, Hjorth,
#  power exponential, skew Laplace
# discrete two parameter distributions:
#  double Poisson, multiplicative Poisson, double binomial, Consul,
#  multiplicative binomial, beta binomial 

### continuous two-parameter distributions
###
### inverse Gaussian distribution
###
pinvgauss <- function(q, m, s){
if(any(q<=0))stop("q must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
t <- q/m
v <- sqrt(q*s)
pnorm((t-1)/v)+exp(2/(m*s))*pnorm(-(t+1)/v)}

dinvgauss <- function(y, m, s, log=FALSE){
if(any(y<=0))stop("y must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
tmp <- -(y-m)^2/(2*y*s*m^2)-(log(2*pi*s)+3*log(y))/2
if(!log)tmp <- exp(tmp)
tmp}

qinvgauss <- function(p, m, s){
h <- function(y){
	t <- y/m[i]
	v <- sqrt(y*s[i])
	pnorm((t-1)/v)+exp(2/(m[i]*s[i]))*pnorm(-(t+1)/v)-p[i]}
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0))stop("m must be positive")
if(any(s<0))stop("s must be positive")
len <- max(length(p),length(m),length(s))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(.Machine$double.xmin,20)
	while(h(interval[1])*h(interval[2])>0)interval <- 2*interval
	tmp[i] <- uniroot(h,interval)$root}
tmp}

rinvgauss <- function(n=1, m, s) qinvgauss(runif(n),m=m,s=s)

### Laplace distribution
###
plaplace <- function(q, m=0, s=1){
if(any(s<=0))stop("s must be positive")
u <- (q-m)/s
t <- exp(-abs(u))/2
ifelse(u<0,t,1-t)}

dlaplace <- function(y, m=0, s=1, log=FALSE){
if(any(s<=0))stop("s must be positive")
tmp <- -abs(y-m)/s-log(2*s)
if(!log)tmp <- exp(tmp)
tmp}

qlaplace <- function(p, m=0, s=1){
h <- function(y){
  #u <- (y-m[i])/s[i]
  #t <- exp(-abs(u))/2
  #ifelse(u<0,t,1-t)-p[i]
## bruce edit:
  u <- (y-m)/s
  t <- exp(-abs(u))/2
  ifelse(u<0,t,1-t)-p
	}
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(s<=0))stop("s must be positive")
ifelse(p<0.5,s*log(2*p)+m,-s*log(2*(1-p))+m)}

rlaplace <- function(n = 1, m = 0, s = 1){
if(any(s<=0))stop("s must be positive")
q <- runif(n)
ifelse(q<0.5,s*log(2*q)+m,-s*log(2*(1-q))+m)}

### Levy distribution
###
plevy <- function(q, m=0, s=1){
if(any(q<=m))stop("some y <= m")
if(any(s<=0))stop("s must be positive")
2*(1-pnorm(1/sqrt((q-m)/s)))}

dlevy <- function(y, m=0, s=1, log=FALSE){
if(any(y<=m))stop("some y <= m")
if(any(s<=0))stop("s must be positive")
#sqrt(s/(2*pi*(y-m)^3))*exp(-s/(2*(y-m)))}
tmp <- log(s/(2*pi))/2-3*log(y-m)/2-s/(2*(y-m))
if(!log)tmp <- exp(tmp)
tmp}

qlevy <- function(p, m=0, s=1){
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(s<0))stop("s must be positive")
s/qnorm(1-p/2)^2+m}

rlevy <- function(n=1, m=0, s=1) {
if(any(s<0))stop("s must be positive")
s/qnorm(1-runif(n)/2)^2+m}

### Pareto distribution
###
ppareto <- function(q, m, s){
if(any(q<0))stop("q must be >= 0")
if(any(m<=0))stop("m must be positive")
if(any(s<=1))stop("s must be > 1")
1-(1+q/(m*(s-1)))^-s}

dpareto <- function(y, m, s, log=FALSE){
if(any(y<0))stop("y must be >= 0")
if(any(m<=0))stop("m must be positive")
if(any(s<=1))stop("s must be > 1")
m <- m*(s-1)
tmp <- log(s)-(s+1)*log(1+y/m)-log(m)
if(!log)tmp <- exp(tmp)
tmp}

qpareto <- function(p, m, s){
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0))stop("m must be positive")
if(any(s<=1))stop("s must be >1")
((1-p)^(-1/s)-1)*m*(s-1)}

rpareto <- function(n=1, m, s){
if(any(m<=0))stop("m must be positive")
if(any(s<=1))stop("s must be >1")
((1-runif(n))^(-1/s)-1)*m*(s-1)}

### simplex distribution
###
psimplex <- function(q, m, s){
if(any(q<=0)||any(q>=1))stop("q must contain values between 0 and 1")
if(any(m<=0)||any(m>=1))stop("m must contain values between 0 and 1")
if(any(s<=0))stop("s must be positive")
len <- max(length(q),length(m),length(s))
if(length(q)!=len){
	if(length(q)==1)q <- rep(q,len)
	else stop("length of q incorrect")}
if(length(m)!=len){
	if(length(m)!=1)stop("m has incorrect length")
	else m <- rep(m,len)}
if(length(s)!=len){
	if(length(s)!=1)stop("s has incorrect length")
	else s <- rep(s,len)}
z <- .C("psimplex_c",
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
	## DUP=FALSE,
	PACKAGE="rmutil")
if(z$err==1)warning("Unable to allocate memory for integration")
if(z$err==2)warning("Division by zero in integration")
else if(z$err==3)warning("No convergence in integration")
z$res}

dsimplex <- function(y, m, s, log=FALSE){
if(any(y<=0)||any(y>=1))stop("y must contain values between 0 and 1")
if(any(m<=0)||any(m>=1))stop("m must contain values between 0 and 1")
if(any(s<=0))stop("s must be positive")
tmp <- -((y-m)/(m*(1-m)))^2/(2*y*(1-y)*s)-(log(2*pi*s)+3*(log(y)+log(1-y)))/2
if(!log)tmp <- exp(tmp)
tmp}

qsimplex <- function(p, m, s){
h <- function(y).C("psimplex_c",
	as.double(y),
	as.double(m[i]),
	as.double(s[i]),
	as.double(1),
	len=as.integer(1),
	eps=as.double(1.0e-6),
	pts=as.integer(5),
	max=as.integer(16),
	err=integer(1),
	res=double(1),
	## DUP=FALSE,
	PACKAGE="rmutil")$res-p[i]
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0)||any(m>=1))stop("m must contain values between 0 and 1")
if(any(s<0))stop("s must be positive")
len <- max(length(p),length(m),length(s))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(0.000001,0.999999)
	while(h(interval[1])*h(interval[2])>0)
		interval <- c(interval[1]/10,1-interval[1]/10)
	tmp[i] <- uniroot(h,interval)$root}
tmp}

rsimplex <- function(n=1, m, s) qsimplex(runif(n),m=m,s=s)

### two-sided power distribution
###
ptwosidedpower <- function(q, m, s=2){
if(any(q<=0)||any(q>=1))stop("q must contain values between 0 and 1")
if(any(m<=0)||any(m>=1))stop("m must contain values between 0 and 1")
if(any(s<=0))stop("s must be positive")
ifelse(q<m,m*(q/m)^s,1-(1-m)*((1-q)/(1-m))^s)}

dtwosidedpower <- function(y, m, s=2, log=FALSE){
if(any(y<=0)||any(y>=1))stop("y must contain values between 0 and 1")
if(any(m<=0)||any(m>=1))stop("m must contain values between 0 and 1")
if(any(s<=0))stop("s must be positive")
tmp <- log(s)+(s-1)*ifelse(y<m,log(y)-log(m),log(1-y)-log(1-m))
if(!log)tmp <- exp(tmp)
tmp}

qtwosidedpower <- function(p, m, s=2){
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0)||any(m>=1))stop("m must contain values between 0 and 1")
if(any(s<0))stop("s must be positive")
ifelse(p<m,m*(p/m)^(1/s),1-(1-m)*((1-p)/(1-m))^(1/s))}

rtwosidedpower <- function(n, m, s=2) qtwosidedpower(runif(n),m=m,s=s)

### continuous three-parameter distributions
###
### Box-Cox distribution
### normed to make it a real distribution with y > 0
###
pboxcox <- function(q, m, s=1, f=1){
if(any(q<=0))stop("q must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
norm <- sign(f)*pnorm(0,m,sqrt(s))
(pnorm(q^f/f,m,sqrt(s))-(f>0)*norm)/(1-(f<0)-norm)}

dboxcox <- function(y, m, s=1, f=1, log=FALSE){
if(any(y<=0))stop("y must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
norm <- sign(f)*pnorm(0,m,sqrt(s))
tmp <- (f-1)*log(y)+dnorm(y^f/f,m,sqrt(s),log=TRUE)-log(1-(f<0)-norm)
if(!log)tmp <- exp(tmp)
tmp}

qboxcox <- function(p, m, s=1, f=1){
h <- function(y){
	norm <- sign(f[i])*pnorm(0,m[i],sqrt(s[i]))
	(pnorm(y^f[i]/f[i],m[i],sqrt(s[i]))-(f[i]>0)*norm)/
		(1-(f[i]<0)-norm)-p[i]}
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0))stop("m must be positive")
if(any(s<0))stop("s must be positive")
len <- max(length(p),length(m),length(s),length(f))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
if(length(f)!=len){
	if(length(f)==1)f <- rep(f,len)
	else stop("length of f incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(.Machine$double.xmin,20)
	while(h(interval[1])*h(interval[2])>0)interval <- 2*interval
	tmp[i] <- uniroot(h,interval)$root}
tmp}

rboxcox <- function(n=1, m, s=1, f=1) qboxcox(runif(n),m=m,s=s,f=f)

### Burr distribution
###
pburr <- function(q, m, s, f){
if(any(q<=0))stop("q must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
1-(1+(q/m)^s)^-f}

dburr <- function(y, m, s, f, log=FALSE){
if(any(y<=0))stop("y must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
y1 <- y/m
tmp <- log(f*s)+(s-1)*log(y1)-log(m)-(f+1)*log(1+y1^s)
if(!log)tmp <- exp(tmp)
tmp}

qburr <- function(p, m, s, f){
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
((1-p)^(-1/f)-1)^(1/s)*m}

rburr <- function(n=1, m, s, f){
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
((1-runif(n))^(-1/f)-1)^(1/s)*m}

### generalized extreme value distribution
### normed to make it a real distribution with y > 0
###
pgextval <- function(q, s, m, f){
if(any(q<=0))stop("q must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
norm <- sign(f)*exp(-m^-s)
ind <- f>0
(pweibull(exp(q^f/f),s,m)-ind+ind*norm)/(1-ind+norm)}

dgextval <- function(y, s, m, f, log=FALSE){
if(any(y<=0))stop("y must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
norm <- sign(f)*exp(-m^-s)
y1 <- exp(y^f/f)
tmp <- (f-1)*log(y)+log(y1)+dweibull(y1,s,m,log=TRUE)-log(1-(f>0)+norm)
if(!log)tmp <- exp(tmp)
tmp}

qgextval <- function(p, s, m, f){
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
norm <- sign(f)*exp(-m^-s)
ind <- f>0
(f*log(qweibull(p*(1-ind+norm)+ind-ind*norm,s,m)))^(1/f)}

rgextval <- function(n=1, s, m, f){
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
norm <- sign(f)*exp(-m^-s)
ind <- f>0
(f*log(qweibull(runif(n)*(1-ind+norm)+ind-ind*norm,s,m)))^(1/f)}

### generalized gamma distribution
###
pggamma <- function(q, s, m, f){
if(any(q<=0))stop("q must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
pgamma(q^f,s,scale=(m/s)^f)}

dggamma <- function(y, s, m, f, log=FALSE){
if(any(y<=0))stop("y must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
tmp <- log(f)+(f-1)*log(y)+dgamma(y^f,s,scale=(m/s)^f,log=TRUE)
if(!log)tmp <- exp(tmp)
tmp}

qggamma <- function(p, s, m, f) {
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
qgamma(p,s,scale=(m/s)^f)^(1/f)}

rggamma <- function(n=1, s, m, f){
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
qgamma(runif(n),s,scale=(m/s)^f)^(1/f)}

### generalized inverse Gaussian distribution
###
pginvgauss <- function(q, m, s, f){
if(any(q<=0))stop("q must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
len <- max(length(q),length(m),length(s))
if(length(q)!=len){
	if(length(q)==1)q <- rep(q,len)
	else stop("length of q incorrect")}
if(length(m)!=len){
	if(length(m)!=1)stop("m has incorrect length")
	else m <- rep(m,len)}
if(length(s)!=len){
	if(length(s)!=1)stop("s has incorrect length")
	else s <- rep(s,len)}
if(length(f)!=len){
	if(length(f)!=1)stop("f has incorrect length")
	else f <- rep(f,len)}
z <- .C("pginvgauss_c",
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
	## DUP=FALSE,
	PACKAGE="rmutil")
if(z$err==1)warning("Unable to allocate memory for integration")
if(z$err==2)warning("Division by zero in integration")
else if(z$err==3)warning("No convergence in integration")
z$res}

dginvgauss <- function(y, m, s, f, log=FALSE){
if(any(y<=0))stop("y must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
tmp <- (f-1)*log(y)-(1/y+y/m^2)/(2*s)-f*log(m)-log(2*besselK(1/(s*m),abs(f)))
if(!log)tmp <- exp(tmp)
tmp}

qginvgauss <- function(p, m, s, f){
h <- function(y).C("pginvgauss_c",
	as.double(y),
	as.double(m[i]),
	as.double(s[i]),
	as.double(f[i]),
	len=as.integer(1),
	eps=as.double(1.0e-6),
	pts=as.integer(5),
	max=as.integer(16),
	err=integer(1),
	res=double(1),
	## DUP=FALSE,
	PACKAGE="rmutil")$res-p[i]
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
len <- max(length(p),length(m),length(s),length(f))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
if(length(f)!=len){
	if(length(f)==1)f <- rep(f,len)
	else stop("length of f incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(.Machine$double.xmin,20)
	while(h(interval[1])*h(interval[2])>0)interval <- 2*interval
	tmp[i] <- uniroot(h,interval)$root}
tmp}

rginvgauss <- function(n=1, m, s, f) qginvgauss(runif(n),m=m,s=s,f=f)

### generalized logistic distribution
###
pglogis <- function(q, m=0, s=1, f=1){
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
(1+exp(-sqrt(3)*(q-m)/(s*pi)))^-f}

dglogis <- function(y, m=0, s=1, f=1, log=FALSE) {
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
y1 <- exp(-sqrt(3)*(y-m)/(s*pi))
tmp <- log(3)/2+log(f*y1)-log(pi*s)-(f+1)*log(1+y1)
if(!log)tmp <- exp(tmp)
tmp}

qglogis <- function(p, m=0, s=1, f=1){
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
-log(p^(-1/f)-1)*s*pi/sqrt(3)+m}

rglogis <- function(n=1, m=0, s=1, f=1){
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
-log(runif(n)^(-1/f)-1)*s*pi/sqrt(3)+m}

### generalized Weibull distribution
###
pgweibull <- function(q, s, m, f){
if(any(q<=0))stop("q must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
(1-exp(-(q/m)^s))^f}

dgweibull <- function(y, s, m, f, log=FALSE){
if(any(y<=0))stop("y must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
y1 <- exp(-(y/m)^s)
tmp <- log(s*f)+(s-1)*log(y)+(f-1)*log(1-y1)+log(y1)-s*log(m)
if(!log)tmp <- exp(tmp)
tmp}

qgweibull <- function(p, s, m, f){
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
m*(-log(1-p^(1/f)))^(1/s)}

rgweibull <- function(n=1, s, m, f){
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
m*(-log(1-runif(n)^(1/f)))^(1/s)}

### Hjorth distribution
###
phjorth <- function(q, m, s, f){
if(any(q<=0))stop("q must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
1-(1+s*q)^(-f/s)*exp(-(q/m)^2/2)}

dhjorth <- function(y, m, s, f, log=FALSE){
if(any(y<=0))stop("y must contain positive values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
tmp <- -(f/s)*log(1+s*y)-(y/m)^2/2+log(y/m^2+f/(1+s*y))
if(!log)tmp <- exp(tmp)
tmp}

qhjorth <- function(p, m, s, f){
h <- function(y) 1-(1+s[i]*y)^(-f[i]/s[i])*exp(-(y/m[i])^2/2)-p[i]
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0))stop("m must be positive")
if(any(s<0))stop("s must be positive")
len <- max(length(p),length(m),length(s),length(f))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
if(length(f)!=len){
	if(length(f)==1)f <- rep(f,len)
	else stop("length of f incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(.Machine$double.xmin,20)
	while(h(interval[1])*h(interval[2])>0)interval <- 2*interval
	tmp[i] <- uniroot(h,interval)$root}
tmp}

rhjorth <- function(n=1, m, s, f) qhjorth(runif(n),m=m,s=s,f=f)

### power exponential distribution
###
ppowexp <- function(q, m=0, s=1, f=1){
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
len <- max(length(q),length(m),length(s))
if(length(q)!=len){
	if(length(q)==1)q <- rep(q,len)
	else stop("length of q incorrect")}
if(length(m)!=len){
	if(length(m)!=1)stop("m has incorrect length")
	else m <- rep(m,len)}
if(length(s)!=len){
	if(length(s)!=1)stop("s has incorrect length")
	else s <- rep(s,len)}
if(length(f)!=len){
	if(length(f)!=1)stop("f has incorrect length")
	else f <- rep(f,len)}
z <- .C("ppowexp_c",
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
	## DUP=FALSE,
	PACKAGE="rmutil")
if(z$err==1)warning("Unable to allocate memory for integration")
if(z$err==2)warning("Division by zero in integration")
else if(z$err==3)warning("No convergence in integration")
ifelse(q-m>0,0.5+z$res,0.5-z$res)}

dpowexp <- function(y, m=0, s=1, f=1, log=FALSE){
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
s <- sqrt(s)
b <- 1+1/(2*f)
tmp <- -(abs(y-m)/s)^(2*f)/2-log(s)-lgamma(b)-b*log(2)
if(!log)tmp <- exp(tmp)
tmp}

qpowexp <- function(p, m=0, s=1, f=1){
h <- function(y) {
	z <- .C("ppowexp_c",
	        as.double(y),
	        as.double(m[i]),
	        as.double(s[i]),
	        as.double(f[i]),
	        len=as.integer(1),
	        eps=as.double(1.0e-6),
	        pts=as.integer(5),
	        max=as.integer(16),
	        err=integer(1),
	        res=double(1),
	        ## DUP=FALSE,
		PACKAGE="rmutil")$res
	if(y-m[i]>0) 0.5+z-p[i] else 0.5-z-p[i]}
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(s<0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
len <- max(length(p),length(m),length(s),length(f))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
if(length(f)!=len){
	if(length(f)==1)f <- rep(f,len)
	else stop("length of f incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- m[i]+s[i]*c(-2,2)
	while(h(interval[1])*h(interval[2])>0)interval <- 2*interval
	tmp[i] <- uniroot(h,interval)$root}
tmp}

rpowexp <- function(n=1, m=0, s=1, f=1) qpowexp(runif(n),m=m,s=s,f=f)

### skew Laplace distribution
###
pskewlaplace <- function(q, m=0, s=1, f=1){
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
u <- (q-m)/s
ifelse(u>0,1-exp(-f*abs(u))/(1+f^2),f^2*exp(-abs(u)/f)/(1+f^2))}

dskewlaplace <- function(y, m=0, s=1, f=1, log=FALSE){
if(any(s<=0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
tmp <- log(f)+ifelse(y>m,-f*(y-m),(y-m)/f)/s-log((1+f^2)*s)
if(!log)tmp <- exp(tmp)
tmp}

qskewlaplace <- function(p, m=0, s=1, f=1){
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(s<0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
ifelse(p<0.5,f*s*log((1+f^2)*p/f^2)+m,-s*log((1+f^2)*(1-p))/f+m)}

rskewlaplace <- function(n=1, m=0, s=1, f=1){
if(any(s<0))stop("s must be positive")
if(any(f<=0))stop("f must be positive")
q <- runif(n)
ifelse(q<0.5,f*s*log((1+f^2)*q/f^2)+m,-s*log((1+f^2)*(1-q))/f+m)}

### discrete (overdispersed) two-parameter distributions
###
### beta-binomial distribution
###
pbetabinom <- function(q, size, m, s){
if(any(q<0))stop("q must contain non-negative values")
if(any(size<0))stop("size must contain non-negative values")
if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
if(any(s<=0))stop("s must be positive")
len <- max(length(q),length(m),length(s),length(size))
if(length(q)!=len){
	if(length(q)==1)q <- rep(q,len)
	else stop("length of q incorrect")}
if(length(size)!=len){
	if(length(size)==1)size <- rep(size,len)
	else stop("size must be the same length as q")}
if(any(q>size))stop("q must be <= size")
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("m and q must have the same length")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("s and q must have the same length")}
t <- s*m
u <- s*(1-m)
res <- vector("numeric",length(q))
for(i in 1:length(q)){
	qq <- 0:q[i]
	res[i] <- sum(exp(lbeta(qq+t[i],size[i]-qq+u[i])-
		lbeta(t[i],u[i])+lchoose(size[i],qq)))}
res}

dbetabinom <- function(y, size, m, s, log=FALSE){
if(any(y<0))stop("y must contain non-negative values")
if(any(size<0))stop("size must contain non-negative values")
if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
if(any(s<=0))stop("s must be positive")
ly <- max(length(y),length(m),length(s),length(size))
if(length(y)!=ly){
	if(length(y)==1)y <- rep(y,ly)
	else stop("length of y incorrect")}
if(length(size)!=ly){
	if(length(size)==1)size <- rep(size,ly)
	else stop("size must be the same length as y")}
if(any(y>size))stop("y must be <= size")
if(length(m)!=ly){
	if(length(m)==1)m <- rep(m,ly)
	else stop("m and y must have the same length")}
if(length(s)!=ly){
	if(length(s)==1)s <- rep(s,ly)
	else stop("s and y must have the same length")}
t <- s*m
u <- s*(1-m)
tmp <- lbeta(y+t,size-y+u)-lbeta(t,u)+lchoose(size,y)
if(!log)tmp <- exp(tmp)
tmp}

qbetabinom <- function(p, size, m, s){
h <- function(y){
	t <- s[i]*m[i]
	u <- s[i]*(1-m[i])
	pp <- 0:y
	sum(exp(lbeta(pp+t,size[i]-pp+u)-
		lbeta(t,u)+lchoose(size[i],pp)))-p[i]}
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
if(any(s<0))stop("s must be positive")
len <- max(length(p),length(m),length(s),length(size))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(size)!=len){
	if(length(size)==1)size <- rep(size,len)
	else stop("length of size incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(0,size[i])
	tmp[i] <- if(h(interval[1])*h(interval[2])>0)0
		else uniroot(h,interval)$root}
round(tmp)}

rbetabinom <- function(n=1, size, m, s) qbetabinom(runif(n),size=size,m=m,s=s)

### double binomial distribution
###
pdoublebinom <- function(q, size, m, s){
if(any(q<0))stop("q must contain non-negative values")
if(any(size<0))stop("n must contain non-negative values")
if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
if(any(s<=0))stop("s must be positive")
len <- max(length(q),length(m),length(s),length(size))
if(length(q)!=len){
	if(length(q)==1)q <- rep(q,len)
	else stop("length of q incorrect")}
if(length(size)!=len){
	if(length(size)==1)size <- rep(size,len)
	else stop("size must be the same length as q")}
if(any(q>size))stop("q must be <= size")
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("m and q must have the same length")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("s and q must have the same length")}
.C("pdb",
	as.integer(q),
	as.integer(size),
	as.double(m),
	as.double(s),
	as.integer(length(q)),
	res=double(length(q)),
	## DUP=FALSE,
	PACKAGE="rmutil")$res}

ddoublebinom <- function(y, size, m, s, log=FALSE){
if(any(y<0))stop("y must contain non-negative values")
if(any(size<0))stop("size must contain non-negative values")
if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
if(any(s<=0))stop("s must be positive")
ly <- max(length(y),length(m),length(s),length(size))
if(length(y)!=ly){
	if(length(y)==1)y <- rep(y,ly)
	else stop("length of y incorrect")}
if(length(size)!=ly){
	if(length(size)==1)size <- rep(size,ly)
	else stop("size must be the same length as y")}
if(any(y>size))stop("y must be <= size")
if(length(m)!=ly){
	if(length(m)==1)m <- rep(m,ly)
	else stop("m and y must have the same length")}
if(length(s)!=ly){
	if(length(s)==1)s <- rep(s,ly)
	else stop("s and y must have the same length")}
tmp <- .C("ddb",
	as.integer(y),
	as.integer(size),
	as.double(m),
	as.double(s),
	as.integer(length(y)),
	as.double(rep(1,length(y))),
	res=double(length(y)),
	## DUP=FALSE,
	PACKAGE="rmutil")$res
if(!log)tmp <- exp(tmp)
tmp}

qdoublebinom <- function(p, size, m, s){
h <- function(y) .C("pdb",
	as.integer(y),
	as.integer(size[i]),
	as.double(m[i]),
	as.double(s[i]),
	as.integer(1),
	res=double(1),
	## DUP=FALSE,
	PACKAGE="rmutil")$res-p[i]
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
if(any(s<=0))stop("s must be positive")
len <- max(length(p),length(m),length(s),length(size))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(size)!=len){
	if(length(size)==1)size <- rep(size,len)
	else stop("length of size incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(0,size[i])
	tmp[i] <- if(h(interval[1])*h(interval[2])>0)0
		else uniroot(h,interval)$root}
round(tmp)}

rdoublebinom <- function(n=1, size, m, s)
	qdoublebinom(runif(n),size=size,m=m,s=s)

### multiplicative binomial distribution
###
pmultbinom <- function(q, size, m, s){
if(any(q<0))stop("q must contain non-negative values")
if(any(size<0))stop("size must contain non-negative values")
if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
if(any(s<=0))stop("s must be positive")
len <- max(length(q),length(m),length(s),length(size))
if(length(q)!=len){
	if(length(q)==1)q <- rep(q,len)
	else stop("length of q incorrect")}
if(length(size)!=len){
	if(length(size)==1)size <- rep(size,len)
	else stop("size must be the same length as q")}
if(any(q>size))stop("q must be <= size")
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("m and q must have the same length")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("s and q must have the same length")}
.C("pmb",
	as.integer(q),
	as.integer(size),
	as.double(m),
	as.double(s),
	as.integer(length(q)),
	res=double(length(q)),
	## DUP=FALSE,
	PACKAGE="rmutil")$res}

dmultbinom <- function(y, size, m, s, log=FALSE){
if(any(y<0))stop("y must contain non-negative values")
if(any(size<0))stop("size must contain non-negative values")
if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
if(any(s<=0))stop("s must be positive")
ly <- max(length(y),length(m),length(s),length(size))
if(length(y)!=ly){
	if(length(y)==1)y <- rep(y,ly)
	else stop("length of y incorrect")}
if(length(size)!=ly){
	if(length(size)==1)size <- rep(size,ly)
	else stop("size must be the same length as y")}
if(any(y>size))stop("y must be <= size")
if(length(m)!=ly){
	if(length(m)==1)m <- rep(m,ly)
	else stop("m and y must have the same length")}
if(length(s)!=ly){
	if(length(s)==1)s <- rep(s,ly)
	else stop("s and y must have the same length")}
tmp <- .C("dmb",
	as.integer(y),
	as.integer(size),
	as.double(m),
	as.double(s),
	as.integer(length(y)),
	as.double(rep(1,length(y))),
	res=double(length(y)),
	## DUP=FALSE,
	PACKAGE="rmutil")$res
if(!log)tmp <- exp(tmp)
tmp}

qmultbinom <- function(p, size, m, s){
h <- function(y).C("pmb",
	as.integer(y),
	as.integer(size[i]),
	as.double(m[i]),
	as.double(s[i]),
	as.integer(1),
	res=double(1),
	## DUP=FALSE,
	PACKAGE="rmutil")$res-p[i]
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0)||any(m>=1))stop("m must lie between 0 and 1")
if(any(s<=0))stop("s must be positive")
len <- max(length(p),length(m),length(s),length(size))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(size)!=len){
	if(length(size)==1)size <- rep(size,len)
	else stop("length of size incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(0,size[i])
	tmp[i] <- if(h(interval[1])*h(interval[2])>0)0
		else uniroot(h,interval)$root}
round(tmp)}

rmultbinom <- function(n=1, size, m, s)
	qmultbinom(runif(n),size=size,m=m,s=s)

### double Poisson distribution
###
pdoublepois <- function(q, m, s){
if(any(q<0))stop("q must contain non-negative values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
len <- max(length(q),length(m),length(s))
if(length(q)!=len){
	if(length(q)==1)q <- rep(q,len)
	else stop("length of q incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("m and q must have the same length")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("s and q must have the same length")}
.C("pdp",
	as.integer(q),
	as.integer(3*max(c(q,100))),
	as.double(m),
	as.double(s),
	as.integer(length(q)),
	res=double(length(q)),
	## DUP=FALSE,
	PACKAGE="rmutil")$res}

ddoublepois <- function(y, m, s, log=FALSE){
if(any(y<0))stop("y must contain non-negative values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
ly <- max(length(y),length(m),length(s))
if(length(y)!=ly){
	if(length(y)==1)y <- rep(y,ly)
	else stop("length of y incorrect")}
if(length(m)!=ly){
	if(length(m)==1)m <- rep(m,ly)
	else stop("m and y must have the same length")}
if(length(s)!=ly){
	if(length(s)==1)s <- rep(s,ly)
	else stop("s and y must have the same length")}
tmp <- .C("ddp",
	as.integer(y),
	as.integer(3*max(c(y,100))),
	as.double(m),
	as.double(s),
	as.integer(length(y)),
	as.double(rep(1,length(y))),
	res=double(length(y)),
	## DUP=FALSE,
	PACKAGE="rmutil")$res
if(!log)tmp <- exp(tmp)
tmp}

qdoublepois <- function(p, m, s){
h <- function(y).C("pdp",
	as.integer(y),
	as.integer(3*max(c(y,100))),
	as.double(m[i]),
	as.double(s[i]),
	as.integer(1),
	res=double(1),
	## DUP=FALSE,
	PACKAGE="rmutil")$res-p[i]
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<0))stop("m must be positive")
if(any(s<0))stop("s must be positive")
len <- max(length(p),length(m),length(s))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(0,20)
	if(h(interval[1])*h(interval[2])>0&&h(interval[1])>0)tmp[i] <- 0
	else {
		while(h(interval[1])*h(interval[2])>0)interval <- 2*interval
		tmp[i] <- uniroot(h,interval)$root}}
round(tmp)}

rdoublepois <- function(n=1, m, s) qdoublepois(runif(n),m=m,s=s)

### multiplicative Poisson distribution
###
pmultpois <- function(q, m, s){
if(any(q<0))stop("q must contain non-negative values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0|s>1))stop("s must be positive and <= 1")
len <- max(length(q),length(m),length(s))
if(length(q)!=len){
	if(length(q)==1)q <- rep(q,len)
	else stop("length of q incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("m and q must have the same length")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("s and q must have the same length")}
.C("pmp",
	as.integer(q),
	as.integer(3*max(c(q,100))),
	as.double(m),
	as.double(s),
	as.integer(length(q)),
	res=double(length(q)),
	## DUP=FALSE,
	PACKAGE="rmutil")$res}

dmultpois <- function(y, m, s, log=FALSE){
if(any(y<0))stop("y must contain non-negative values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0|s>1))stop("s must be positive and <= 1")
ly <- max(length(y),length(m),length(s))
if(length(y)!=ly){
	if(length(y)==1)y <- rep(y,ly)
	else stop("length of y incorrect")}
if(length(m)!=ly){
	if(length(m)==1)m <- rep(m,ly)
	else stop("m and y must have the same length")}
if(length(s)!=ly){
	if(length(s)==1)s <- rep(s,ly)
	else stop("s and y must have the same length")}
tmp <- .C("dmp",
	as.integer(y),
	as.integer(3*max(c(y,100))),
	as.double(m),
	as.double(s),
	as.integer(length(y)),
	as.double(rep(1,length(y))),
	res=double(length(y)),
	## DUP=FALSE,
	PACKAGE="rmutil")$res
if(!log)tmp <- exp(tmp)
tmp}

qmultpois <- function(p, m, s){
h <- function(y).C("pmp",
	as.integer(y),
	as.integer(3*max(c(y,100))),
	as.double(m[i]),
	as.double(s[i]),
	as.integer(1),
	res=double(1),
	## DUP=FALSE,
	PACKAGE="rmutil")$res-p[i]
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<0))stop("m must be positive")
if(any(s<=0|s>1))stop("s must be positive and <= 1")
len <- max(length(p),length(m),length(s))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(0,20)
	if(h(interval[1])*h(interval[2])>0&&h(interval[1])>0)tmp[i] <- 0
	else {
		while(h(interval[1])*h(interval[2])>0)interval <- 2*interval
		tmp[i] <- uniroot(h,interval)$root}}
round(tmp)}

rmultpois <- function(n=1, m, s) qmultpois(runif(n),m=m,s=s)

### power variance function Poisson distribution
###
ppvfpois <- function(q, m, s, f){
if(any(q<0))stop("q must contain non-negative values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
len <- max(length(q),length(m),length(s))
if(length(q)!=len){
	if(length(q)==1)q <- rep(q,len)
	else stop("length of q incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("m and q must have the same length")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("s and q must have the same length")}
if(length(f)!=len){
	if(length(f)==1)f <- rep(f,len)
	else stop("f and q must have the same length")}
.C("ppvfp",
	as.integer(q),
	as.double(m),
	as.double(s),
	as.double(f),
	as.integer(length(q)),
	res=double(length(q)),
	## DUP=FALSE,
	PACKAGE="rmutil")$res}

dpvfpois <- function(y, m, s, f, log=FALSE){
if(any(y<0))stop("y must contain non-negative values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
if(any(f>=1))stop("f must be < 1")
ly <- max(length(y),length(m),length(s))
if(length(y)!=ly){
	if(length(y)==1)y <- rep(y,ly)
	else stop("length of y incorrect")}
if(length(m)!=ly){
	if(length(m)==1)m <- rep(m,ly)
	else stop("m and y must have the same length")}
if(length(s)!=ly){
	if(length(s)==1)s <- rep(s,ly)
	else stop("s and y must have the same length")}
if(length(f)!=ly){
	if(length(f)==1)f <- rep(f,ly)
	else stop("f and y must have the same length")}
tmp <- log(.C("dpvfp",
	as.integer(y),
	as.double(m),
	as.double(s),
	as.double(f),
	as.integer(length(y)),
	as.double(rep(1,length(y))),
	res=double(length(y)),
	## DUP=FALSE,
	PACKAGE="rmutil")$res)
if(!log)tmp <- exp(tmp)
tmp}

qpvfpois <- function(p, m, s, f){
h <- function(y).C("ppvfp",
	as.integer(y),
	as.double(m[i]),
	as.double(s[i]),
	as.double(f[i]),
	as.integer(1),
	res=double(1),
	## DUP=FALSE,
	PACKAGE="rmutil")$res-p[i]
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<0))stop("m must be positive")
if(any(s<0))stop("s must be positive")
len <- max(length(p),length(m),length(s))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
if(length(f)!=len){
	if(length(f)==1)f <- rep(f,len)
	else stop("f and p must have the same length")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	if(f[i]==0)tmp[i] <- qnbinom(p[i],s[i],m[i]/(m[i]+s[i]))
	else {
		interval <- c(0,20)
		if(h(interval[1])*h(interval[2])>0&&h(interval[1])>0)
			tmp[i] <- 0
		else {
			while(h(interval[1])*h(interval[2])>0)
				interval <- 2*interval
			tmp[i] <- uniroot(h,interval)$root}}}
round(tmp)}

rpvfpois <- function(n=1, m, s, f) qpvfpois(runif(n),m=m,s=s,f=f)

### gamma count distribution
###
pgammacount <- function(q, m, s){
if(any(q<0))stop("q must contain non-negative values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
1-pgamma(m*s,(q+1)*s,1)}

dgammacount <- function(y, m, s, log=FALSE){
if(any(y<0))stop("y must contain non-negative values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
tmp <- ifelse(y==0,pgamma(m*s,(y+1)*s,1,log.p=TRUE,lower.tail=FALSE),
	log(pgamma(m*s,y*s+(y==0),1)-pgamma(m*s,(y+1)*s,1)))
if(!log)tmp <- exp(tmp)
tmp}

qgammacount <- function(p, m, s){
h <- function(y) 1-pgamma(m[i]*s[i],(y+1)*s[i],1)-p[i]
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0))stop("m must be positive")
if(any(s<0))stop("s must be positive")
len <- max(length(p),length(m),length(s))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(0,20)
	if(h(interval[1])*h(interval[2])>0&&h(interval[1])>0)tmp[i] <- 0
	else {
		while(h(interval[1])*h(interval[2])>0)interval <- 2*interval
		tmp[i] <- uniroot(h,interval)$root}}
round(tmp)}

rgammacount <- function(n=1, m, s) qgammacount(runif(n),m=m,s=s)

### Consul generalized Poisson distribution
###
pconsul <- function(q, m, s){
if(any(q<0))stop("q must contain non-negative values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
len <- max(length(q),length(m),length(s))
if(length(q)!=len){
	if(length(q)==1)q <- rep(q,len)
	else stop("length of q incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(any(s<0))stop("s must be positive")
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
res <- vector("numeric",length(q))
for(i in 1:length(q)){
	qq <- 0:q[i]
	res[i] <- sum(m[i]*exp(-(m[i]+qq*(s[i]-1))/s[i])*
		(m[i]+qq*(s[i]-1))^(qq-1)/(s[i]^qq*gamma(qq+1)))}
res}

dconsul <- function(y, m, s, log=FALSE){
if(any(y<0))stop("y must contain non-negative values")
if(any(m<=0))stop("m must be positive")
if(any(s<=0))stop("s must be positive")
tmp <- log(m)-(m+y*(s-1))/s+(y-1)*log(m+y*(s-1))-y*log(s)-lgamma(y+1)
if(!log)tmp <- exp(tmp)
tmp}

qconsul <- function(p, m, s){
h <- function(y) {
	pp <- 0:y
	sum(m[i]*exp(-(m[i]+pp*(s[i]-1))/s[i])*(m[i]+pp*(s[i]-1))^(pp-1)/
		(s[i]^pp*gamma(pp+1)))-p[i]}
if(any(p<0|p>1))stop("p must lie between 0 and 1")
if(any(m<=0))stop("m must be positive")
if(any(s<0))stop("s must be positive")
len <- max(length(p),length(m),length(s))
if(length(p)!=len){
	if(length(p)==1)p <- rep(p,len)
	else stop("length of p incorrect")}
if(length(m)!=len){
	if(length(m)==1)m <- rep(m,len)
	else stop("length of m incorrect")}
if(length(s)!=len){
	if(length(s)==1)s <- rep(s,len)
	else stop("length of s incorrect")}
tmp <- vector(mode="numeric",len)
for (i in 1:len){
	interval <- c(0,20)
	if(h(interval[1])*h(interval[2])>0&&h(interval[1])>0)tmp[i] <- 0
	else {
		while(h(interval[1])*h(interval[2])>0)interval <- 2*interval
		tmp[i] <- uniroot(h,interval)$root}}
round(tmp)}

rconsul <- function(n=1, m, s) qconsul(runif(n),m=m,s=s)
