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
#     gettvc(response, times=NULL, tvcov=NULL, tvctimes=NULL,
#	oldtvcov=NULL, ties=TRUE)
#
#  DESCRIPTION
#
#    Function to find the most recent value of a time-varying
# covariate not recorded at the same time as the response.

gettvc <- function(response, times=NULL, tvcov=NULL, tvctimes=NULL,
	oldtvcov=NULL, ties=TRUE){
if(is.list(response)&&inherits(response,"response"))zr <- response
else zr <- restovec(response,times)
if(any(is.na(zr$y)))stop("NAs in response; matching cannot be done")
zt <- restovec(tvcov,tvctimes)
if(any(is.na(zt$y))){
	tmp <- NULL
	rna <- !is.na(zt$y)
	j <- c(0,cumsum(zt$nobs))
	for(i in 1:length(zt$nobs))tmp <- c(tmp,sum(rna[(j[i]+1):j[i+1]]))
	zt$nobs <- tmp
	zt$times <- zt$times[rna]
	zt$y <- zt$y[rna]}
if(length(zr$nobs)!=length(zt$nobs))stop("response and covariate do not have the same number of individuals")
nind <- length(zr$nobs)
nld <- max(c(zr$nobs,zt$nobs))
z2 <- .Fortran("gettvc",
	x=as.double(zr$times),
	y=as.double(zr$y),
	xtvc=as.double(zt$times),
	tvcov=as.double(zt$y),
	nobs=as.integer(zr$nobs),
	nind=as.integer(nind),
	nknt=as.integer(zt$nobs),
	ties=as.logical(ties),
	xu=matrix(0,nrow=nind,ncol=2*nld),
	ndelta=logical(2*nld*nind),
	tvcov2=matrix(0,nrow=nind,ncol=2*nld),
	nu=integer(nind),
	wu=double(2*nld),
	nld=as.integer(nld),
	tvcov3=double(length(zr$y)),
	DUP=F)
tvcov3 <- z2$tvcov3
z2 <- NULL
new <- missing(oldtvcov)
if(!new&!is.list(oldtvcov)){
	warning("oldtvcov must form a list - ignored")
	new <- T}
cname <- paste(deparse(substitute(tvcov)))
if(new)oldtvcov <- vector(mode="list",nind)
else if(!inherits(oldtvcov,"tvcov")){
	if(length(oldtvcov)!=nind)
		stop(paste("Previous time-varying covariate list must have length",nind))
	else if(!is.null(colnames(oldtvcov[[1]])))
		cname <- c(colnames(oldtvcov[[1]]),cname)
	else cname <- NULL}
else if(inherits(oldtvcov,"tvcov"))cname <- c(colnames(oldtvcov$tvcov),cname)
if(!inherits(oldtvcov,"tvcov")){
	nm <- 0
	for(i in 1:nind){
		oldtvcov[[i]] <- cbind(oldtvcov[[i]],tvcov3[(nm+1):(nm+zr$nobs[i])])
		nm <- nm+zr$nobs[i]}
	if(!is.null(cname))colnames(oldtvcov[[1]]) <- cname
	oldtvcov <- tvctomat(oldtvcov)}
else {
	oldtvcov$tvcov <- cbind(oldtvcov$tvcov,tvcov3)
	colnames(oldtvcov$tvcov) <- cname}
invisible(oldtvcov)}
