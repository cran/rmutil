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
#     gettvc(response, times=NULL, tvcov=NULL, tvctimes=NULL,
#	oldtvcov=NULL, ties=TRUE)
#
#  DESCRIPTION
#
#    Function to find the most recent value of a time-varying
# covariate not recorded at the same time as the response.

gettvc <- function(response, times=NULL, tvcov=NULL, tvctimes=NULL,
	oldtvcov=NULL, ties=TRUE){
#
# if necessary, make into response object and remove NAs
#
if(inherits(response,"response"))izr <- response
else izr <- restovec(response,times)
if(any(is.na(izr$y))){
	isna <- TRUE
	nna <- length(izr$y)
	tmp <- NULL
	rna <- !is.na(izr$y)
	irna <- (1:nna)*rna
	j <- c(0,cumsum(nobs(izr)))
	for(i in 1:length(nobs(izr)))tmp <- c(tmp,sum(rna[(j[i]+1):j[i+1]]))
	zr <- list()
	zr$nobs <- tmp
	zr$times <- izr$times[rna]
	zr$y <- izr$y[rna]
	class(zr) <- "response"}
else {
	isna <- FALSE
	zr <- izr}
#
# remove NAs from time-varying covariate
#
zt <- restovec(tvcov,tvctimes)
if(any(is.na(zt$y))){
	tmp <- NULL
	rna <- !is.na(zt$y)
	j <- c(0,cumsum(nobs(zt)))
	for(i in 1:length(nobs(zt)))tmp <- c(tmp,sum(rna[(j[i]+1):j[i+1]]))
	zt$nobs <- tmp
	zt$times <- zt$times[rna]
	zt$y <- zt$y[rna]}
if(length(nobs(izr))!=length(nobs(zt)))
	stop("response and covariate do not have the same number of individuals")
#
# obtain new aligned times
#
nind <- length(nobs(zr))
nld <- max(c(nobs(zr),nobs(zt)))
z2 <- .Fortran("gettvc",
	x=as.double(zr$times),
	y=as.double(zr$y),
	xtvc=as.double(zt$times),
	tvcov=as.double(zt$y),
	nobs=as.integer(nobs(zr)),
	nind=as.integer(nind),
	nknt=as.integer(nobs(zt)),
	ties=as.logical(ties),
	xu=matrix(0,nrow=nind,ncol=2*nld),
	ndelta=logical(2*nld*nind),
	tvcov2=matrix(0,nrow=nind,ncol=2*nld),
	nu=integer(nind),
	wu=double(2*nld),
	nld=as.integer(nld),
	tvcov3=double(length(zr$y)),
	## DUP=FALSE,
	PACKAGE="rmutil")
if(isna){
	tvcov3 <- rep(NA,nna)
	tvcov3[irna] <- z2$tvcov3}
else tvcov3 <- z2$tvcov3
rm(z2)
#
# check if new covariate or to be combined with others
#
new <- missing(oldtvcov)
if(!new&!is.list(oldtvcov)){
	warning("oldtvcov must form a list - ignored")
	new <- TRUE}
cname <- paste(deparse(substitute(tvcov)))
if(new)oldtvcov <- vector(mode="list",nind)
else if(!inherits(oldtvcov,"tvcov")){
	if(length(oldtvcov)!=nind)
		stop(paste("Previous time-varying covariate list must have length",nind))
	else if(!is.null(colnames(oldtvcov[[1]])))
		cname <- c(colnames(oldtvcov[[1]]),cname)
	else cname <- NULL}
else if(inherits(oldtvcov,"tvcov"))cname <- c(colnames(oldtvcov$tvcov),cname)
if(inherits(oldtvcov,"tvcov")){
#
# combine tvcov objects
#
	oldtvcov$tvcov <- cbind(oldtvcov$tvcov,tvcov3)
	colnames(oldtvcov$tvcov) <- cname}
else {
#
# create a new tvcov object
#
	nm <- 0
	for(i in 1:nind){
		oldtvcov[[i]] <- cbind(oldtvcov[[i]],tvcov3[(nm+1):
			(nm+nobs(izr)[i])])
		nm <- nm+nobs(izr)[i]
		if(!is.null(cname))colnames(oldtvcov[[i]]) <- cname}
	oldtvcov <- tvctomat(oldtvcov)}
invisible(oldtvcov)}
