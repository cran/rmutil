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
#     read.list(file="", skip=0, nlines=2, order=NULL)
#     read.surv(file="", skip=0, nlines=1, cumulative=T, all=T)
#     read.rep(file, header=TRUE, skip=0, col.names=NULL,
#	response, id=NULL, times=NULL, censor=NULL,
#	totals=NULL, weights=NULL, nest=NULL, delta=NULL,
#	coordinates=NULL, type=NULL, ccov=NULL, tvcov=NULL, na.rm=T)
#
#  DESCRIPTION
#
#    Utility functions for reading repeated measurements data

### read unbalanced repeated measurements from a file into a list
###
read.list <- function(file="", skip=0, nlines=2, order=NULL){
#
# check if order of columns is to be changed
if(!is.null(order)){
	if(length(order)!=nlines)stop("order must have length",nlines,"\n")
	else if(any(range(order)!=c(1,nlines)))
		stop("order must have values in (",c(1,nlines),")\n")}
#
# scan in the data until EOF
#
continue <- TRUE
result <- list()
while(continue){
	x <- scan(file,skip=skip,nlines=nlines,quiet=TRUE)
	skip <- skip+nlines
	if(length(x)==0)continue <- FALSE
	else {
		tmp <- matrix(x,ncol=nlines)
		if(!is.null(order))tmp <- tmp[,order]
		result <- c(result,list(tmp))}}
invisible(result)}

### read unbalanced event history data from a file into a list
###
read.surv <- function(file="", skip=0, nlines=1, cumulative=TRUE, all=TRUE){
#
# scan in the data until EOF
#
continue <- TRUE
result <- list()
censor <- NULL
while(continue){
	x <- scan(file,skip=skip,nlines=nlines,quiet=TRUE)
	skip <- skip+nlines
	if(length(x)==0)continue <- FALSE
	else {
		# find response times (if all==TRUE, times and censor
		# indicator alternate on the line
		if(all)mm <- matrix(x,ncol=2,byrow=TRUE)[,1]
		else mm <- x[1:(length(x)-1)]
		# if cumulative times, find times between events
		if(cumulative)mm <- c(mm[1],diff(mm))
		result <- c(result,list(mm))
		# create vector of censor indicators for last time of
		# each individual only
		censor <- c(censor,x[length(x)])}}
invisible(list(result,censor))}

### read a rectangular data set from a file and create a repeated data object
###
read.rep <- function(file, header=TRUE, skip=0, sep = "",
	na.strings="NA", response, id=NULL, times=NULL, censor=NULL,
	totals=NULL, weights=NULL, nest=NULL, delta=NULL,
	coordinates=NULL, type=NULL, ccov=NULL, tvcov=NULL, na.rm=TRUE){
if(missing(response)||!is.character(response))
	stop("name(s) of response variables must be supplied")
if(missing(file)||!is.character(file))
	stop("a file name must be supplied")
dataframe <- read.table(file,header=header,skip=skip,
	na.strings=na.strings,sep=sep)
#
# find response information and construct object
#
cn <- colnames(dataframe)
nc <- match(response,cn)
if(any(is.na(nc)))stop(paste("response",response[is.na(nc)],"not found"))
if(!is.numeric(z <- as.matrix(dataframe[,nc,drop=FALSE])))
	stop("response must be numeric")
z <- list(response=list(y=z,nobs=NULL,times=NULL,nest=NULL,coordinates=NULL,
	censor=NULL,n=NULL,wt=NULL,delta=NULL,units=NULL,type=NULL),
	ccov=NULL,tvcov=NULL)
class(z) <- "repeated"
class(z$response) <- "response"
tobs <- dim(z$response$y)[1]
nrcol <- dim(z$response$y)[2]
if(is.null(type))z$response$type <- rep("unknown",nrcol)
else if(length(type)!=nrcol)stop("a type must be supplied for each response")
else {
	for(i in 1:length(type))
		z$response$type[i] <- match.arg(type[i],
			c("nominal","ordinal","discrete","duration","continuous","unknown"))
	if(any(is.na(z$response$type)))
		z$response$type[is.na(z$response$type)] <- "unknown"}
rna <- rep(TRUE,tobs)
for(i in 1:nrcol)rna <- rna&!is.na(z$response$y[,i])
if(is.null(id))
	z$response$nobs <- if(is.null(times)) rep(1,tobs) else tobs
else {
	if(!is.character(id)||length(id)>1)
		stop("id must be the name of one variable")
	nc <- match(id,cn)
	if(is.na(nc))stop("id not found")
	id <- as.vector(dataframe[,nc])
	if(is.character(id)||is.factor(id))id <- as.numeric(as.factor(id))
	else if(any(diff(id)!=0&diff(id)!=1,na.rm=TRUE))
		warning("id not consecutively numbered")
	nobs <- table(id)
	z$response$nobs <- as.vector(nobs[match(unique(id),names(nobs))])}
if(any(z$response$nobs!=1)&&length(z$response$nobs)>1)for(i in unique(id)){
	if(any(diff((1:tobs)[id==i])>1,na.rm=TRUE))
		stop(paste("observations for individual",i,"not together in table"))}
if(!is.null(nest)){
	if(all(z$response$nobs==1))
		stop("these are not repeated measurements - nest not valid")
	if(!is.character(nest)||length(nest)>1)
		stop("nest must be the name of one variable")
	nc <- match(nest,cn)
	if(is.na(nc))stop("nest not found")
	z$response$nest <- as.vector(dataframe[,nc])
	if(is.character(z$response$nest))
		z$response$nest <- as.numeric(as.factor(z$response$nest))
	else if(!is.numeric(z$response$nest))stop("nest must be numeric")
	rna <- rna&!is.na(z$response$nest)}
if(!is.null(times)){
	if(!is.character(times)||length(times)>1)
		stop("times must be the name of one variable")
	nc <- match(times,cn)
	if(is.na(nc))stop("times not found")
	z$response$times <- as.vector(dataframe[,nc])
	if(!is.numeric(z$response$times))stop("times must be numeric")
	rna <- rna&!is.na(z$response$times)}
if(!is.null(times)||!is.null(nest))for(i in unique(id)){
	if(!is.null(times)&&any(diff(z$response$times[id==i])<0,na.rm=TRUE))
		stop(paste("negative time step for individual",i))
	if(!is.null(nest)&&any(diff(z$response$nest[id==i])!=0&
		diff(z$response$nest[id==i])!=1,na.rm=TRUE))
		stop(paste("nest for individual",i,"not consecutively numbered"))}
if(!is.null(censor)){
	if(!is.character(censor)||length(censor)!=nrcol)
		stop("censor must have one name per response variable")
	nc <- match(censor,cn)
	if(any(is.na(nc)))stop("censor",censor[is.na(nc)],"not found")
	z$response$censor <- as.matrix(dataframe[,nc,drop=FALSE])
	if(!is.numeric(z$response$censor))stop("censor must be numeric")
	if(any(z$response$censor!=1&z$response$censor!=0&
		z$response$censor!=-1,na.rm=TRUE))
		stop("censor can only have values, -1, 0, 1")
	for(i in 1:nrcol)if(!all(is.na(z$response$censor[,i]))){
		rna <- rna&!is.na(z$response$censor[,i])
		if(z$response$type[i]=="unknown")
			z$response$type[i] <- "duration"}}
if(!is.null(totals)){
	if(!is.character(totals)||length(totals)!=nrcol)
		stop("totals must have one name per response variable")
	nc <- match(totals,cn)
	if(any(is.na(nc)))stop("totals",totals[is.na(nc)],"not found")
	z$response$n <- as.matrix(dataframe[,nc,drop=FALSE])
	if(!is.numeric(z$response$n))stop("totals must be numeric")
	if(any(z$response$y<0|z$response$n<z$response$y,na.rm=TRUE))
		stop("responses must lie between 0 and totals")
	for(i in 1:nrcol)if(!all(is.na(z$response$n[,i]))){
		rna <- rna&!is.na(z$response$n[,i])
		if(z$response$type[i]=="unknown")
			z$response$type[i] <- "nominal"}}
if(!is.null(delta)){
	if(is.numeric(delta)){
		if(length(delta)==1)
			z$response$delta <- matrix(delta,ncol=nrcol,nrow=tobs)
		else if(length(delta)==nrcol)
			z$response$delta <- matrix(rep(delta,tobs),ncol=nrcol,nrow=tobs,byrow=TRUE)
		else stop("delta must contain one value for each response")
		if(any(z$response$type=="unknown"))for(i in 1:nrcol)
			if(z$response$type[i]=="unknown")
				z$response$type[i] <- "continuous"}
	else {
		if(!is.character(delta)||length(delta)!=nrcol)
			stop("delta must have one name per response variable")
		nc <- match(delta,cn)
		if(any(is.na(nc)))stop("delta",delta[is.na(nc)],"not found")
		z$response$delta <- as.matrix(dataframe[,nc,drop=FALSE])
		if(!is.numeric(z$response$delta))stop("delta must be numeric")
		if(any(z$response$y<=0,na.rm=TRUE))
			stop("delta must be strictly positive")
		for(i in 1:nrcol)if(!all(is.na(z$response$delta[,i]))){
			rna <- rna&!is.na(z$response$delta[,i])
			if(z$response$type[i]=="unknown")
				z$response$type[i] <- "continuous"}}}
if(!is.null(weights)){
	if(!is.character(weights)||length(times)>1)
		stop("weights must be the name of one variable")
	nc <- match(weights,cn)
	if(is.na(nc))stop("weights not found")
	z$response$wt <- as.vector(dataframe[,nc])
	if(!is.numeric(z$response$wt))stop("weights must be numeric")
	rna <- rna&!is.na(z$response$wt)}
if(!is.null(coordinates)){
	if(!is.character(coordinates)||(length(coordinates)!=2&&
		length(coordinates)!=3))
		stop("coordinates must be the name of 2 or 3 variables")
	nc <- match(coordinates,cn)
	if(any(is.na(nc)))
		stop("coordinates",coordinates[is.na(nc)],"not found")
	z$response$coordinates <- as.matrix(dataframe[,nc,drop=FALSE])
	if(!is.numeric(z$response$coordinates))
		stop("coordinates must be numeric")
	for(i in 1:length(coordinates))
		rna <- rna&!is.na(z$response$coordinates[,i])}
#
# find time-varying covariates
#
if(!is.null(tvcov)){
	if(all(z$response$nobs==1))
		stop("these are not repeated measurements - tvcov not valid")
	z$tvcov <- list(tvcov=NULL,nobs=z$response$nobs)
	class(z$tvcov) <- "tvcov"
	nc <- match(tvcov,cn)
	if(any(is.na(nc)))stop("tvcov",tvcov[is.na(nc)],"not found")
	z$tvcov$tvcov <- dataframe[,nc,drop=FALSE]
	for(i in 1:length(tvcov))rna <- rna&!is.na(z$tvcov$tvcov[,i])
	# if no factor variables present, return a matrix anyway
	fac <- FALSE
	for(i in 1:dim(z$tvcov$tvcov)[2])
		if(!is.vector(z$tvcov$tvcov[,i],mode="numeric")){
		fac <- TRUE
		break}
	if(!fac)z$tvcov$tvcov <- as.matrix(z$tvcov$tvcov)}
#
# find time-constant covariates
#
if(!is.null(ccov)){
	z$ccov <- list(ccov=NULL)
	class(z$ccov) <- "tccov"
	nc <- match(ccov,cn)
	if(any(is.na(nc)))stop("ccov",ccov[is.na(nc)],"not found")
	z$ccov$ccov <- dataframe[,nc,drop=FALSE]
	for(i in unique(id))for(j in 1:length(ccov))
		if(sum(!is.na(unique(z$ccov$ccov[id==i,j])))>1)
		stop(paste("ccov",ccov[j],"for individual",i,"not constant"))
	for(i in 1:length(ccov))rna <- rna&!is.na(z$ccov$ccov[,i])
	j <- c(0,cumsum(z$response$nobs)[-length(z$response$nobs)])+1
	z$ccov$ccov <- z$ccov$ccov[j,,drop=FALSE]
	# if no factor variables present, return a matrix anyway
	fac <- FALSE
	for(i in 1:dim(z$ccov$ccov)[2])
		if(!is.vector(z$ccov$ccov[,i],mode="numeric")){
		fac <- TRUE
		break}
	if(!fac)z$ccov$ccov <- as.matrix(z$ccov$ccov)}
#
# remove NAs
#
if(na.rm&&any(!rna)){
	# remove NAs from variables associated with response
	z$response$y <- z$response$y[rna,,drop=FALSE]
	if(!is.null(z$response$times))z$response$times <- z$response$times[rna]
	if(!is.null(z$response$nest))z$response$nest <- z$response$nest[rna]
	if(!is.null(z$response$coordinates))
		z$response$coordinates <- z$response$coordinates[rna,]
	if(!is.null(z$response$n))z$response$n <- z$response$n[rna,,drop=FALSE]
	if(!is.null(z$response$censor)){
		z$response$censor <- z$response$censor[rna,,drop=FALSE]
		if(all(z$response$censor==1))z$response$censor <- NULL}
	if(!is.null(z$response$delta)&&length(z$response$delta)>1)
		z$response$delta <- z$response$delta[rna,,drop=FALSE]
	if(!is.null(z$tvcov))z$tvcov$tvcov <- z$tvcov$tvcov[rna,,drop=FALSE]
	# correct nobs
	tmp <- NULL
	j <- c(0,cumsum(z$response$nobs))
	for(i in 1:length(z$response$nobs)){
		tmp <- c(tmp,sum(rna[(j[i]+1):j[i+1]]))
		if(tmp[i]==0)
			warning(paste("Individual",i,"has no observations"))}
	z$response$nobs <- tmp[tmp>0]
	# remove NAs from ccov
	if(!is.null(z$ccov)){
		z$ccov$ccov <- z$ccov$ccov[tmp>0,,drop=FALSE]
		for(i in 1: dim(z$ccov$ccov)[2])
			if(length(unique(z$ccov$ccov[,i]))==1)
			warning(paste("covariate",colnames(z$ccov$ccov)[i],"has only one value\n"))}
	# remove NAs from tvcov
	if(!is.null(z$tvcov)){
		z$tvcov$nobs <- z$response$nobs
		for(i in 1: dim(z$tvcov$tvcov)[2])
			if(length(unique(z$tvcov$tvcov[,i]))==1)
			warning(paste("covariate",colnames(z$tvcov$tvcov)[i],"has only one value\n"))}}
if(!na.rm&&any(!rna))z$NAs <- !rna
#
# if independent observations, reset nobs
#
if(all(z$response$nobs==1))z$response$nobs <- 1
z}
