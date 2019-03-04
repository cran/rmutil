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
#     restovec(response, times=NULL, nest=NULL, coordinates=NULL,
#	censor=NULL, totals=NULL, weights=NULL, delta=NULL, type=NULL,
#	names=NULL, units=NULL, oldresponse=NULL, description=NULL)
#     tcctomat(ccov, names=NULL, units=NULL, oldtccov=NULL, dataframe=TRUE,
#     description=NULL)
#     tvctomat(tvcov, names=NULL, units=NULL, interaction=NULL, ccov=NULL,
#	oldtvcov=NULL, dataframe=TRUE, description=NULL)
#     rmna(response, tvcov=NULL, ccov=NULL)
#     lvna(response, tvcov=NULL, ccov=NULL)
#     dftorep(dataframe, response, id=NULL, times=NULL, censor=NULL,
#	totals=NULL, weights=NULL, nest=NULL, delta=NULL,
#	coordinates=NULL, type=NULL, ccov=NULL, tvcov=NULL, na.rm=TRUE)
#     as.data.frame(z)
#     as.matrix(z)
#     covariates(z, nind=NULL, names=NULL)
#     covind(z)
#     delta(z, nind=NULL, names=NULL)
#     formula(z)
#     names(z)
#     nesting(z, nind=NULL)
#     nobs(z)
#     plot.response(z, nind=NULL, name=NULL, nest=1, ccov=NULL, add=FALSE,
#	lty=NULL, pch=NULL, main=NULL, ylim=NULL, xlim=NULL, xlab=NULL,
#	ylab=NULL, ...)
#     plot.repeated(z, name=NULL, nind=NULL, ccov=NULL, add=FALSE, lty=NULL,
#	main=NULL, ylim=NULL, xlim=NULL, xlab=NULL, ylab=NULL, ...)
#     print(z)
#     response(z, nind=NULL, names=NULL)
#     times(z, nind=NULL)
#     transform(z, ...)
#     weights(z, nind=NULL)
#
#  DESCRIPTION
#
#    Utility functions for converting repeated measurements data to R objects

### function to create a response object
###
restovec <- function(response=NULL, times=NULL, nest=NULL, coordinates=NULL,
	censor=NULL, totals=NULL, weights=NULL, delta=NULL, type=NULL,
	names=NULL, units=NULL, oldresponse=NULL, description=NULL){
if(is.null(response))stop("A response must be supplied")
#
# check type
#
mvr <- length(names)>1
if(!is.null(type)){
	tmp <- NULL
	for(i in 1:length(type))
		tmp <- c(tmp,match.arg(type[i],
			c("nominal","ordinal","discrete","duration","continuous","unknown")))
	type <- tmp
	if(!mvr)mvr <- length(type)>1}
#
# initial values
#
nind <- 0
tnest <- nobs <- y <- NULL
#
# check if times are required
#
ttime <- !is.logical(times)||times
if((ttime&&is.logical(times))||!ttime)times <- NULL
#
# handle various forms of response
#
if(is.null(names))names <- paste(deparse(substitute(response)))
if(is.vector(response,mode="numeric")){
#
# numerical vector supplied: either univariate independent
#   observations or one time series
#
	y <- response
	# check if independent observations or one time series
	nobs <- if(is.null(times)) 1 else length(response)
	# check times
	if(is.vector(times,mode="numeric")){
		if(length(times)!=length(y))
			stop("times must be the same length as the response")}
	else if(!is.null(times))stop("times must be a vector")
	# check censor indicator
	if(is.vector(censor,mode="numeric")){
		if(length(censor)!=length(y)){
			if(length(censor)==1)
				censor <- c(rep(1,length(y)-1),censor)
			else stop("censor must be the same length as the response")}}
	else if(!is.null(censor))stop("censor must be a scalar or vector")
	# check coordinates
	if(!is.null(coordinates)&&(!is.matrix(coordinates)||(is.matrix(coordinates)&&length(dim(coordinates))!=2&&dim(coordinates)[2]!=2&&dim(coordinates)[1]!=length(y))))
		stop("coordinates must be a matrix with two columns and the same number of rows as the length of response")
	# check totals
	if(is.vector(totals,mode="numeric")){
		if(length(totals)!=length(y)){
			if(length(totals)==1)totals <- rep(totals,length(y))
			else stop("totals must be the same length as the response")}}
	else if(!is.null(totals)) stop("totals must be a vector")
	# check weights
	if(is.vector(weights,mode="numeric")){
		if(length(weights)!=length(y))
		stop("weights must be the same length as the response")}
	else if(!is.null(weights))stop("weights must be a vector")
	# check delta
	if(is.vector(delta,mode="numeric")){
		if(length(delta)!=length(y)){
			if(length(delta)==1)delta <- rep(delta,length(y))
			else stop("delta must be the same length as the response")}}
	else if(!is.null(delta))stop("delta must be a scalar or vector")}
else if(is.array(response)&&length(dim(response))==3){
#
# one multivariate 3-dim array of balanced repeated measurements supplied
#
	nobs <- rep(dim(response)[2],dim(response)[1])
	if(length(names)==1)names <- paste(names,1:dim(response)[3],sep="")
	# check times
	if(is.null(times)){
		if(is.null(censor)&&ttime)
			times <- as.double(rep(1:dim(response)[2],dim(response)[1]))}
	else if(is.matrix(times)){
		if(dim(times)[2]!=dim(response)[2]||dim(times)[1]!=dim(response)[1])
			stop("times matrix must be the same size as first two dimensions of array of responses")
		# time steps can be 0 but not negative
		for(i in 1:dim(response)[1])
			if(any(diff(times[i,])<0,na.rm=TRUE))
			stop(paste("negative time step for individual ",i))
		times <- as.vector(t(times))}
	else if(is.vector(times,mode="numeric")) {
		if(is.null(nest)&&any(diff(times)<0,na.rm=TRUE))
			stop("times must be increasing")
		if(length(times)!=dim(response)[2])
			stop("number of times must equal number of response columns")
		times <- rep(times,dim(response)[1])}
	else stop("times must be a vector or matrix")
	# check weights
	if(is.matrix(weights)){
		if(dim(weights)[1]!=dim(response)[1]||dim(weights)[2]!=dim(response)[2])
			stop("weights matrix must be the same size as first two dimensions of array of responses")
		else weights <- as.vector(t(weights))}
	else if(is.vector(weights,mode="numeric")){
		if(length(weights)!=dim(response)[1])
			stop("weights vector must have same length as number of individuals")
		else weights <- rep(weights,rep(dim(response)[2],dim(response)[1]))}
	else if(!is.null(weights))stop("weights must be a vector or matrix")
	# check nesting
	if(is.matrix(nest)){
		if(dim(weights)[1]!=dim(response)[1]||dim(weights)[2]!=dim(response)[2])
			stop("nest matrix must be the same size as first two dimensions of array of responses")
		for(i in 1:dim(nest)[1])
			if(any(diff(nest[i,])!=0&diff(nest[i,])!=1,na.rm=TRUE))
				stop("nest categories must be consecutive increasing integers")
		tnest <- as.vector(t(nest))}
	else if(is.vector(nest,mode="numeric")){
		if(length(nest)!=dim(response)[2])
			stop("nest vector must have same length as number of individuals")
		if(any(diff(nest)!=0&diff(nest)!=1,na.rm=TRUE))
			stop("nest categories must be consecutive increasing integers")
		tnest <- rep(nest,dim(response)[1])}
	else if(!is.null(nest))stop("nest must be a vector or matrix")
	# check censor indicator
	if(is.array(censor)&&length(dim(censor))==3){
		if(any(dim(censor)!=dim(response)))
			stop("censor array must be the same size as array of responses")
		tmp <- NULL
		for(i in 1:dim(response)[3])
			tmp <- cbind(tmp,as.vector(t(censor[,,i])))
		censor <- tmp
		rm(tmp)}
	else if(is.matrix(censor)){
	# if a matrix, corresponds to last observation of each subject
		if(dim(censor)[1]!=dim(response)[1]||dim(censor)[2]!=dim(response)[3])
			stop("censor matrix must have dimensions, number of individuals by number of variables")
		tmp <- array(1,dim(response))
		for(i in 1:dim(response)[3])
			tmp[,dim(tmp)[2],i] <- censor[,i]
		censor <- NULL
		for(i in 1:dim(response)[3])
			censor <- cbind(censor,as.vector(t(tmp[,,i])))
		rm(tmp)}
	else if(!is.null(censor))stop("censor must be a matrix or array")
	# check totals
	if(is.array(totals)&&length(dim(totals))==3){
		if(any(dim(totals)!=dim(response)))
			stop("totals array must be the same size as array of responses")
		tmp <- NULL
		for(i in 1:dim(response)[3])
			tmp <- cbind(tmp,as.vector(t(totals[,,i])))
		totals <- tmp
		rm(tmp)}
	else if(is.matrix(totals)){
		if(dim(totals)[1]!=dim(response)[1]||dim(totals)[2]!=dim(response)[3])
			stop("totals matrix must have dimensions, number of individuals by number of variables")
		for(i in 1:dim(response)[3])
		tmp <- NULL
		for(i in 1:dim(response)[3])
			tmp <- cbind(tmp,rep(totals[,i],rep(dim(response)[2],dim(response)[1])))
		totals <- tmp
		rm(tmp)}
	else if(!is.null(totals))stop("totals must be a matrix or array")
	# check delta
	if(is.array(delta)&&length(dim(delta))==3){
		if(any(dim(delta)!=dim(response)))
			stop("delta array must be the same size as array of responses")
		tmp <- NULL
		for(i in 1:dim(response)[3])
			tmp <- cbind(tmp,as.vector(t(delta[,,i])))
		delta <- tmp
		rm(tmp)}
	else if(is.matrix(delta)){
		if(dim(delta)[1]!=dim(response)[1]||dim(delta)[2]!=dim(response)[3])
			stop("delta matrix must have dimensions, number of individuals by number of variables")
		for(i in 1:dim(response)[3])
		tmp <- NULL
		for(i in 1:dim(response)[3])
			tmp <- cbind(tmp,rep(delta[,i],rep(dim(response)[2],dim(response)[1])))
		delta <- tmp
		rm(tmp)}
	else if(!is.null(delta))stop("delta must be a matrix or array")
	y <- NULL
	for(i in 1:dim(response)[3])y <- cbind(y,as.vector(t(response[,,i])))
	mvr <- TRUE}
else if(mvr&&(is.matrix(response)||is.data.frame(response))){
#
# multivariate matrix or dataframe of independent observations or time
#   series supplied
#
	y <- as.matrix(response)
	# check if independent observations or one time series
	nobs <- if(is.null(times)) 1 else dim(response)[1]
	# check times
	if(is.vector(times,mode="numeric")){
		if(length(times)!=dim(y)[1])
			stop("times must be the same length as the number of responses/variable")}
	else if(!is.null(times))stop("times must be a vector")
	# check weights
	if(is.vector(weights,mode="numeric")){
		if(length(weights)!=dim(y)[1])
		stop("weights must be the same length as the number of responses/variable")}
	else if(!is.null(weights))stop("weights must be a vector")
	# check censor indicator
	if(is.vector(censor,mode="numeric")){
		if(is.null(times))stop("censor must be a matrix")
		if(length(censor)!=dim(y)[2])
			stop("censor must be the same length as the number of variables")
		censor <- rbind(matrix(1,ncol=dim(y)[2],nrow=dim(y)[1]-1),censor)}
	else if(is.matrix(censor)){
		if(any(dim(censor)!=dim(y)))stop("censor and response must have the same dimensions")}
	else if(!is.null(censor))stop("censor must be a vector or matrix")
	# check totals
	if(is.vector(totals,mode="numeric")){
		if(length(totals)!=dim(y)[2])
			stop("totals must be the same length as the number of variables")
		totals <- matrix(rep(totals,dim(y)[1]),ncol=dim(y)[2],byrow=TRUE)}
	else if(is.matrix(totals)){
		if(any(dim(totals)!=dim(y)))stop("totals and response must have the same dimensions")}
	else if(!is.null(totals)) stop("totals must be a vector or matrix")
	# check delta
	if(is.vector(delta,mode="numeric")){
		if(length(delta)!=dim(y)[2]){
			if(length(delta)==1)
				delta <- matrix(delta,nrow=dim(y)[1],ncol=dim(y)[2])
			else if(length(delta)==dim(y)[2])
				delta <- matrix(rep(delta,dim(y)[1]),ncol=dim(y)[2],byrow=TRUE)
			else stop("delta must be the same length as the number of variables")}}
	else if(is.matrix(delta)){
		if(any(dim(delta)!=dim(y)))stop("delta and response must have the same dimensions")}
	else if(!is.null(delta))stop("delta must be a scalar, vector, or matrix")}
else if(is.matrix(response)||is.data.frame(response)){
#
# one matrix or dataframe with independent observations or several time
#   series supplied
#
	# transform to a matrix
	if(is.data.frame(response))response <- as.matrix(response)
	if(is.data.frame(totals))totals <- as.matrix(totals)
	if(is.data.frame(delta))delta <- as.matrix(delta)
	# create vector of observations/individual
	nobs <- rep(dim(response)[2],dim(response)[1])
	# check times
	if(is.null(times)){
		if(is.null(censor)&&ttime)
			times <- as.double(rep(1:dim(response)[2],dim(response)[1]))}
	else if(is.matrix(times)){
		if(any(dim(times)!=dim(response)))
			stop("times matrix must have the same dimensions as response")
		# time steps can be 0 but not negative
		for(i in 1:dim(response)[1])
			if(any(diff(times[i,])<0,na.rm=TRUE))
			stop(paste("negative time step for individual ",i))
		times <- as.vector(t(times))}
	else if(is.vector(times,mode="numeric")) {
		if(is.null(nest)&&any(diff(times)<0,na.rm=TRUE))
			stop("times must be increasing")
		if(length(times)!=dim(response)[2])
			stop("number of times must equal number of response columns")
		times <- rep(times,dim(response)[1])}
	else stop("times must be a vector or matrix")
	# check weights
	if(is.matrix(weights)){
		if(any(dim(weights)!=dim(response)))
			stop("weights matrix must have the same dimensions as response")
		else weights <- as.vector(t(weights))}
	else if(is.vector(weights,mode="numeric")){
		if(length(weights)!=dim(response)[1])
			stop("weights vector must have same length as number of individuals")
		else weights <- rep(weights,rep(dim(response)[2],dim(response)[1]))}
	else if(!is.null(weights))stop("weights must be a vector or matrix")
	# check nesting
	if(is.matrix(nest)){
		if(any(dim(nest)!=dim(response)))
			stop("nest matrix must have the same dimensions as response")
		for(i in 1:dim(nest)[1])
			if(any(diff(nest[i,])!=0&diff(nest[i,])!=1,na.rm=TRUE))
				stop("nest categories must be consecutive increasing integers")
		tnest <- as.vector(t(nest))}
	else if(is.vector(nest,mode="numeric")){
		if(length(nest)!=dim(response)[2])
			stop("nest vector must have same length as number of individuals")
		if(any(diff(nest)!=0&diff(nest)!=1,na.rm=TRUE))
			stop("nest categories must be consecutive increasing integers")
		tnest <- rep(nest,dim(response)[1])}
	else if(!is.null(nest))stop("nest must be a vector or matrix")
	# check censor indicator
	if(is.matrix(censor)){
		if(any(dim(censor)!=dim(response)))
			stop("censor matrix must have the same dimensions as response")
		censor <- as.vector(t(censor))}
	else if(is.vector(censor,mode="numeric")){
	# if a vector, corresponds to last observation of each subject
		if(length(censor)!=dim(response)[1])
			stop("censor must be the same length as the number of individuals")
		else {
			tmp <- matrix(1,nrow=dim(response)[1],ncol=dim(response)[2])
			tmp[,dim(tmp)[2]] <- censor
			censor <- tmp}}
	else if(!is.null(censor))stop("censor must be a vector or matrix")
	# check totals
	if(is.matrix(totals)){
		if(any(dim(totals)!=dim(response)))
			stop("totals matrix must have the same dimensions as response")
		totals <- as.vector(t(totals))}
	else if(is.vector(totals,mode="numeric")){
		if(length(totals)!=dim(response)[1])
			stop("totals vector must have same length as number of individuals")
		else totals <- rep(totals,rep(dim(response)[2],dim(response)[1]))}
	else if(!is.null(totals))stop("totals must be a vector or matrix")
	# check delta
	if(is.matrix(delta)){
		if(any(dim(delta)!=dim(response)))
			stop("delta matrix must have the same dimensions as response")
		delta <- as.vector(t(delta))}
	else if(is.vector(delta,mode="numeric")){
		if(length(delta)>1){
			if(length(delta)!=dim(response)[2])
				stop("delta vector must have same length as number of individuals")
			else delta <- rep(delta,dim(response)[1])}}
	else if(!is.null(delta))stop("delta must be a vector or matrix")
	y <- as.vector(t(response))}
else if(mvr&&is.list(response)){
#
# a list of unbalanced multivariate repeated measurements supplied
#
	# obtain responses and possibly times
	y <- times <- nobs <- NULL
	nind <- 0
	nmv <- dim(as.matrix(response[[1]]))[2]
	for(i in response){
		nind <- nind+1
		if(!is.matrix(i))stop("response must be a list of matrices")
		nobs <- c(nobs,dim(i)[1])
		if(dim(i)[2]!=nmv)
			stop("all matrices must have the same number of columns")
		if(ttime){
			y <- rbind(y,i[,1:(nmv-1)])
			if(is.null(nest)&&any(diff(i[,nmv])<0,na.rm=TRUE))
				stop(paste("negative time step for individual ",nind))
			times <- c(times,i[,nmv])}
		else y <- rbind(y,i)}
	if(ttime)nmv <- nmv-1
	if(!is.null(nest)){
	# obtain nesting
		nind <- 0
		nes <- NULL
		if(!is.list(nest))stop("nest must be a list")
		for(i in nest){
			if(!is.vector(i,mode="numeric")&&!(is.matrix(i)&&dim(i)[2]==1))
				stop("nest must be a list of vectors")
			if(any(diff(i)!=0&diff(i)!=1,na.rm=TRUE))
				stop("nest categories must be consecutive increasing integers")
			nind <- nind+1
			if(length(i)!=nobs[nind])
				stop(paste("nest for individual",nind,"should have length",nobs[nind]))
			nes <- c(nes,i)}
		tnest <- nes
		rm(nes)}
	if(!is.null(weights)){
	# obtain weights
		nind <- 0
		wt <- NULL
		if(!is.list(weights))stop("weights must be a list")
		for(i in weights){
			if(!is.vector(i,mode="numeric")&&!is.vector(i,mode="logical")&&!(is.matrix(i)&&dim(i)[2]==1))
				stop("weights must be a list of vectors")
			nind <- nind+1
			if(length(i)!=nobs[nind])
				stop(paste("weights for individual",nind,"should have length",nobs[nind]))
			wt <- c(wt,i)}
		weights <- wt
		rm(wt)}
	if(!is.null(censor)){
	# obtain censoring
		nind <- 0
		cen <- NULL
		if(!is.list(censor))stop("censor must be a list")
		for(i in censor){
			if(!is.matrix(i))
				stop("censor must be a list of matrices")
			nind <- nind+1
			if(any(dim(i)!=c(nobs[nind],nmv)))
				stop(paste("censor for individual",nind,"should be",nobs[nind],"x",nmv,"matrix"))
			cen <- rbind(cen,i)}
		censor <- cen
		rm(cen)}
	if(!is.null(totals)){
	# obtain totals
		nind <- 0
		tot <- NULL
		if(!is.list(totals))stop("totals must be a list")
		for(i in totals){
			if(!is.matrix(i))
				stop("totals must be a list of matrices")
			nind <- nind+1
			if(any(dim(i)!=c(nobs[nind],nmv)))
				stop(paste("totals for individual",nind,"should be",nobs[nind],"x",nmv,"matrix"))
			tot <- rbind(tot,i)}
		totals <- tot
		rm(tot)}
	if(!is.null(delta)){
	# obtain delta
		nind <- 0
		del <- NULL
		if(!is.list(delta))stop("delta must be a list")
		for(i in delta){
			if(!is.matrix(i))
				stop("delta must be a list of matrices")
			nind <- nind+1
			if(any(dim(i)!=c(nobs[nind],nmv)))
				stop(paste("delta for individual",nind,"should be",nobs[nind],"x",nmv,"matrix"))
			del <- rbind(del,i)}
		delta <- del
		rm(del)}}
else if(is.list(response)){
#
# a list of unbalanced repeated measurements supplied
#
	# check if nest, delta, and/or totals are supplied
	# separately
	delv <- !is.null(delta)
	totv <- !is.null(totals)
	nestv <- !is.null(nest)
	if(is.null(censor)){
		# initialize
		times <- NULL
		tot <- del <- cen <- nes <- 0
		ncols <- dim(as.matrix(response[[1]]))[2]
		nc <- ttime+1
		if(ncols<nc)stop("matrices must have at least 2 columns: responses and times")
		else if(ncols>nc)for(j in response){
			# find columns for censor, totals, delta, nest
			j <- as.matrix(j)
			for(k in (nc+1):ncols){
				if(any(j[,k]<=0,na.rm=TRUE))cen <- k
				else if(any(j[,k]>1,na.rm=TRUE)&&all(j[,k]==trunc(j[,k]),na.rm=TRUE)){
					if(all(j[,1]>=0,na.rm=TRUE)&&all(j[,1]==trunc(j[,1]),na.rm=TRUE)&&all(j[,k]>=j[,1],na.rm=TRUE)&&any(diff(j[,k])<0,na.rm=TRUE)&&!totv)tot <- k
					else if(all(diff(j[,k])==0|diff(j[,k])==1,na.rm=TRUE)&&!nestv)nes <- k
					else stop(paste("column",k,"of unknown type"))}
				else if(!delv&&all(j[,k]>0,na.rm=TRUE))del <- k}
			if((((ncols==2&&!ttime)||ncols==3)&&(nes>0||cen>0||del>0||tot>0))||(((ncols==3&&!ttime)||ncols==4)&&((nes>0&&cen>0)||(nes>0&&del>0)||(nes>0&&tot>0)||(cen>0&&del>0)))||(((ncols>=4&&!ttime)||ncols>=5)&&((nes>0&&cen>0&&del>0)||(nes>0&&tot>0))))break}
		for(i in response){
			# obtain responses and times
			i <- as.matrix(i)
			nind <- nind+1
			if(dim(i)[2]!=ncols)
				stop(paste("Individual ",nind,"does not have a",ncols,"column matrix"))
			if(!nestv&&nes==0&&ttime&&any(diff(i[,2])<0,na.rm=TRUE))
				stop(paste("negative time step for individual ",nind))
			nobs <- c(nobs,dim(i)[1])
			y <- c(y,i[,1])
			if(ttime)times <- c(times,i[,2])
			# obtain censor, totals, delta, nest
			# from their columns
			if(nes>0){
				if(any(diff(i[,nes])!=0&diff(i[,nes])!=1,na.rm=TRUE))
					stop("nest categories for individual ",nind,"are not consecutive increasing integers")
				tnest <- c(tnest,i[,nes])}
			if(!totv&&tot>0)totals <- c(totals,i[,tot])
			if(cen>0)censor <- c(censor,i[,cen])
			if(!delv&&del>0)delta <- c(delta,i[,del])}}
	else if(is.vector(censor,mode="numeric")){
		# initialize (cannot be totals and censor available)
		del <- nes <- 0
		ncols <- dim(as.matrix(response[[1]]))[2]
		if(ncols>1){
			for(j in response){
			# find columns for delta, nest
				j <- as.matrix(j)
				for(k in 2:ncols){
					if(is.null(censor)&&any(j[,k]<=0,na.rm=TRUE))cen <- k
					else if(any(j[,k]>1,na.rm=TRUE)&&all(j[,k]==trunc(j[,k]),na.rm=TRUE))nes <- k
					else if(is.null(delta)&&all(j[,k]>0,na.rm=TRUE))del <- k}
				if((ncols==3&&(nes>0||cen>0||del>0))||(ncols==4&&((nes>0&&cen>0)||(nes>0&&del>0)||(cen>0&&del>0)))||(ncols>=5&&((nes>0&&cen>0&&del>0))))break}
			tmp <- NULL
			j <- 0
			for(i in response){
			# obtain responses
				i <- as.matrix(i)
				nind <- nind+1
				if(dim(i)[2]!=ncols)stop(paste("Individual ",nind,"does not have a",ncols,"column matrix"))
				nobs <- c(nobs,dim(i)[1])
				y <- c(y,i[,1])
				# construct censor
				tmp <- c(tmp,rep(1,dim(i)[1]-1),censor[j <- j+1])
				# obtain delta, nest from their columns
				if(nes>0){
					if(any(diff(i[,nes])!=0&diff(i[,nes])!=1,na.rm=TRUE))
					stop("nest categories for individual ",nind,"are not consecutive increasing integers")
					tnest <- c(tnest,i[,nes])}
				if(del>0)delta <- c(delta,i[,del])}
			censor <- tmp}
		else {
			tmp <- NULL
			j <- 0
			for(i in response){
			# obtain responses
				nind <- nind+1
				if(!is.vector(i,mode="numeric")&&!(is.matrix(i)&&dim(i)[2]==1))
					stop(paste("Individual ",nind,"does not have a vector or one column matrix"))
				# construct censor
				tmp <- c(tmp,rep(1,length(i)-1),censor[j <- j+1])
				y <- c(y,i)
				nobs <- c(nobs,length(i))}
			censor <- tmp}}
	else stop("If response is a list, censor must be a vector")
	if(nestv){
	# obtain nesting
		nind <- 0
		nes <- NULL
		if(!is.list(nest))stop("nest must be a list")
		for(i in nest){
			if(!is.vector(i,mode="numeric")&&!(is.matrix(i)&&dim(i)[2]==1))
				stop("nest must be a list of vectors")
			if(any(diff(i)!=0&diff(i)!=1,na.rm=TRUE))
				stop("nest categories must be consecutive increasing integers")
			nind <- nind+1
			if(length(i)!=nobs[nind])
				stop(paste("nest for individual",nind,"should have length",nobs[nind]))
			nes <- c(nes,i)}
		tnest <- nes
		rm(nes)}
	if(totv){
	# obtain totals
		nind <- 0
		tot <- NULL
		if(!is.list(totals))stop("totals must be a list")
		for(i in totals){
			if(!is.vector(i,mode="numeric")&&!(is.matrix(i)&&dim(i)[2]==1))
				stop("totals must be a list of vectors")
			nind <- nind+1
			if(length(i)!=nobs[nind])
				stop(paste("totals for individual",nind,"should have length",nobs[nind]))
			tot <- rbind(tot,i)}
		totals <- tot
		rm(tot)}
	if(delv){
	# obtain delta
		if(is.vector(delta)&&length(delta)==1)
			delta <-  rep(delta,length(y))
		else {
			nind <- 0
			del <- NULL
			if(!is.list(delta))stop("delta must be a list")
			for(i in delta){
				if(!is.vector(i,mode="numeric")&&!(is.matrix(i)&&dim(i)[2]==1))
					stop("delta must be a list of vectors")
				nind <- nind+1
				if(length(i)!=nobs[nind])
					stop(paste("delta for individual",nind,"should have length",nobs[nind]))
				del <- c(del,i)}
			delta <- del
			rm(del)}}
	if(!is.null(weights)){
	# obtain weights
		nind <- 0
		wt <- NULL
		if(!is.list(weights))stop("weights must be a list")
		for(i in weights){
			if(!is.vector(i,mode="numeric")&&!is.vector(i,mode="logical")&&!(is.matrix(i)&&dim(i)[2]==1))
				stop("weights must be a list of vectors")
			nind <- nind+1
			if(length(i)!=nobs[nind])
				stop(paste("weights for individual",nind,"should have length",nobs[nind]))
			wt <- c(wt,i)}
		weights <- wt
		rm(wt)}
	# check that totals, delta, and weights are now
	# correct if supplied separately
	if(totv){
		if(length(totals)==1)totals <- rep(totals,length(y))
		else if(length(totals)==length(nobs))
			totals <- totals[rep(1:length(nobs),nobs)]
		else if(length(totals)!=length(y))
			stop("totals must have one value per response")}
	if(delv&&length(delta)>1&&length(delta)!=length(y))
		stop("delta must have one value per response")
	if(!is.null(weights)&&length(y)!=length(weights))stop("weights must have one value per response")}
else stop("Responses must be supplied as a vector, matrix, dataframe, array, or list of matrices")
#
# make sure that the response is a matrix
#
if(!is.numeric(y))stop("response must be numeric")
if(is.matrix(y)){
	if(is.null(colnames(y)))colnames(y) <- if(length(names)==1)
		paste(names,1:dim(y)[2],sep="")
		else if(length(names)==dim(y)[2]) names
		else stop("incorrect number of names")}
else {
	y <- matrix(y,ncol=1)
	colnames(y) <- if(length(names)==1)names else "y"}
rownames(y) <- 1:dim(y)[1]
if(!is.null(units)){
	if(!is.character(units))stop("units must be a character vector")
	if(length(units)!=dim(y)[2])
		stop("units must be given for all responses")}
#
# check that type has the right length
#
if(is.null(type))type <- rep("unknown",dim(y)[2])
else {
	if(length(type)!=dim(y)[2])
		stop("a type must be supplied for each response")
	if(length(type)==1)type <- rep(type,dim(y)[2])
	if(any(is.na(type)))type[is.na(type)] <- "unknown"}
if(!is.null(censor)){
#
# check that censor has correct values
#
	if(!is.numeric(censor))stop("censor must be numeric")
	if(any(censor!=-1&censor!=0&censor!=1,na.rm=TRUE))
		stop("censor must only contain -1, 0, and 1")
	if(!is.null(censor)&&!is.matrix(censor))censor <- matrix(censor,ncol=1)
	# construct (cumulative) times from responses if univariate
	if(!mvr&&length(nobs)>1&&is.null(times)&&type=="duration"){
		j <- 1
		na <- is.na(y[,1])
		y[na,1] <- 0
		for(i in 1:length(nobs)){
			times <- c(times,cumsum(y[j:(j+nobs[i]-1),1]))
			j <- j+nobs[i]}
		y[na,1] <- NA}
	for(i in 1:dim(y)[2])
		if(type[i]=="unknown"&&any(!is.na(censor[,i])))
			type[i] <- "duration"
	# remove censor if unnecessary
	if(all(censor==1,na.rm=TRUE))censor <- NULL}
if(!is.null(totals)){
#
# check that totals has correct values
#
	if(!is.numeric(totals))stop("totals must be numeric")
	if(!is.matrix(totals))totals <- matrix(totals,ncol=1)
	for(i in 1:dim(y)[2])if(any(!is.na(totals[,i]))){
		if(any(y[,i]<0,na.rm=TRUE))
			stop("all responses must be non-negative for binomial data")
		if(any(totals[,i]<y[,i],na.rm=TRUE)||any(totals[,i]<0,na.rm=TRUE))
			stop("all totals must be non-negative and >= to responses")
		if(type[i]=="unknown"&&any(!is.na(totals[,i])))
			type[i] <- "nominal"}}
#
# check that nest, delta, weights have correct values
#
if(!is.null(tnest)){
	if(!is.numeric(tnest))stop("tnest must be numeric")
	if((any(tnest<1,na.rm=TRUE)||any(tnest!=trunc(tnest),na.rm=TRUE)))
		stop("nest must contain integers starting at 1")}
if(!is.null(delta)){
	if(!is.numeric(delta))stop("delta must be numeric")
	if(any(delta<=0,na.rm=TRUE))stop("delta must be strictly positive")
	if(!is.matrix(delta))delta <- matrix(delta,ncol=1)
	for(i in 1:dim(y)[2])
		if(type[i]=="unknown"&&any(!is.na(delta[,i])))
			type[i] <- "continuous"}
if(!is.null(weights)){
	if(!is.numeric(weights)&&!is.logical(weights))stop("weights must be numeric or logical")
	if(any(weights<0,na.rm=TRUE))stop("weights must be non-negative")}
#
# check that ordinal data have correct values
#
if(!is.na(match("ordinal",type))){
	for(i in 1:dim(y)[2])if(!is.na(match("ordinal",type[i]))){
		if(min(y[,i],na.rm=TRUE)!=0)
			stop("ordinal responses must be indexed from zero")
		if(any(as.integer(y[!is.na(y[,i]),i])!=y[!is.na(y[,i]),i]))
			stop("ordinal responses must be integers")}}
#
# combine with old object, if provided
#
if(!is.null(oldresponse)){
	if(any(nobs!=oldresponse$nobs))
		stop("old and new objects do not have the same number of observations per individual")
	if(dim(y)[1]!=dim(oldresponse$y)[1])
		stop("old and new objects do not have the same number of observations")
	if((!is.null(times)&&(is.null(oldresponse$times)||any(times!=oldresponse$times,na.rm=TRUE)))||(is.null(times)&&!is.null(oldresponse$times)))
		stop("old and new objects do not have the same times")
	if((!is.null(tnest)&&(is.null(oldresponse$tnest)||any(tnest!=oldresponse$tnest,na.rm=TRUE)))||(is.null(tnest)&&!is.null(oldresponse$tnest)))
		stop("old and new objects do not have the same nesting")
	if((!is.null(weights)&&(is.null(oldresponse$wt)||any(weights!=oldresponse$wt,na.rm=TRUE)))||(is.null(weights)&&!is.null(oldresponse$wt)))
		stop("old and new objects do not have the same weights")
	y <- cbind(oldresponse$y,y)
	if(!is.null(censor)){
		if(is.null(oldresponse$censor))
			censor <- cbind(matrix(NA,nrow=dim(oldresponse$y)[1],ncol=dim(oldresponse$y)[2]),censor)
		else censor <- cbind(oldresponse$censor,censor)}
	else if(!is.null(oldresponse$censor))
		censor <- cbind(oldresponse$censor,matrix(NA,nrow=dim(y)[1],ncol=dim(y)[2]))
	if(!is.null(totals)){
		if(is.null(oldresponse$n))
			totals <- cbind(matrix(NA,nrow=dim(oldresponse$y)[1],ncol=dim(oldresponse$y)[2]),totals)
		else totals <- cbind(oldresponse$n,totals)}
	else if(!is.null(oldresponse$n))
		totals <- cbind(oldresponse$n,matrix(NA,nrow=dim(y)[1],ncol=dim(y)[2]))
	if(!is.null(delta)){
		if(is.null(oldresponse$delta))
			delta <- cbind(matrix(1,nrow=dim(oldresponse$y)[1],ncol=dim(oldresponse$y)[2]),delta)
		else delta <- cbind(oldresponse$delta,delta)}
	else if(!is.null(oldresponse$delta))
		delta <- cbind(oldresponse$delta,matrix(1,nrow=dim(y)[1],ncol=dim(y)[2]))
	if(!is.null(units)){
		if(is.null(oldresponse$units))
			units <- c(rep(NA,dim(oldresponse$y)[2]),units)
		else units <- c(oldresponse$units,units)}
	else if(!is.null(oldresponse$units))
		units <- c(oldresponse$units,rep(NA,dim(y)[2]))
	type <- c(oldresponse$type,type)}
#
# put names on matrices
#
if(!is.null(totals)&&is.null(colnames(totals)))colnames(totals) <- colnames(y)
if(!is.null(censor)&&is.null(colnames(censor)))colnames(censor) <- colnames(y)
if(!is.null(delta)&&is.null(colnames(delta)))colnames(delta) <- colnames(y)
if(!is.null(units))names(units) <- colnames(y)
names(type) <- colnames(y)
#
# check for variable descriptions
#
if(!is.null(description)){
	if(!is.list(description))stop("description must be a list")
	if(!all(tmp <- names(description)%in%colnames(y)))
		stop(paste("variable(s)",names(description)[!tmp],"not found"))
	for(i in description)if(!is.character(i))
		stop("description list must contain character vectors")
	if(!is.null(oldresponse))
		description <- c(oldresponse$description,description)}
z <- list(y=y, nobs=nobs, times=times, nest=tnest, censor=censor, n=totals,
	coordinates=coordinates, wt=weights, delta=delta, units=units,
	type=type, description=description)
class(z) <- "response"
z}

### function to create a time-constant covariate (tccov) object
###
tcctomat <- function(ccov, names=NULL, units=NULL, oldccov=NULL,
	dataframe=TRUE, description=NULL){
if(inherits(ccov,"tccov")&&inherits(oldccov,"tccov")){
#
# check for compatibility
#
	if(dim(ccov$ccov)[1]!=dim(oldccov$ccov)[1])
		stop("incompatible tccov objects")
	if(!is.null(oldccov$units)||!is.null(ccov$units)){
		if(is.null(oldccov$units))
			oldccov$units <- rep(NA,dim(oldccov$ccov)[2])
		if(is.null(ccov$units))ccov$units <- rep(NA,dim(ccov$ccov)[2])
		oldccov$units <- c(oldccov$units,ccov$units)}
	oldccov$ccov <- cbind(oldccov$ccov,ccov$ccov)
	oldccov$linear <- NULL
	return(oldccov)}
linear <- NULL
if(is.language(ccov)){
#
# if provided as a formula, transform to a matrix
#
	linear <- ccov
	mt <- terms(ccov)
	mf <- model.frame(mt,parent.frame(),na.action=NULL)
	ccov <- model.matrix(mt,mf)[,-1,drop=FALSE]}
else if(is.factor(ccov)||is.vector(ccov,mode="character")){
#
# if provided as factor variables, make a dataframe
#
	if(is.null(names))names <- paste(deparse(substitute(ccov)))
	ccov <- data.frame(ccov)
	colnames(ccov) <- names}
if(is.vector(ccov,mode="numeric")){
#
# if a vector, get a name
#
	if(is.null(names))names <- paste(deparse(substitute(ccov)))
	ccov <- matrix(ccov,ncol=1)}
else if(is.data.frame(ccov)){
	if(!dataframe){
	# if a dataframe supplied, but should not be stored as one
                rm(names)
                units2 <- mt <- tmp3 <- tmp <- NULL
                j <- 0
                for(i in ccov){
                	j <- j+1
                	if(is.vector(i,mode="numeric")){
			# handle numeric vectors
                		tmp2 <- as.matrix(i)
				units2 <- c(units2,units[j])
                		colnames(tmp2) <- names(ccov)[j]}
                	else {
			# handle factor variables
                		mt <- terms(~i)
                		tmp2 <- model.matrix(mt,model.frame(mt,na.action=NULL))[,-1,drop=FALSE]
				units2 <- c(units2,rep(units[j],length(levels(i)[-1])))
				tmp3 <- dimnames(get(getOption("contrasts")[[if(is.ordered(i))2 else 1]])(levels(i),contrasts =TRUE))[[2]]
				if(is.null(tmp3))tmp3 <- 1:(length(levels(i))-1)
                		colnames(tmp2) <- paste(names(ccov)[j],tmp3,sep="")}
#                		colnames(tmp2) <- paste(names(ccov)[j],levels(i)[-1],sep="")}
                	tmp <- cbind(tmp,tmp2)}
		units <- units2
                ccov <- tmp
                rm(tmp,tmp2,tmp3,mt)}}
else if(!is.matrix(ccov))
	stop("Inter-unit (time-constant) covariates must be a vector, matrix, dataframe, or model formula")
if(is.null(colnames(ccov))){
#
# create names if not available
#
	if(is.null(names))names <- paste(deparse(substitute(ccov)))
	if(length(names)==1&&dim(ccov)[2]>1)
		names <- paste(names,1:dim(ccov)[2],sep="")
	colnames(ccov) <- names}
#
# check units
#
if(!is.null(units)){
	if(!is.character(units))stop("units must be a character vector")
	if(length(units)!=dim(ccov)[2])
		stop("units must be supplied for all covariates")}
#
# check for variable descriptions
#
if(!is.null(description)){
	if(!is.list(description))stop("description must be a list")
	if(!all(tmp <- names(description)%in%colnames(ccov)))
		stop(paste("variable(s)",names(description)[!tmp],"not found"))
	for(i in description)if(!is.character(i))
		stop("description list must contain character vectors")}
if(!is.null(oldccov)){
#
# combine new data with old data when available
#
	if(!inherits(oldccov,"tccov"))
		stop("oldccov must have class, tccov")
	oldccov$description <- c(oldccov$description,description)
	if(!is.null(oldccov$units)||!is.null(units)){
		if(is.null(oldccov$units))
			oldccov$units <- rep(NA,dim(oldccov$ccov)[2])
		if(is.null(units))units <- rep(NA,dim(ccov)[2])
		oldccov$units <- c(oldccov$units,units)}
	else if(dim(oldccov$ccov)[1]==dim(ccov)[1]){
		if(dataframe)oldccov$ccov <- data.frame(oldccov$ccov,ccov)
		else oldccov$ccov <- cbind(oldccov$ccov,ccov)}
	else stop("old and new covariates do not have the same number of individuals")}
else {
	if(dataframe)ccov <- as.data.frame(ccov)
	oldccov <- list(ccov=ccov,linear=linear,units=units,description=description)
	class(oldccov) <- "tccov"}
if(!is.null(oldccov$units))names(oldccov$units) <- colnames(oldccov$ccov)
if(is.data.frame(oldccov$ccov)){
#
# if no factor variables, store as a matrix anyway
#
	fac <- FALSE
	for(i in 1:dim(oldccov$ccov)[2])if(!is.vector(oldccov$ccov[,i],mode="numeric")){
		fac <- TRUE
		break}
	if(!fac)oldccov$ccov <- as.matrix(oldccov$ccov)}
oldccov}

### function to create a time-varying covariate (tvcov) object
###
tvctomat <- function(tvcov, names=NULL, units=NULL, interaction=NULL,
	ccov=NULL, oldtvcov=NULL, dataframe=TRUE, description=NULL){
#
# check for compatibility
#
if(inherits(tvcov,"tvcov")&&inherits(oldtvcov,"tvcov")){
	if(length(tvcov$nobs)!=length(oldtvcov$nobs)||
		any(tvcov$nobs!=oldtvcov$nobs))
		stop("incompatible tvcov objects")
	if(!is.null(oldtvcov$units)||!is.null(tvcov$units)){
		if(is.null(oldtvcov$units))
			oldtvcov$units <- rep(NA,dim(oldtvcov$tvcov)[2])
		if(is.null(tvcov$units))
			tvcov$units <- rep(NA,dim(tvcov$tvcov)[2])
		oldtvcov$units <- c(oldtvcov$units,tvcov$units)}
	oldtvcov$tvcov <- cbind(oldtvcov$tvcov,tvcov$tvcov)
	return(oldtvcov)}
nbs <- tvcv <- NULL
if(is.data.frame(tvcov)){
	if(is.null(names))names <- paste(deparse(substitute(tvcov)))
	if(length(names)!=1)stop("too many names")
	if(dataframe){
	# make new one-column dataframe
		nbs <- rep(dim(tvcov)[2],dim(tvcov)[1])
		tvcv <- as.data.frame(as.vector(t(as.matrix(tvcov))))
		colnames(tvcv) <- names}
	# if factors, as.matrix transforms to character for next step
	else tvcov <- as.matrix(tvcov)}
if(is.matrix(tvcov)&&is.character(tvcov)){
	nbs <- rep(dim(tvcov)[2],dim(tvcov)[1])
	if(is.null(names))names <- paste(deparse(substitute(tvcov)))
	if(length(names)!=1)stop("too many names")
	tvcv <- as.factor(as.vector(t(tvcov)))
	if(dataframe){
	# make new one-column dataframe
		tvcv <- as.data.frame(as.vector(t(as.matrix(tvcov))))
		colnames(tvcv) <- names}
	else {
	# make indicator matrix from factor
		mt <- terms(~tvcv)
		tmp3 <- dimnames(get(getOption("contrasts")[[if(is.ordered(tvcv))2 else 1]])(levels(tvcv),contrasts =TRUE))[[2]]
		if(is.null(tmp3))tmp3 <- 1:(length(levels(tvcv))-1)
		names <- paste(names,tmp3,sep="")
#		names <- paste(names,levels(tvcv)[-1],sep="")
		tvcv <- model.matrix(mt,model.frame(mt,na.action=NULL))[,-1,drop=FALSE]
		colnames(tvcv) <- names
		if(!is.null(units))units <- rep(units,dim(tvcv)[2])}}
else if(is.matrix(tvcov)){
#
# transform to a vector
#
	nbs <- rep(dim(tvcov)[2],dim(tvcov)[1])
	tvcv <- matrix(as.vector(t(tvcov)),ncol=1)
	if(!is.null(names)&&length(names)!=1)stop("too many names")
	colnames(tvcv) <- if(!is.null(names)) names
		else paste(deparse(substitute(tvcov)))}
else if(is.list(tvcov)&&!is.data.frame(tvcov)){
	if(inherits(tvcov,"tvcov")){
		nbs <- tvcov$nobs
		tvcv <- tvcov$tvcov}
	else {
		ncols <- dim(as.data.frame(tvcov[[1]]))[2]
		# create names
		if(is.null(names)){
			if(is.null(colnames(tvcov[[1]]))){
				names <- if(is.matrix(tvcov[[1]]))paste(paste(deparse(substitute(tvcov))),1:ncols,sep="")
				else paste(deparse(substitute(tvcov)))}
			else names <- colnames(tvcov[[1]])}
		# concatenate list elements into a matrix
		ff <- TRUE
		for(i in tvcov){
		# create one big dataframe
			i <- as.data.frame(i)
			if(dim(i)[2]!=ncols)
				stop("all elements of the list must have the same number of columns")
			nbs <- c(nbs,dim(i)[1])
			if(ff){
				tvcv <- i
				ff <- FALSE}
			else tvcv <- rbind(tvcv,i)}
		if(is.character(tvcv)||!dataframe){
		# if necessary, transform to indicator matrix
			tmp <- tmp2 <- mt <- NULL
			for(i in 1:ncols){
				if(is.numeric(tvcv[,i])){
					tmp <- cbind(tmp,tvcv[,i])
					tmp2 <- c(tmp2,names[i])}
				else {
					mt <- terms(~tvcv[,i])
					tmp <- cbind(tmp,model.matrix(mt,model.frame(mt,na.action=NULL))[,-1,drop=FALSE])
					tmp3 <- dimnames(get(getOption("contrasts")[[if(is.ordered(tvcv[,i]))2 else 1]])(levels(tvcv[,i]),contrasts =TRUE))[[2]]
					if(is.null(tmp3))tmp3 <- 1:(length(levels(tvcv[,i]))-1)
					tmp2 <- c(tmp2,paste(names[i],tmp3,sep=""))}}
#					tmp2 <- c(tmp2,paste(names[i],levels(as.factor(tvcv[,i]))[-1],sep=""))}}
			tvcv <- tmp
			names <- tmp2
			ncols <- length(names)
			rm(tmp,tmp2,mt)}
		# create colnames if necessary
		if(is.null(colnames(tvcov[[1]]))){
			if(is.null(names))names <- paste(deparse(substitute(tvcov)))}
		else if(length(names)!=ncols)names <- colnames(tvcov[[1]])
		if((length(names)==1&&ncols>1))
			names <- paste(names,1:ncols,sep="")
		if(length(names)!=ncols)
			stop(paste(ncols,"variable names required"))
		colnames(tvcv) <- names}}
else if(!is.data.frame(tvcov))
	stop("Intra-unit (time-varying) covariates must be a matrix, dataframe, or list")
if(!is.null(interaction)){
#
# if necessary, create interactions
#
	name <- colnames(tvcv)
	units <- tvcov$units
	if(is.character(interaction)){
	# if names supplied, find corresponding columns
		mat <- match(interaction,name)
		if(any(is.na(mat)))
			stop(paste("Intra-unit (time-varying) covariate(s)",ccov[is.na(mat)],"not found"))
		interaction <- mat}
	if(is.vector(interaction,mode="numeric")){
		if(length(interaction)>length(name))
			stop("too many interactions")
		if(!is.data.frame(tvcv))colnames(tvcv) <- NULL
		if(!is.null(ccov)){
		# if interactions with time-constant covariates
			if(inherits(ccov,"tccov")){
			  ## bruce swihart edit:
			  ## switch the next two lines to avoid
			  ## R CMD Check error
			  ## in examples of tvctomat.Rd
			  units2 <- ccov$units
			  ccov <- ccov$ccov
			  if(!is.null(units)&&is.null(units2))
					units2 <- rep("NA",dim(ccov)[2])
				if(is.null(units)&&!is.null(units2))
					units <- rep("NA",dim(tvcov$tvcov)[2])}
			if(!is.matrix(ccov)&&!is.data.frame(ccov)){
			# if a vector, transform to matrix
				tmp <- paste(deparse(substitute(ccov)))
				ccov <- matrix(ccov)
				colnames(ccov) <- tmp}
			if(dim(ccov)[1]!=length(nbs))
				stop("ccov does not have one observation per individual")
			# find desired covariates
			if(is.null(names))names <- colnames(ccov)
			mat <- match(names,colnames(ccov))
			if(any(is.na(mat)))stop(paste("covariates",names[is.na(mat)],"not found"))
			oldtvcov <- tvcv
			if(!is.data.frame(oldtvcov))colnames(oldtvcov) <- name
			if(!is.data.frame(oldtvcov)&&!is.data.frame(ccov)){
			# calculate interactions for ordinary matrices
				for(i in 1:length(interaction))for(j in 1:length(mat)){
					oldtvcov <- cbind(oldtvcov,tvcv[,interaction[i]]*rep(ccov[,mat[j]],nbs))
					name <-  c(name,paste(name[interaction[i]],".",names[j],sep=""))
					if(!is.null(units))units <- c(units,paste(units[interaction[i]],".",units2[j],sep=""))}}
			else {
			# calculate interactions when expansion required
				for(i in 1:length(interaction))
					for(j in 1:length(mat)){
					mt <- terms(~tvcv[,interaction[i]]:rep(ccov[,mat[j]],nbs))
					tmp <- model.matrix(mt,model.frame(mt,na.action=NULL))[,-1,drop=FALSE]
			                if(!is.vector(tvcv[,interaction[i]],mode="numeric")){
			                	if(!is.vector(ccov[,mat[j]],mode="numeric")){
			                		nam <- NULL
			                		tmp2 <- paste(name[interaction[i]],levels(tvcv[,interaction[i]])[-1],".",sep="")
			                		for(k in 1:length(levels(ccov[,mat[j]])[-1]))
			                		nam <- c(nam,paste(tmp2,names[j],levels(ccov[,mat[j]])[-1][k],sep=""))
			                		tmp <- tmp[,-c(1:length(levels(tvcv[,interaction[i]])),seq(1,length(levels(tvcv[,interaction[i]]))*length(levels(ccov[,mat[j]])),by=length(levels(tvcv[,interaction[i]])))),drop=FALSE]}
			                	else {
			                		nam <- paste(paste(name[interaction[i]],levels(tvcv[,interaction[i]])[-1],sep=""),".",names[j],sep="")
			                		tmp <- tmp[,-1,drop=FALSE]}
							if(!is.null(units))units <- c(units,rep(units2[j],length(levels(tvcv[,interaction[i]])[-1])))}
			                else {
			                	if(!is.vector(ccov[,mat[j]],mode="numeric")){
			                		nam <- paste(name[interaction[i]],".",paste(names[j],levels(ccov[,mat[j]])[-1],sep=""),sep="")
			                		tmp <- tmp[,-1,drop=FALSE]
							if(!is.null(units))units <- c(units,rep(units[interaction[i]],length(levels(ccov[,mat[j]])[-1])))}
			                	else {
							nam <- paste(name[interaction[i]],".",names[j],sep="")
							if(!is.null(units))units <- c(units,paste(units[interaction[i]],".",units2[j],sep=""))}}
			                colnames(tmp) <- nam
					name <- c(name,nam)
					oldtvcov <- cbind(oldtvcov,tmp)}}
			if(!is.data.frame(oldtvcov))colnames(oldtvcov) <- name
			oldtvcov <- list(tvcov=oldtvcov,nobs=nbs,units=units)}
		else if(length(interaction)==2){
		# one pair of interactions
			if(is.data.frame(tvcv)){
			# expand dataframe
				mt <- terms(~tvcv[,interaction[1]]:tvcv[,interaction[2]])
				tmp <- model.matrix(mt,model.frame(mt,na.action=NULL))[,-1,drop=FALSE]
				if(!is.vector(tvcv[,interaction[1]],mode="numeric")){
					if(!is.vector(tvcv[,interaction[2]],mode="numeric")){
						names <- NULL
						tmp2 <- paste(name[interaction[1]],levels(tvcv[,interaction[1]])[-1],".",sep="")
						for(i in 1:length(levels(tvcv[,interaction[2]])[-1]))
						names <- c(names,paste(tmp2,name[interaction[2]],levels(tvcv[,interaction[2]])[-1][i],sep=""))
						tmp <- tmp[,-c(1:length(levels(tvcv[,interaction[1]])),seq(1,length(levels(tvcv[,interaction[1]]))*length(levels(tvcv[,interaction[2]])),by=length(levels(tvcv[,interaction[1]])))),drop=FALSE]}
					else {
						names <- paste(paste(name[interaction[1]],levels(tvcv[,interaction[1]])[-1],sep=""),".",name[interaction[2]],sep="")
						tmp <- tmp[,-1,drop=FALSE]
						if(!is.null(units))units <- c(units,paste(rep(units[interaction[1]],length(levels(tvcv[,interaction[1]])[-1])),".",units[interaction[2]],sep=""))}}
				else {
					if(!is.vector(tvcv[,interaction[2]],mode="numeric")){
						names <- paste(name[interaction[1]],".",paste(name[interaction[2]],levels(tvcv[,interaction[2]])[-1],sep=""),sep="")
						tmp <- tmp[,-1,drop=FALSE]
						if(!is.null(units))units <- c(units,paste(units[interaction[1]],".",rep(units[interaction[2]],length(levels(tvcv[,interaction[2]])[-1])),sep=""))}
					else {
						names <- paste(name[interaction[1]],".",name[interaction[2]],sep="")
						if(!is.null(units))units <- c(units,paste(units[interaction[1]],".",units[interaction[2]],sep=""))}}
				colnames(tmp) <- names
				oldtvcov <- list(tvcov=cbind(tvcv,tmp),nobs=nbs,units=units)}
			else {
				units <- if(is.null(tvcov$units))NULL else c(tvcov$units,paste(tvcov$units[interaction[1]],".",tvcov$units[interaction[2]],sep=""))
				oldtvcov <- list(tvcov=cbind(tvcv,tvcv[,interaction[1]]*tvcv[,interaction[2]]),nobs=nbs,units=units)}
			if(!is.data.frame(oldtvcov$tvcov))
				colnames(oldtvcov$tvcov) <- c(name,paste(name[interaction[1]],".",name[interaction[2]],sep=""))}
		else stop("interaction must be a vector containing column numbers or variable names")
		class(oldtvcov) <- "tvcov"}}
else if(!is.null(oldtvcov)){
#
# check for variable descriptions
#
if(!is.null(description)){
	if(!is.list(description))stop("description must be a list")
	if(!all(tmp <- names(description)%in%colnames(tvcv)))
		stop(paste("variable(s)",names(description)[!tmp],"not found"))
	for(i in description)if(!is.character(i))
		stop("description list must contain character vectors")}
#
# if old tvcov, combine with new one
#
	if(!inherits(oldtvcov,"tvcov"))
		stop("oldtvcov must have class, tvcov")
	oldtvcov$description <- c(oldtvcov$description,description)
	if(!is.null(oldtvcov$units)||!is.null(units)){
		if(is.null(oldtvcov$units))
			oldtvcov$units <- rep(NA,dim(oldtvcov$tvcov)[2])
		if(is.null(units))units <- rep(NA,dim(tvcv)[2])
		oldtvcov$units <- c(oldtvcov$units,units)}
	if((dim(oldtvcov$tvcov)[1]==dim(tvcv)[1])&&all(oldtvcov$nobs==nbs)){
		if(dataframe)oldtvcov$tvcov <- data.frame(oldtvcov$tvcov,tvcv)
		else oldtvcov$tvcov <- cbind(oldtvcov$tvcov,tvcv)}
	else stop("old and new covariates do not have the same numbers of observations")}
else {
	# check units
	if(!is.null(units)){
		if(!is.character(units))
			stop("units must be a character vector")
		if(length(units)!=dim(tvcv)[2])
			stop("units must be supplied for all covariates")}
	if(dataframe)tvcv <- as.data.frame(tvcv)
	oldtvcov <- list(tvcov=tvcv,nobs=nbs,units=units,description=description)
	class(oldtvcov) <- "tvcov"}
if(!is.null(oldtvcov$units))names(oldtvcov$units) <- colnames(oldtvcov$tvcov)
#
# if no factor variables present, return a matrix anyway
#
if(is.data.frame(oldtvcov$tvcov)){
	fac <- FALSE
	for(i in 1:dim(oldtvcov$tvcov)[2])if(!is.vector(oldtvcov$tvcov[,i],mode="numeric")){
		fac <- TRUE
		break}
	if(!fac)oldtvcov$tvcov <- as.matrix(oldtvcov$tvcov)}
oldtvcov}

### functions to create a repeated object
###
### method to create repeated object removing NAs
###
rmna <- function(response, ccov=NULL, tvcov=NULL){
#
# if necessary, convert response
#
if(!inherits(response,"response"))response <- restovec(response)
#
# if necessary, expand nobs
#
if(length(response$nobs)==1&&response$nobs==1)
	response$nobs <- rep(1,length(response$y))
#
# if necessary, convert ccov
#
if(!is.null(ccov)){
	if(!inherits(ccov,"tccov"))ccov <- tcctomat(ccov)
	if(length(response$nobs)!=dim(ccov$ccov)[1])
		stop("Numbers of individuals for response and for inter-unit (time-constant) covariates do not agree.")}
#
# if necessary, convert tvcov
#
if(!is.null(tvcov)){
	if(!inherits(tvcov,"tvcov"))tvcov <- tvctomat(tvcov)
	if(length(response$nobs)!=length(tvcov$nobs)||
		any(response$nobs!=tvcov$nobs))
		stop("Numbers of observations for response and intra-unit (time-varying) covariates do not agree.")}
#
# create NA indicators
#
rna <- rep(TRUE,dim(response$y)[1])
for(i in 1:dim(response$y)[2])rna <- rna&!is.na(response$y[,i])
if(!is.null(response$times))rna <- rna&!is.na(response$times)
if(!is.null(response$nest))rna <- rna&!is.na(response$nest)
if(!is.null(response$wt))rna <- rna&!is.na(response$wt)
if(!is.null(response$coordinates))
	rna <- rna&!is.na(response$coordinates[,1])&!is.na(response$coordinates[,2])
if(!is.null(response$n)){
	for(i in 1:dim(response$y)[2])if(any(!is.na(response$n[,i])))
		rna <- rna&!is.na(response$n[,i])}
for(i in 1:length(response$nobs))
	if(!is.null(ccov)&&any(is.na(ccov$ccov[i,])))
		rna[covind(response)==i] <- FALSE
if(!is.null(tvcov))
	for(i in 1:dim(tvcov$tvcov)[2])rna <- rna&!is.na(tvcov$tvcov[,i])
#
# remove NAs
#
if(any(!rna)){
	# remove NAs from variables associated with response
        response$y <- response$y[rna,,drop=FALSE]
        if(!is.null(response$times))response$times <- response$times[rna]
        if(!is.null(response$nest))response$nest <- response$nest[rna]
        if(!is.null(response$wt))response$wt <- response$wt[rna]
        if(!is.null(response$coordinates))
        	response$coordinates <- response$coordinates[rna,]
        if(!is.null(response$n))response$n <- response$n[rna,,drop=FALSE]
        if(!is.null(response$censor)){
        	response$censor <- response$censor[rna,,drop=FALSE]
        	if(all(response$censor==1))response$censor <- NULL}
        if(!is.null(response$delta)&&length(response$delta)>1)
        	response$delta <- response$delta[rna,,drop=FALSE]
        if(!is.null(tvcov))tvcov$tvcov <- tvcov$tvcov[rna,,drop=FALSE]
        # correct nobs
        tmp <- NULL
        j <- c(0,cumsum(response$nobs))
        for(i in 1:length(response$nobs)){
        	tmp <- c(tmp,sum(rna[(j[i]+1):j[i+1]]))
        	if(tmp[i]==0)
        		warning(paste("Individual",i,"has no observations"))}
        response$nobs <- tmp[tmp>0]
        # remove NAs from ccov
        if(!is.null(ccov)){
        	ccov$ccov <- ccov$ccov[tmp>0,,drop=FALSE]
        	for(i in 1: dim(ccov$ccov)[2])
        		if(length(unique(ccov$ccov[,i]))==1)
        			warning(paste("covariate",colnames(ccov$ccov)[i],"has only one value\n"))}
        # remove NAs from tvcov
	if(!is.null(tvcov)){
		tvcov$nobs <- response$nobs
		for(i in 1: dim(tvcov$tvcov)[2])
			if(length(unique(tvcov$tvcov[,i]))==1)
			warning(paste("covariate",colnames(tvcov$tvcov)[i],"has only one value\n"))}}
#
# if independent observations, reset nobs
#
if(all(response$nobs==1))response$nobs <- 1
z <- list(response=response,tvcov=tvcov,ccov=ccov)
class(z) <- "repeated"
z}

### method to create repeated object leaving NAs
###
lvna <- function(response, ccov=NULL, tvcov=NULL){
#
# if necessary, convert response
#
if(!inherits(response,"response"))response <- restovec(response)
#
# if necessary, expand nobs
#
if(length(response$nobs)==1&&response$nobs==1)
	response$nobs <- rep(1,length(response$y))
#
# if necessary, convert ccov
#
if(!is.null(ccov)){
	if(!inherits(ccov,"tccov"))ccov <- tcctomat(ccov)
	if(length(response$nobs)!=dim(ccov$ccov)[1])
		stop("Numbers of individuals for response and for intra-unit (time-constant) covariates do not agree.")}
#
# if necessary, convert tvcov
#
if(!is.null(tvcov)){
	if(!inherits(tvcov,"tvcov"))tvcov <- tvctomat(tvcov)
	if(any(response$nobs!=tvcov$nobs))
		stop("Numbers of observations for response and intra-unit (time-varying) covariates do not agree.")}
#
# create NA indicators
#
rna <- rep(TRUE,dim(response$y)[1])
for(i in 1:dim(response$y)[2])rna <- rna&!is.na(response$y[,i])
if(!is.null(response$times))rna <- rna&!is.na(response$times)
if(!is.null(response$nest))rna <- rna&!is.na(response$nest)
if(!is.null(response$coordinates))
	rna <- rna&!is.na(response$coordinates[,1])&!is.na(response$coordinates[,2])
if(!is.null(response$n)){
	for(i in 1:dim(response$y)[2])if(any(!is.na(response$n[,i])))
		rna <- rna&!is.na(response$n[,i])}
for(i in 1:length(response$nobs))
	if(!is.null(ccov)&&any(is.na(ccov$ccov[i,])))rna[covind(response)==i] <- FALSE
if(!is.null(tvcov))
	for(i in 1:dim(tvcov$tvcov)[2])rna <- rna&!is.na(tvcov$tvcov[,i])
#
# if independent observations, reset nobs
#
if(all(response$nobs==1))response$nobs <- 1
z <- list(response=response,tvcov=tvcov,ccov=ccov,
	NAs=if(any(!rna))!rna else NULL)
class(z) <- "repeated"
z}

### method to create repeated object from a dataframe
###
dftorep <- function(dataframe, response, id=NULL, times=NULL, censor=NULL,
	totals=NULL, weights=NULL, nest=NULL, delta=NULL,
	coordinates=NULL, type=NULL, ccov=NULL, tvcov=NULL, na.rm=TRUE){
if(missing(dataframe)||!is.data.frame(dataframe))
	stop("a dataframe must be supplied")
if(missing(response)||!is.character(response))
	stop("name(s) of response variables must be supplied")
#
# find response information and construct object
#
cn <- colnames(dataframe)
nc <- match(response,cn)
if(any(is.na(nc)))stop(paste("response",response[is.na(nc)],"not found"))
tot <- matrix(NA,ncol=nc,nrow=dim(dataframe)[1])
for(i in nc)if(is.matrix(dataframe[[i]])&&dim(dataframe[[i]])[2]==2){
		tot[,i] <- dataframe[[i]][,1]+dataframe[[i]][,2]
		dataframe[[i]] <- dataframe[[i]][,1]}
if(all(is.na(tot)))tot <- NULL
if(!is.numeric(z <- as.matrix(dataframe[,nc,drop=FALSE])))
	stop("response must be numeric")
z <- list(response=list(y=z,nobs=NULL,times=NULL,nest=NULL,coordinates=NULL,
	censor=NULL,n=tot,wt=NULL,delta=NULL,units=NULL,type=NULL),
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
if(!is.null(times))for(i in unique(id))
	if(!is.null(nest))for(j in unique(z$response$nest)){
		if(any(diff(z$response$times[id==i&z$response$nest==j])<0,na.rm=TRUE))
			stop(paste("negative time step for individual",i))}
	else
		if(any(diff(z$response$times[id==i])<0,na.rm=TRUE))
			stop(paste("negative time step for individual",i))
if(!is.null(nest))for(i in unique(id))
	if(any(diff(z$response$nest[id==i])!=0&
		diff(z$response$nest[id==i])!=1,na.rm=TRUE))
		stop(paste("nest for individual",i,"not consecutively numbered"))
if(!is.null(censor)){
	if(!is.character(censor)||length(censor)!=nrcol)
		stop("censor must have one name per response variable")
	nc <- match(censor,cn)
	if(any(is.na(nc)))stop("censor",censor[is.na(nc)],"not found")
	z$response$censor <- as.matrix(dataframe[,nc,drop=FALSE])
	if(!is.numeric(z$response$censor))stop("censor must be numeric")
	if(any(z$response$censor!=1&z$response$censor!=0&
		z$response$censor!=-1,na.rm=TRUE))
		stop("censor indicator can only have values, -1, 0, 1")
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

### methods for objects created by these functions

### print methods
###
print.response <- function(x, nindmax=50, ...){
  z <- x; rm(x)
nobs <- nobs(z)
nind <- length(nobs)
cn <- colnames(z$y)
tmp <- rbind(colnames(z$y),z$type)
rntmp <- c(if(length(cn)>1)"Names:" else "Name: ","Type:")
if(!is.null(z$units)){
	tmp <- rbind(tmp,z$units)
	rntmp <- c(rntmp,"Units:")}
rownames(tmp) <- rntmp
colnames(tmp) <- rep("",length(cn))
print(tmp,quote=FALSE)
if(any(nobs>1))cat("Number of individuals:                ",nind,"\n")
cat("Total number of observations:         ",dim(z$y)[1],"\n")
if(nind>1&&any(nobs>1)){
	if(all(diff(nobs)==0))
		cat("Number of observations per individual:",nobs[1],"\n")
 	else {
		if(nind>nindmax)cat("Range of observations per individual: ",
			range(nobs),"\n")
		else cat("Number of observations per individual:\n",
			nobs,"\n")}}
av <- rg <- NULL
for(i in 1:dim(z$y)[2]){
      y <- if(is.null(z$n)||all(is.na(z$n[,i])))z$y[,i] else z$y[,i]/z$n[,i]
      av <- c(av,if((z$type[i]=="nominal"&&(is.null(z$n)||all(is.na(z$n[,i]))))
	||z$type[i]=="ordinal")NA else mean(y,na.rm=TRUE))
      rg <- c(rg,range(y,na.rm=TRUE))}
names(av) <- cn
rg <- matrix(rg,nrow=2)
dimnames(rg) <- list(c("lower","upper"),cn)
if(!all(is.na(av))){
	cat("Mean response:")
	if(dim(z$y)[2]>1){
		cat("\n")
		print(av)
		cat("\n")}
	else cat("                        ",av,"\n")}
cat("Range of responses:")
if(dim(z$y)[2]>1){
	cat("\n")
	print(rg)
	cat("\n")}
else cat("                   ",rg,"\n")
if(any(is.na(z$y)))
	cat("Number of NAs:                        ",sum(is.na(z$y)),"\n")
if(!is.null(z$wt))
	cat("Number of positive weights:           ",sum(z$wt>0),"\n")
if(!is.null(z$times)){
	cat("Mean time:                            ",mean(z$times,na.rm=TRUE),"\n")
	cat("Range of times:                       ",range(z$times,na.rm=TRUE),"\n")
	if(nind>1&&is.null(z$nest)){
		mn <- if(any(z$times<0,na.rm=TRUE))z$times[cumsum(c(1,nobs[1:(nind-1)]))]
			else 0
		cat("Mean total time:                      ",mean(z$times[cumsum(nobs)]-mn,na.rm=TRUE),"\n")
		cat("Range of total times:                 ",range(z$times[cumsum(nobs)]-mn,na.rm=TRUE),"\n")}}
if(!is.null(z$nest))
	cat("Maximum number of clusters:           ",max(z$nest),"\n")
if(!is.null(z$censor)){
	cens0 <- cens1 <- n0 <- n1 <- NULL
	for(i in 1:dim(z$y)[2]){
		if(sum(z$censor[,i]==0,na.rm=TRUE)>0){
			cens0 <- c(cens0,sum(z$censor[,i]==0,na.rm=TRUE))
			n0 <- c(n0,cn[i])}
		if(sum(z$censor[,i]==-1,na.rm=TRUE)>0){
			cens1 <- c(cens1,sum(z$censor[,i]==-1,na.rm=TRUE))
			n1 <- c(n1,cn[i])}}
	names(cens0) <- n0
	names(cens1) <- n1
	if(!is.null(cens0)){
		cat("Number of right-censored observations:")
		if(length(cens0)>1){
			cat("\n")
			print(cens0)
			cat("\n")}
		else cat(" ",cens0," (",n0,")\n",sep="")}
	if(!is.null(cens1)){
		cat("Number of left-censored observations:")
		if(length(cens1)>1){
			cat("\n")
			print(cens1)
			cat("\n")}
		else cat("  ",cens1," (",n1,")\n",sep="")}}
if(!is.null(z$delta)&&length(z$delta)==1)
	cat("Unit of measurement:                  ",z$delta,"\n")
if(!is.null(z$description))for(i in 1:length(z$description)){
	cat(names(z$description)[i],": ",sep="")
	cat(z$description[[i]],"\n")}}

print.tccov <- function(x, ...){
  z <- x; rm(x)
if(is.function(z))print.default(unclass(z))
else {
	tmp <- matrix(colnames(z$ccov),nrow=1)
	rn <- "Names:"
	if(!is.null(z$units)){
		tmp <- rbind(tmp,z$units)
		rn <- c(rn,"Units:")}
	dimnames(tmp) <- list(rn,rep("",dim(z$ccov)[2]))
	print(tmp,quote=FALSE)
	cat("Number of individuals:             ",dim(z$ccov)[1],"\n")
	if(!is.null(z$description))for(i in 1:length(z$description)){
		cat(names(z$description)[i],": ",sep="")
		cat(z$description[[i]],"\n")}}}

print.tvcov <- function(x, nindmax=50, ...){
  z <- x; rm(x)
if(is.function(z))print.default(unclass(z))
else {
	tmp <- matrix(colnames(z$tvcov),nrow=1)
	rn <- "Names:"
	if(!is.null(z$units)){
		tmp <- rbind(tmp,z$units)
		rn <- c(rn,"Units:")}
	dimnames(tmp) <- list(rn,rep("",dim(z$tvcov)[2]))
	print(tmp,quote=FALSE)
	cat("Number of individuals:               ",length(nobs(z)),"\n")
	cat("Number of observations:              ",sum(nobs(z)),"\n")
	if(length(nobs(z))>nindmax)cat("Range of observations per individual:",
		range(nobs(z)),"\n")
	else cat("Number of observations per individual:\n",
		nobs(z),"\n")
	if(!is.null(z$description))for(i in 1:length(z$description)){
		cat(names(z$description)[i],": ",sep="")
		cat(z$description[[i]],"\n")}}}

print.repeated <- function(x, nindmax=50, ...){
  z <- x; rm(x)
if(is.function(z))print.default(unclass(z))
else {
	cat("\nResponse variable:\n")
	print.response(z$response,nindmax=nindmax)
	if(!is.null(z$ccov)){
		cat("\nInter-unit (time-constant) covariates:\n")
		print.tccov(z$ccov)}
	if(!is.null(z$tvcov)){
		cat("\nIntra-unit (time-varying) covariates:\n")
		print.tvcov(z$tvcov,nindmax=nindmax)}}}

### plot methods
###
plot.response <- function(x, name=NULL, nind=NULL, nest=1, ccov=NULL,
	add=FALSE, lty=NULL, pch=NULL, main=NULL, ylim=NULL, xlim=NULL,
	xlab=NULL, ylab=NULL, ...){
if(is.null(name)){
	if(dim(x$y)[2]>1)stop("please specify which variable to plot")
	name <- colnames(x$y)}
else if(length(name)>1)stop("only one variable can be plotted")
col <- match(name,colnames(x$y))
if(is.na(col))stop(paste(name,"not found"))
if(x$type[col]=="ordinal"){
	# special case: ordinal response
	#return(plot.ordinal(z=x,ccov=ccov,main=main,xlab=xlab,ylab=ylab,
		#xlim=xlim,ylim=ylim,lty=lty,add=add,...))
  stop(paste("email Bruce to email Jim for plot.ordinal() implementation"))
  }
if(is.null(ylab))ylab <- name
if(is.null(x$times)){
#
# when no times, set up for index plot
#
	if(all(nobs(x)==1)){
		x$times <- 1:length(nobs(x))
		# set to value different from 1 so not time series
		x$nobs <- 5}
	else x$times <- sequence(nobs(x))
	if(is.null(xlab)) xlab <- "Index number"}
else if(is.null(xlab)) xlab <- "Time"
if(is.null(ylim))ylim <- range(x$y[,col],na.rm=TRUE)
if(is.null(xlim))xlim <- range(x$times,na.rm=TRUE)
tnest <- if(!is.null(x$nest)) x$nest else 1
#
# initialize
#
nm <- covind(x)
j <- 1
lt <- 0
#
# if no individuals chosen, plot them all
#
if(is.null(nind))nind <- 1:length(nobs(x))
#
# if binomial, plot proportions
#
y <- if(is.null(x$n)||all(is.na(x$n[,col])))x$y[,col] else x$y[,col]/x$n[,col]
if(!is.null(lty)){
#
# set up line types
#
	if(length(lty)==1)lty <- rep(lty,length(nind))
	else if(length(lty)!=length(nind))
		stop("lty must have one value for each item in nind")}
if(!is.null(pch)){
#
# set up symbol types
#
	if(length(pch)==1)pch <- rep(pch,length(nind))
	else if(length(pch)!=length(nind))
		stop("pch must have one value for each item in nind")}
for(i in 1:length(nobs(x)))if(any(i==nind))for(k in nest){
	lt <- if(is.null(lty))lt%%6+1 else lty[j]
	if(!add&&j==1)plot(x$times[nm==i&k==tnest],y[nm==i&k==tnest],lty=lt,
		type="l",ylim=ylim,xlim=xlim,main=main,ylab=ylab,xlab=xlab,...)
	else lines(x$times[nm==i&k==tnest],y[nm==i&k==tnest],lty=lt)
	if(!is.null(pch))points(x$times[nm==i&k==tnest],
		y[nm==i&k==tnest],pch=pch[j])
	j <- j+1}}

plot.repeated <- function(x, name=NULL, nind=NULL, nest=1, ccov=NULL,
	add=FALSE, lty=NULL, pch=NULL, main=NULL, ylim=NULL, xlim=NULL,
	xlab=NULL, ylab=NULL, ...){
if(is.null(name)){
	if(dim(x$response$y)[2]>1)stop("please specify which variable to plot")
	name <- colnames(x$response$y)}
else if(length(name)>1)stop("only one variable name can be supplied")
col <- match(name,colnames(x$response$y))
if(is.na(col)){
	col <- match(name,colnames(x$tvcov$tvcov))
	if(is.na(col))stop(paste(name,"not found"))
	variable <- "tvc"}
else variable <- "response"
#
# check individuals to plot
# if no individuals chosen, plot them all
#
uncov <- NULL
if(!is.null(nind)&&!is.null(ccov))
	stop("only one of nind and ccov can be specified")
if(is.null(nind)&&is.null(ccov))nind <- 1:length(nobs(x))
else if(!is.null(ccov)&&x$response$type[col]!="ordinal"){
	if(is.numeric(ccov)){
		if(length(ccov)!=dim(x$ccov$ccov)[2])
			stop("a covariate value must be given for each covariate")
	        tccov <- x$ccov$ccov
	        if(is.data.frame(tccov))for(i in 1:length(ccov)){
	        	if(is.factor(tccov[[i]]))
	        		tccov[[i]] <- as.numeric(tccov[[i]])}
	        tccov <- as.matrix(tccov)
		nind <- NULL
		for(i in 1:length(nobs(x)))
			if(all(ccov==tccov[i,]))nind <- c(nind,i)}
	else if(is.character(ccov)){
		if(is.null(x$ccov$ccov))stop("no covariates found")
		if(length(ccov)>1)
			stop("only one variable name can be given in ccov")
		col2 <- match(ccov,colnames(x$ccov$ccov))
		if(is.na(col2))stop(paste("covariate",ccov,"not found"))
		uncov <- unique(x$ccov$ccov[,col2])
		if(is.factor(uncov))uncov <- as.character(uncov)
		if(length(uncov)>6)
			stop(paste(ccov,"has too many distinct values to plot"))
		if(length(uncov)==1)uncov <- NULL}
	else stop("ccov must be a vector of covariate values or a covariate name")}
if(variable=="response"){
	# special case: ordinal response
	if(x$response$type[col]=="ordinal")
		#plot.ordinal(z=x,ccov=ccov,main=main,xlab=xlab,ylab=ylab,
		#	xlim=xlim,ylim=ylim,lty=lty,add=add,...)
	  stop(paste("email Bruce to email Jim for plot.ordinal() implementation"))
	else {
		if(is.null(uncov))plot.response(x$response,name=name,nind=nind,
			nest=nest,add=add,lty=lty,pch=pch,main=main,
			ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,...)
		else {
			mfrow <- if(length(uncov)==2)c(2,1)
				else if(length(uncov)<5)c(2,2) else c(3,2)
			oldpar <- par(mfrow=mfrow,no.readonly=TRUE)
			for(i in uncov){
				nind <- NULL
				main <- paste(ccov,"=",i)
				for(j in 1:length(nobs(x)))
					if(i==x$ccov$ccov[j,col2])
					nind <- c(nind,j)
				plot.response(x$response,name=name,nind=nind,
				nest=nest,add=add,lty=lty,pch=pch,main=main,
				ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,...)}
			par(oldpar)}}}
else if(variable=="tvc"){
	# set up covariate as a "response"
	if(is.null(ylab))ylab <- colnames(x$tvcov$tvcov)[col]
	zz <- list(times=x$response$times,y=x$tvcov$tvcov[,col,drop=FALSE],
		nobs=x$tvcov$nobs,n=NULL,type="unknown")
	class(zz) <- "response"
	if(is.null(uncov))plot.response(zz,name=name,nind=nind,nest=nest,
		add=add,lty=lty,pch=pch,main=main,ylim=ylim,xlim=xlim,
		xlab=xlab,ylab=ylab,...)
	else {
		mfrow <- if(length(uncov)==2)c(2,1)
			else if(length(uncov)<5)c(2,2) else c(3,2)
		oldpar <- par(mfrow=mfrow,no.readonly=TRUE)
		for(i in uncov){
			nind <- NULL
			main <- paste(ccov,"=",i)
			for(j in 1:length(nobs(x)))
				if(i==x$ccov$ccov[j,col2])
				nind <- c(nind,j)
			plot.response(zz,name=name,nind=nind,nest=nest,add=add,
			lty=lty,pch=pch,main=main,ylim=ylim,xlim=xlim,
			xlab=xlab,ylab=ylab,...)}
		par(oldpar)}}}

### methods to find the response variable
###
response <- function(z, ...) UseMethod("response")

response.response <- function(z, nind=NULL, names=NULL, ...){
if(is.null(nind))nind <- 1:dim(z$y)[1]
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nind <- !is.na(match(covind(z),nind))
if(all(!nind))stop("No such individuals")
if(is.null(names))col <- 1:dim(z$y)[2]
else {
	col <- match(names,colnames(z$y))
	if(any(is.na(col)))stop(paste(names[is.na(col)],"not found"))}
cn <- tmp <- NULL
for(i in col){
	if(!is.null(z$n)&&!any(is.na(z$n[,i]))){
	# binomial response
		tmp <- cbind(tmp,z$y[,i],z$n[,i]-z$y[,i])
		cn <- c(cn,colnames(z$y)[i],paste("n-",colnames(z$y)[i],sep=""))}
	else if(!is.null(z$censor)&&!any(is.na(z$censor[,i]))){
	# censored response with indicator
		tmp <- cbind(tmp,z$y[,i],z$censor[,i])
		cn <- c(cn,colnames(z$y)[i],paste(colnames(z$y)[i],".cens",sep=""))}
	else {
		tmp <- cbind(tmp,z$y[,i])
		cn <- c(cn,colnames(z$y)[i])}}
colnames(tmp) <- cn
tmp[nind,]}

response.repeated <- function(z, nind=NULL, names=NULL, ...){
if(is.null(nind))nind <- 1:dim(z$response$y)[1]
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nind <- !is.na(match(covind(z),nind))
if(all(!nind))stop("No such individuals")
if(is.null(names))col <- 1:dim(z$response$y)[2]
else {
	col <- match(names,colnames(z$response$y))
	if(any(is.na(col)))stop(paste(names[is.na(col)],"not found"))}
cn <- tmp <- NULL
for(i in col){
	if(!is.null(z$response$n)&&!any(is.na(z$response$n[,i]))){
	# binomial response
		tmp <- cbind(tmp,z$response$y[,i],z$response$n[,i]-z$response$y[,i])
		cn <- c(cn,colnames(z$response$y)[i],paste("n-",colnames(z$response$y)[i],sep=""))}
	else if(!is.null(z$response$censor)&&!any(is.na(z$response$censor[,i]))){
	# censored response with indicator
		tmp <- cbind(tmp,z$response$y[,i],z$response$censor[,i])
		cn <- c(cn,colnames(z$response$y)[i],paste(colnames(z$response$y)[i],".cens",sep=""))}
	else {
		tmp <- cbind(tmp,z$response$y[,i])
		cn <- c(cn,colnames(z$response$y)[i])}}
colnames(tmp) <- cn
tmp[nind,]}

### methods for indexing of time-constant covariates for individuals
###
covind <- function(z, ...) UseMethod("covind")

covind.default <- function(z, ...) rep(1:length(nobs(z)),nobs(z))

### methods to find numbers of observations/individual
###
nobs <- function(z, ...) UseMethod("nobs")

nobs.default <- function(z, ...) {
if(is.null(z$response)||is.null(z$response$nobs))return(NULL)
if(length(z$response$nobs)>1||z$response$nobs>1)z$response$nobs
else rep(1,dim(z$response$y)[1])}

nobs.response <- function(z, ...) {
if(length(z$nobs)>1||z$nobs>1)z$nobs
else rep(1,length(z$y))}

nobs.tvcov <- function(z, ...) z$nobs

nobs.data.frame <- function(z, ...) rep(1,dim(z)[1])

### methods to find times
###
times <- function(z, ...) UseMethod("times")

times.default <- function(z, nind=NULL, ...){
if(is.null(nind)||is.null(z$response$times))return(z$response$times)
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nind <- !is.na(match(covind(z),nind))
if(all(!nind))stop("No such individuals")
z$response$times[nind]}

times.response <- function(z, nind=NULL, ...){
if(is.null(nind)||is.null(z$times))return(z$times)
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nind <- !is.na(match(covind(z),nind))
if(all(!nind))stop("No such individuals")
z$times[nind]}

### methods to find unit of measurement/Jacobian
###
delta <- function(z, ...) UseMethod("delta")

delta.response <- function(z, nind=NULL, names=NULL, ...){
#
# find individuals
#
if(is.null(z$delta))return(NULL)
if(is.null(nind))nind <- 1:dim(z$y)[1]
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nind <- !is.na(match(covind(z),nind))
if(all(!nind))stop("No such individuals")
#
# find variables
#
col <- if(is.null(names))1:dim(z$y)[2] else match(names,colnames(z$y))
if(any(is.na(col)))stop(paste(names[is.na(col)],"not found"))
z$delta[nind,col]}

delta.repeated <- function(z, nind=NULL, names=NULL, ...){
#
# find individuals
#
if(is.null(z$response$delta))return(NULL)
if(is.null(nind))nind <- 1:dim(z$response$y)[1]
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nind <- !is.na(match(covind(z),nind))
if(all(!nind))stop("No such individuals")
#
# find variables
#
col <- if(is.null(names))1:dim(z$response$y)[2]
	else match(names,colnames(z$response$y))
print(col)
if(any(is.na(col)))stop(paste(names[is.na(col)],"not found"))
z$response$delta[nind,col]}

### methods to find weights
###
#weights <- function(object, ...) UseMethod("weights")

weights.response <- function(object, nind=NULL, ...){
if(is.null(nind)||is.null(object$wt))return(object$wt)
else if(length(nind)>length(nobs(object))||any(nind>length(nobs(object))))
	stop("Individual not found")
else nind <- !is.na(match(covind(object),nind))
if(all(!nind))stop("No such individuals")
object$wt[nind]}

weights.repeated <- function(object, nind=NULL, ...){
if(is.null(nind)||is.null(object$response$wt))return(object$response$wt)
else if(length(nind)>length(nobs(object))||any(nind>length(nobs(object))))
	stop("Individual not found")
else nind <- !is.na(match(covind(object),nind))
if(all(!nind))stop("No such individuals")
object$response$wt[nind]}

### methods to find nesting indicators
###
nesting <- function(z, ...) UseMethod("nesting")

nesting.response <- function(z, nind=NULL, ...){
if(is.null(nind))nind <- 1:length(z$y)
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nind <- !is.na(match(covind(z),nind))
if(all(!nind))stop("No such individuals")
#
# if nest is NULL, return individual index otherwise both
#
if(length(nobs(z))==1||all(nobs(z)==1))return(NULL)
else if(is.null(z$nest))return(covind(z)[nind])
else {
	z <- cbind(covind(z),z$nest)
	colnames(z) <- c("Individual","Cluster")
	return(z)[nind,,drop=FALSE]}}

nesting.repeated <- function(z, nind=NULL, ...){
if(is.null(nind))nind <- 1:dim(z$response$y)[1]
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nind <- !is.na(match(covind(z),nind))
if(all(!nind))stop("No such individuals")
#
# if nest is NULL, return individual index otherwise both
#
if(length(nobs(z))==1||all(nobs(z)==1))return(NULL)
else if(is.null(z$response$nest))return(covind(z)[nind])
else {
	z <- cbind(covind(z),z$response$nest)
	colnames(z) <- c("Individual","Cluster")
	return(z[nind,,drop=FALSE])}}

### methods to find covariates
###
covariates <- function(z, ...) UseMethod("covariates")

covariates.tccov <- function(z, nind=NULL, names=NULL, expand=FALSE, ...){
if(is.null(nind))nind <- 1:dim(z$ccov)[1]
else if(length(nind)>dim(z$ccov)[1]||any(nind>dim(z$ccov)[1])||any(nind<1))
	stop("Individual not found")
if(is.null(names))return(z$ccov[nind,])
else {
	col <- match(names,colnames(z$ccov))
	if(any(is.na(col)))
		stop(paste("covariate(s)",names[is.na(col)],"not found"))
	return(z$ccov[nind,col])}}

covariates.tvcov <- function(z, nind=NULL, names=NULL, expand=FALSE, ...){
if(is.null(nind))nind <- 1:dim(z$tvcov)[1]
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nind <- !is.na(match(covind(z),nind))
if(all(!nind))stop("No such individuals")
if(is.null(names))return(z$tvcov[nind,])
else {
	col <- match(names,colnames(z$tvcov))
	if(any(is.na(col)))
		stop(paste("covariate(s)",names[is.na(col)],"not found"))
	return(z$tvcov[nind,col])}}

covariates.repeated <- function(z, nind=NULL, names=NULL, expand=FALSE, ...){
if(expand&&!is.null(nind))stop("can only expand for all individuals")
ind <- covind(z)
if(is.null(nind)){
	nindv <- 1:dim(z$response$y)[1]
	nind <- if(expand)ind else 1:length(nobs(z))}
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nindv <- !is.na(match(ind,nind))
if(all(!nindv))stop("No such individuals")
if(is.null(names)){
	if(is.null(z$tvcov$tvcov))return(z$ccov$ccov[nind,])
	else if(is.null(z$ccov$ccov))return(z$tvcov$tvcov[nindv,])
	else {
		if(expand)return(cbind(z$ccov$ccov[nind,],
			z$tvcov$tvcov[nindv,]))
		else return(list(ccov=z$ccov$ccov[nind,],
			tvcov=z$tvcov$tvcov[nindv,]))}}
else {
	mat1 <- match(names,colnames(z$ccov$ccov))
	mat2 <- match(names,colnames(z$tvcov$tvcov))
	if(any(is.na(mat1)&is.na(mat2)))
		stop(paste("covariate(s)",names[is.na(mat1)&is.na(mat2)],"not found"))
	if(all(is.na(mat1)))return(z$tvcov$tvcov[nindv,mat2])
	else if(all(is.na(mat2)))return(z$ccov$ccov[nind,mat1])
	else {
		if(expand)return(cbind(z$ccov$ccov[nind,mat1],
			z$tvcov$tvcov[nindv,mat2]))
		else return(list(ccov=z$ccov$ccov[nind,mat1],
			tvcov=z$tvcov$tvcov[nindv,mat2]))}}}

### methods to find names
###
names.response <- function(x, ...) colnames(x$y)

names.tccov <- function(x, ...) colnames(x$ccov)

names.tvcov <- function(x, ...) colnames(x$tvcov)

names.repeated <- function(x, ...)
	list(response=colnames(x$response$y),ccov=colnames(x$ccov$ccov),
		tvcov=colnames(x$tvcov$tvcov))

### methods to find units of measurements
###
units <- function(x, ...) UseMethod("units")

units.default <- function(x, ...) x$units

units.repeated <- function(x, ...)
	list(response=units(x$response),ccov=units(x$ccov),
		tvcov=units(x$tvcov))

### methods to find description of variables
###
description <- function(z, ...) UseMethod("description")

description.default <- function(z, ...) z$description

description.repeated <- function(z, ...)
	list(response=description(z$response),ccov=description(z$ccov),
		tvcov=description(z$tvcov))

### methods to find response type(s)
###
resptype <- function(z, ...) UseMethod("resptype")

resptype.response <- function(z, ...) z$type

resptype.repeated <- function(z, ...) z$response$type

### methods to find formula used in tccov
###
formula.tccov <- function(x, ...) x$linear

formula.repeated <- function(x, ...) x$ccov$linear

### methods to transform response, times, or covariates
###
transform.response <- function(`_data`, times=NULL, units=NULL, ...){
  z <- `_data`; rm(`_data`)
if(is.call(substitute(times)))times <- substitute(times)
tran <- as.list(substitute(list(...)))[-1]
if(!is.null(tran)){
#
# transform response
#
	if(is.null(z$delta))
		z$delta <- matrix(1,nrow=dim(z$y)[1],ncol=dim(z$y)[2])
	else if(length(delta)==1)
		z$delta <- matrix(z$delta,nrow=dim(z$y)[1],ncol=dim(z$y)[2])
	cn <- colnames(z$y)
	cn2 <- names(tran)
	col <- match(cn2,cn)
	col2 <- NULL
	if(!is.null(units)){
		if(length(units)!=length(tran))
			stop(paste(length(tran),"units required"))
		if(!is.character(units))
			stop("units must be a character vector")}
	for(i in tran){
		tmp <- NULL
		for(j in 1:length(cn)){
			if(length(grep(cn[j],as.character(i)))>0){
				tmp <- j
				break}}
		if(is.null(tmp))stop("variable to transform not found")
		col2 <- c(col2,tmp)}
	j <- 0
	for(i in tran){
		j <- j+1
		if(z$type[col2[j]]!="continuous"&&z$type[col2[j]]!="duration"&&z$type[col2[j]]!="unknown")
			stop(paste("transformations do not make sense with",z$type[col2[j]],"responses"))
		# transform
		tran <- eval(deriv(i,cn[col2[j]]),as.data.frame(z$y),NULL)
		# calculate Jacobian
		jacob <- as.vector(abs(attr(tran,"gradient")))
		if(any(is.na(jacob)&!is.na(tran)))
			stop("NAs in jacobian - invalid transformation")
		if(any(abs(jacob)==Inf,na.rm=TRUE))
			stop("infinite jacobian - invalid transformation")
		if(any(jacob<=0,na.rm=TRUE))
			stop("nonpositive value in jacobian - invalid transformation")
		# store in the response object, updating delta, units, type
		if(is.na(col[j])){
			if(any(is.na(tran)&!is.na(z$y[,col2[j]]),na.rm=TRUE))
				stop("NAs created by transformation")
			tran <- as.matrix(tran)
			colnames(tran) <- cn2[j]
			z$y <- cbind(z$y,tran)
			z$delta <- if(all(is.na(z$delta[,col2[j]])))
				cbind(z$delta,jacob)
				else cbind(z$delta,z$delta[,col2[j]]*jacob)
			colnames(z$delta) <- colnames(z$y)
			if(!is.null(z$n)){
				z$n <- cbind(z$n,rep(NA,dim(z$n)[1]))
				colnames(z$n) <- colnames(z$y)}
			if(!is.null(z$censor)){
				z$censor <- cbind(z$censor,z$censor[,col2[j]])
				colnames(z$censor) <- colnames(z$y)}
			z$type <- c(z$type,z$type[col2[j]])
			names(z$type) <- colnames(z$y)
			if(!is.null(z$units)){
				z$units <- c(z$units,if(is.null(units))NA else
					units[j])
				names(z$units) <- colnames(z$y)}}
		else {
			if(any(is.na(tran)&!is.na(z$y[,col[j]]),na.rm=TRUE))
				stop("NAs created by transformation")
			z$y[,col[j]] <- as.vector(tran)
			z$delta[,col[j]] <- if(all(is.na(z$delta[,col[j]])))
				jacob else z$delta[,col[j]]*jacob
			if(!is.null(z$units))
				z$units[col[j]] <- if(is.null(units))NA
					else units[j]}}}
if(!is.null(z$units)&&all(is.na(z$units)))z$units <- NULL
if(!is.null(times)){
#
# transform times
#
	z$times <- eval(times,z,NULL)
	for(i in 1:length(nobs(z)))if(any(diff(z$times[covind(z)==i])<0))
		stop("transformation produces negative time steps")}
z}

transform.repeated <- function(`_data`, times=NULL, ...){
    z <- `_data`; rm(`_data`)
if(is.call(substitute(times)))times <- substitute(times)
z$response <- transform.response(z$response,times,...)
z}

transform.tccov <- function(`_data`, ...){
    z <- `_data`; rm(`_data`)
isf <- is.data.frame(z$ccov)
if(!isf)z$ccov <- as.data.frame(z$ccov)
#
# find transformations
#
e <- eval(substitute(list(...)),z$ccov,NULL)
tags <- names(e)
for(i in 1:length(e))if(all(is.na(e[[i]])))
	stop(paste(tags[i],"defines an invalid tranformation\n or attempts to transform a factor variable"))
#
# find covariates to transform
#
inx <- match(tags,colnames(z$ccov))
matched <- !is.na(inx)
if(any(matched))z$ccov[inx[matched]] <- e[matched]
if(!all(matched))z$ccov <- data.frame(z$ccov,e[!matched])
if(!isf)z$ccov <- as.matrix(z$ccov)
z}

transform.tvcov <- function(`_data`, ...){
    z <- `_data`; rm(`_data`)
isf <- is.data.frame(z$tvcov)
if(!isf)z$tvcov <- as.data.frame(z$tvcov)
#
# find transformations
#
e <- eval(substitute(list(...)),z$tvcov,NULL)
tags <- names(e)
for(i in 1:length(e))if(all(is.na(e[[i]])))
	stop(paste(tags[i],"defines an invalid tranformation\n or attempts to transform a factor variable"))
#
# find covariates to transform
#
inx <- match(tags,colnames(z$tvcov))
matched <- !is.na(inx)
if(any(matched))z$tvcov[inx[matched]] <- e[matched]
if(!all(matched))z$tvcov <- data.frame(z$tvcov,e[!matched])
if(!isf)z$tvcov <- as.matrix(z$tvcov)
z}

### as.data.frame methods
###
as.data.frame <- function(x, ...) UseMethod("as.data.frame")

as.data.frame.response <- function(x,row.names=NULL,optional=FALSE, ...){
	tmp <- data.frame(x$y)
	if(!is.null(x$n))for(i in 1:dim(x$y)[2])if(any(!is.na(x$n[,i])))
		tmp[[i]] <- I(cbind(x$y[,i],x$n[,i]-x$y[,i]))
	if(!is.null(x$censor))for(i in 1:dim(x$y)[2])
		if(any(!is.na(x$censor[,i])))
			tmp[[i]] <- I(cbind(x$y[,i],x$censor[,i]))
	cn <- colnames(x$y)
	if(length(x$nobs)!=1){
		cn <- c(cn,"individuals")
		tmp <- data.frame(tmp,as.factor(rep(1:length(x$nobs),x$nobs)))}
	if(!is.null(x$nest)){
		cn <- c(cn,"nesting")
		tmp <- data.frame(tmp,as.factor(x$nest))}
	if(!is.null(x$times)){
		cn <- c(cn,"times")
		tmp <- data.frame(tmp,x$times)}
	colnames(tmp) <- cn
	data.frame(tmp)}

as.data.frame.tccov <- function(x,row.names=NULL,optional=FALSE, ...)
	as.data.frame(x$ccov)

as.data.frame.tvcov <- function(x,row.names=NULL,optional=FALSE, ...)
	as.data.frame(x$tvcov)

as.data.frame.repeated <- function(x,row.names=NULL,optional=FALSE, ...){
	tmp <- data.frame(x$response$y)
	if(!is.null(x$response$n))for(i in 1:dim(x$response$y)[2])
		if(any(!is.na(x$response$n[,i])))
			tmp[[i]] <- I(cbind(x$response$y[,i],x$response$n[,i]-x$response$y[,i]))
	if(!is.null(x$response$censor))for(i in 1:dim(x$response$y)[2])
		if(any(!is.na(x$response$censor[,i])))
			tmp[[i]] <- I(cbind(x$response$y[,i],x$response$censor[,i]))
	cn <- colnames(x$response$y)
	if(!(length(x$response$nobs)==1)){
		cn <- c(cn,"individuals")
		tmp <- data.frame(tmp,as.factor(rep(1:length(x$response$nobs),
			x$response$nobs)))}
	if(!is.null(x$response$nest)){
		cn <- c(cn,"nesting")
		tmp <- data.frame(tmp,as.factor(x$response$nest))}
	if(!is.null(x$response$times)){
		cn <- c(cn,"times")
		tmp <- data.frame(tmp,x$response$times)}
	tmp <- data.frame(tmp)
	colnames(tmp) <- cn
	if(!is.null(x$ccov$ccov))
		tmp <- data.frame(tmp,x$ccov$ccov[covind(x),,drop=FALSE])
	if(!is.null(x$tvcov$tvcov))
		tmp <- data.frame(tmp,x$tvcov$tvcov)
	tmp}

### as.matrix methods
###
as.matrix <- function(x, ...) UseMethod("as.matrix")

as.matrix.response <- function(x, ...){
	tmp <- x$y
	cn <- colnames(x$y)
	if(!is.null(x$times)){
		cn <- c(cn,"times")
		tmp <- cbind(tmp,x$times)}
	colnames(tmp) <- cn
	tmp}

as.matrix.tccov <- function(x, ...) as.matrix(x$ccov)

as.matrix.tvcov <- function(x, ...) as.matrix(x$tvcov)

as.matrix.repeated <- function(x, ...){
	tmp <- x$response$y
	cn <- colnames(x$response$y)
	if(!is.null(x$response$times)){
		cn <- c(cn,"times")
		tmp <- cbind(tmp,x$response$times)}
	colnames(tmp) <- cn
	if(!is.null(x$ccov$ccov))
		tmp <- cbind(tmp,x$ccov$ccov[covind(x),,drop=FALSE])
	if(!is.null(x$tvcov$tvcov))
		tmp <- cbind(tmp,x$tvcov$tvcov)
	tmp}
