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
#     restovec(response, times=NULL, nest=NULL, coordinates=NULL,
#	 censor=NULL, totals=NULL, weights=NULL, delta=NULL)
#     tcctomat(ccov, names=NULL, oldtccov=NULL)
#     rmna(response, tvcov=NULL, ccov=NULL)
#     response(z, ...)
#     times(z, ...)
#     nesting(z, ...)
#     covariates(z, ...)
#     parameters(z, ...)
#     covind(z, ...)
#
#  DESCRIPTION
#
#    Utility functions for converting repeated measurements data to R objects

restovec <- function(response, times=NULL, nest=NULL, coordinates=NULL,
	censor=NULL, totals=NULL, weights=NULL, delta=NULL, type=NULL){
	if(missing(response))stop("A response must be supplied")
	nind <- 0
	tnest <- nobs <- y <- NULL
	ttime <- !is.logical(times)
	if(!ttime)times <- NULL
	if(is.vector(response,mode="numeric")){
		y <- response
		nobs <- if(is.null(times)) 1 else length(response)
		if(is.vector(censor,mode="numeric")){
			if(length(censor)!=length(y)){
				if(length(censor)==1)censor <- c(rep(1,length(y)-1),censor)
				else stop("Censoring vector must be the same length as the response")}}
		else if(!is.null(censor))stop("Censor must be a scalar or vector")
		if(!is.null(coordinates)&&(!is.matrix(coordinates)||(is.matrix(coordinates)&&length(dim(coordinates))!=2&&ncol(coordinates)!=2&&nrow(coordinates)!=length(y))))
			stop("coordinates must be a 2-dimensional matrix with two columns and the same number of rows as the length of response")
		if(is.vector(totals,mode="numeric")){
			if(length(totals)!=length(y)){
				if(length(totals)==1)totals <- rep(totals,length(y))
				else stop("Totals vector must be the same length as the response")}}
		else if(!is.null(totals)) stop("Totals must be a vector")
		if(is.vector(times,mode="numeric")){
			if(length(times)!=length(y))stop("Times vector must be the same length as the response")}
		else if(!is.null(times))stop("Times must be a vector")
		if(is.vector(weights,mode="numeric")){
			if(length(weights)!=length(y))stop("Weights vector must be the same length as the response")}
		else if(!is.null(weights))stop("Weights must be a vector")
		if(is.vector(delta,mode="numeric")){
			if(length(delta)!=length(y)){
				if(length(delta)==1)delta <- rep(delta,length(y))
				else stop("Delta vector must be the same length as the response")}}
		else if(!is.null(delta))stop("Delta must be a scalar or vector")}
	else if(is.matrix(response)||is.data.frame(response)){
		if(is.data.frame(response))response <- as.matrix(response)
		if(is.data.frame(totals))totals <- as.matrix(totals)
		if(is.vector(censor,mode="numeric")){
			if(length(censor)!=nrow(response))
				stop("Censoring vector must be the same length as the number of individuals")
			else {
				tmp <- matrix(1,nrow=nrow(response),ncol=ncol(response))
				tmp[,ncol(tmp)] <- censor
				censor <- tmp}}
		if(is.matrix(censor)){
			if(nrow(censor)!=nrow(response)||ncol(censor)!=ncol(response))
				stop("Censoring matrix must have the same dimensions as response")
			else censor <- as.vector(t(censor))}
		if(is.matrix(totals)){
			if(nrow(totals)!=nrow(response)||ncol(totals)!=ncol(response))
				stop("totals matrix must have the same dimensions as response")
			else totals <- as.vector(t(totals))}
		else if(!is.null(totals)&&is.vector(totals,mode="numeric")){
		     if(length(totals)!=nrow(response))stop("totals vector must have same length as number of individuals")
		     else totals <- rep(totals,rep(ncol(response),nrow(response)))}
		if(is.matrix(weights)){
			if(nrow(weights)!=nrow(response)||ncol(weights)!=ncol(response))
				stop("weights matrix must have the same dimensions as response")
			else weights <- as.vector(t(weights))}
		else if(!is.null(weights)&&is.vector(weights,mode="numeric")){
		     if(length(weights)!=nrow(response))stop("weights vector must have same length as number of individuals")
		     else weights <- rep(weights,rep(ncol(response),nrow(response)))}
		if(is.null(times)){
			if(is.null(censor)&&ttime)times <- as.double(rep(1:ncol(response),nrow(response)))}
		else if(is.vector(times,mode="numeric")) {
			if(is.null(nest)&&any(diff(times)<0,na.rm=T))stop("Times must be increasing")
			if(length(times)!=ncol(response))stop("Number of times must equal number of response columns")
			times <- rep(times,nrow(response))}
		else if(is.matrix(times)){
			if(ncol(times)!=ncol(response)|nrow(times)!=nrow(response))stop("Matrix of times must be the same size as matrix of responses")
			for(i in 1:nrow(response))
				if(any(diff(times[i,])<0,na.rm=T))
				stop(paste("Negative time step for individual ",i))
			times <- as.vector(t(times))}
		nobs <- rep(ncol(response),nrow(response))
		if(!is.null(nest)){
			if(is.vector(nest,mode="numeric")){
				if(length(nest)!=ncol(response))
					stop("Length of nest vector not equal to number of observations per individual")
				else if(any(diff(nest)!=0&diff(nest)!=1,na.rm=T))
					stop("Nest categories must be consecutive increasing integers")
				else tnest <- rep(nest,nrow(response))}
			else if(is.matrix(nest)){
				if(any(dim(nest)!=dim(response)))
					stop("Dimensions of nest not the same as response")
				for(i in nrow(nest))if(any(diff(nest[i,])!=0&diff(nest[i,])!=1,na.rm=T))
					stop("Nest categories must be consecutive increasing integers")
				tnest <- as.vector(t(nest))}
			else stop("nest must be a vector or matrix")}
		if(!is.null(delta)){
			if(is.data.frame(delta))delta <- as.matrix(delta)
			if(is.vector(delta,mode="numeric")){
				if(length(delta)>1){
					if(length(delta)!=ncol(response))
						stop("Length of delta not equal to number of observations per individual")
					else delta <- rep(delta,nrow(response))}}
			else if(is.matrix(delta)){
				if(any(dim(delta)!=dim(response)))
					stop("Dimensions of delta not the same as response")
				delta <- as.vector(t(delta))}
			else stop("delta must be a vector or matrix")}
		y <- as.vector(t(response))}
	else if(is.list(response)){
		if(!is.null(delta)&&!is.vector(delta,mode="numeric")||length(delta)>1)
			stop("delta must be a scalar when response is a list")
		if(!is.null(totals)){
			if(!is.vector(totals,mode="numeric"))stop("totals must be a scalar or vector when response is a list")
			else totv <- T}
		else totv <- F
		if(!is.null(weights)){
			if(is.list(weights))weights <- unlist(weights)
			else if(!is.vector(weights,mode="numeric"))stop("weights must be a list or vector when response is a list")}
		if(is.null(censor)){
			times <- NULL
			tot <- del <- cen <- nes <- 0
			ncols <- ncol(as.matrix(response[[1]]))
			nc <- ttime+1
			if(ncols<nc)stop("matrices must have at least 2 columns: responses and times")
			else if(ncols>nc)for(j in response){
				j <- as.matrix(j)
				for(k in (nc+1):ncols){
					if(is.null(censor)&&any(j[,k]<=0,na.rm=T))
						cen <- k
					else if(any(j[,k]>1,na.rm=T)&&all(j[,k]==trunc(j[,k]),na.rm=T)){
					     if(all(j[,k]>=j[,1],na.rm=T)&&any(diff(j[,k])<0))tot <- k
						else nes <- k}
					else if(is.null(delta)&&all(j[,k]>0,na.rm=T))del <- k}
				if((((ncols==2&&!ttime)||ncols==3)&&(nes>0||cen>0||del>0||tot>0))||(((ncols==3&&!ttime)||ncols==4)&&((nes>0&&cen>0)||(nes>0&&del>0)||(nes>0&&tot>0)||(cen>0&&del>0)))||(((ncols>=4&&!ttime)||ncols>=5)&&((nes>0&&cen>0&&del>0)||(nes>0&&tot>0))))break}
			for(i in response){
				i <- as.matrix(i)
				nind <- nind+1
				if(ncol(i)!=ncols)
					stop(paste("Individual ",nind,"does not have a",ncols,"column matrix"))
				if(nes==0&&ttime&&any(diff(i[,2])<0,na.rm=T))
					stop(paste("Negative time step for individual ",nind))
				nobs <- c(nobs,nrow(i))
				y <- c(y,i[,1])
				if(ttime)times <- c(times,i[,2])
				if(nes>0){
					if(any(diff(i[,nes])!=0&diff(i[,nes])!=1,na.rm=T))
						stop("nest categories for individual ",nind,"are not consecutive increasing integers")
					tnest <- c(tnest,i[,nes])}
				if(!totv&&tot>0)totals <- c(totals,i[,tot])
				if(cen>0)censor <- c(censor,i[,cen])
				if(del>0)delta <- c(delta,i[,del])}}
		else if(!is.vector(censor,mode="numeric"))
			stop("If response is a list, censor must be a vector")
		else {
			del <- nes <- 0
			ncols <- ncol(as.matrix(response[[1]]))
			if(ncols>1){
				for(j in response){
					j <- as.matrix(j)
					for(k in 2:ncols){
						if(is.null(censor)&&any(j[,k]<=0,na.rm=T))cen <- k
						else if(any(j[,k]>1,na.rm=T)&&all(j[,k]==trunc(j[,k]),na.rm=T))nes <- k
						else if(is.null(delta)&&all(j[,k]>0,na.rm=T))del <- k}
					if((ncols==3&&(nes>0||cen>0||del>0))||(ncols==4&&((nes>0&&cen>0)||(nes>0&&del>0)||(cen>0&&del>0)))||(ncols>=5&&((nes>0&&cen>0&&del>0))))break}
				tmp <- NULL
				j <- 0
				for(i in response){
					i <- as.matrix(i)
					nind <- nind+1
					if(ncol(i)!=ncols)stop(paste("Individual ",nind,"does not have a",ncols,"column matrix"))
					nobs <- c(nobs,nrow(i))
					y <- c(y,i[,1])
					tmp <- c(tmp,rep(1,length(i)-1),censor[j <- j+1])
					if(nes>0){
						if(any(diff(i[,nes])!=0&diff(i[,nes])!=1,na.rm=T))stop("nest categories for individual ",nind,"are not consecutive increasing integers")
						tnest <- c(tnest,i[,nes])}
					if(del>0)delta <- c(delta,i[,del])}
				censor <- tmp}
			else {
				tmp <- NULL
				j <- 0
				for(i in response){
					nind <- nind+1
					if(!is.vector(i,mode="numeric")&&!(is.matrix(i)&&ncol(i)==1))
						stop(paste("Individual ",nind,"does not have a vector or one column matrix"))
					tmp <- c(tmp,rep(1,length(i)-1),censor[j <- j+1])
					y <- c(y,i)
					nobs <- c(nobs,length(i))}
				censor <- tmp}}
		if(totv){
			if(length(totals)==1)totals <- rep(totals,length(y))
			else if(length(totals)!=length(y))
			     stop("totals must have one value per response")}}
	else
		stop("Responses must be supplied as a vector, matrix, dataframe, or list of matrices")
	if(!is.null(censor)){
		if(any(censor!=-1&censor!=0&censor!=1,na.rm=T))
			stop("censor must only contain -1, 0, and 1")
		if(is.null(times)){
			j <- 1
			na <- is.na(y)
			y[na] <- 0
			for(i in 1:length(nobs)){
				times <- c(times,cumsum(y[j:(j+nobs[i]-1)]))
				j <- j+nobs[i]}
			y[na] <- NA}
		if(all(censor==1,na.rm=T))censor <- NULL}
	if(!is.null(totals)){
		if(any(y<0,na.rm=T))stop("all responses must be non-negative")
		if(any(totals<y,na.rm=T)||any(totals<0,na.rm=T))
			stop("all totals must be non-negative and >= to responses")}
	if(!is.null(tnest)&&(any(tnest<1,na.rm=T)||any(tnest!=trunc(tnest),na.rm=T)))
		stop("nest must contain integers starting at 1")
	if(!is.null(delta)&&any(delta<=0,na.rm=T))
		stop("delta must be strictly positive")
	if(!is.null(weights)&&any(weights<0,na.rm=T))
		stop("weights must be non-negative")
	type <- if(!is.null(type))match.arg(type,c("nominal","ordinal","discrete","duration","continuous","unknown"))
	else if(!is.null(totals))"nominal"
	else if(!is.null(censor))"duration"
	else if(all(as.integer(y)==y,na.rm=T))"discrete"
	else "unknown"
	z <- list(y=y, nobs=nobs, times=times, nest=tnest,
		coordinates=coordinates, censor=censor, n=totals,
		wt=weights, delta=delta)
	attr(z,"type") <- type
	class(z) <- "response"
	z}

tcctomat <- function(ccov, names=NULL, oldccov=NULL, dataframe=TRUE){
	if(inherits(ccov,"tccov")&&inherits(oldccov,"tccov")){
		if(nrow(ccov$ccov)!=nrow(oldccov$ccov))stop("incompatible tccov objects")
		oldccov$ccov <- cbind(oldccov$ccov,ccov$ccov)
		oldccov$linear <- NULL
		return(oldccov)}
	linear <- NULL
	if(is.language(ccov)){
		linear <- ccov
		mt <- terms(ccov)
		mf <- model.frame(mt,sys.frame(sys.parent()))
		ccov <- model.matrix(mt,mf)[,-1,drop=F]}
	else if(is.factor(ccov)||is.vector(ccov,mode="character")){
		if(is.null(names))names <- paste(deparse(substitute(ccov)))
		ccov <- data.frame(ccov)
		colnames(ccov) <- names}
	if(is.vector(ccov,mode="numeric")){
		if(is.null(names))names <- paste(deparse(substitute(ccov)))
		ccov <- matrix(ccov,ncol=1)}
	else if(is.data.frame(ccov)){
		if(!dataframe){
	                rm(names)
	                tmp <- NULL
	                j <- 0
	                for(i in ccov){
	                	j <- j+1
	                	if(is.vector(i,mode="numeric")){
	                		tmp2 <- as.matrix(i)
	                		colnames(tmp2) <- names(ccov)[j]}
	                	else {
	                		mt <- terms(~i)
	                		tmp2 <- model.matrix(mt,model.frame(mt))[,-1,drop=F]
	                		colnames(tmp2) <- paste(names(ccov)[j],levels(i)[-1],sep="")}
	                	tmp <- cbind(tmp,tmp2)}
	                ccov <- tmp
	                rm(tmp,tmp2,mt)}}
	else if(!is.matrix(ccov))
		stop("Time-constant covariates must be a vector, matrix, dataframe, or model formula")
	if(is.null(colnames(ccov))){
		if(is.null(names))names <- paste(deparse(substitute(ccov)))
		if(length(names)==1&&ncol(ccov)>1)
			names <- paste(names,1:ncol(ccov),sep="")
		colnames(ccov) <- names}
	if(!is.null(oldccov)){
		if(!inherits(oldccov,"tccov"))
			stop("oldccov must have class, tccov")
		else if(nrow(oldccov$ccov)==nrow(ccov)){
			if(dataframe)oldccov$ccov <- data.frame(oldccov$ccov,ccov)
			else oldccov$ccov <- cbind(oldccov$ccov,ccov)}
		else stop("old and new covariates do not have the same number of individuals")}
	else {
		if(dataframe)ccov <- as.data.frame(ccov)
		oldccov <- list(ccov=ccov, linear=linear)
		class(oldccov) <- "tccov"}
	fac <- F
	for(i in 1:ncol(oldccov$ccov))if(!is.vector(oldccov$ccov[,i],mode="numeric")){
		fac <- T
		break}
	if(!fac)oldccov$ccov <- as.matrix(oldccov$ccov)
	oldccov}

tvctomat <- function(tvcov, names=NULL, interaction=NULL, ccov=NULL,
	oldtvcov=NULL, dataframe=TRUE){
	if(inherits(tvcov,"tvcov")&&inherits(oldtvcov,"tvcov")){
		if(length(tvcov$nobs)!=length(oldtvcov$nobs)||any(tvcov$nobs!=oldtvcov$nobs))stop("incompatible tvcov objects")
		oldtvcov$tvcov <- cbind(oldtvcov$tvcov,tvcov$tvcov)
		return(oldtvcov)}
	nbs <- tvcv <- NULL
	if(is.data.frame(tvcov)){
		if(is.null(names))names <- paste(deparse(substitute(tvcov)))
		if(dataframe){
			nbs <- rep(ncol(tvcov),nrow(tvcov))
			tvcv <- as.data.frame(as.vector(t(as.matrix(tvcov))))
			colnames(tvcv) <- names}
		else tvcov <- as.matrix(tvcov)}
	if(is.matrix(tvcov)&&is.character(tvcov)){
		nbs <- rep(ncol(tvcov),nrow(tvcov))
		if(is.null(names))names <- paste(deparse(substitute(tvcov)))
		if(length(names)!=1)stop("too many names")
		tvcv <- as.factor(as.vector(t(tvcov)))
		mt <- terms(~tvcv)
		names <- paste(names,levels(tvcv)[-1],sep="")
		tvcv <- model.matrix(mt,model.frame(mt))[,-1,drop=F]
		colnames(tvcv) <- names}
	else if(is.matrix(tvcov)){
		nbs <- rep(ncol(tvcov),nrow(tvcov))
		tvcv <- matrix(as.vector(t(tvcov)),ncol=1)
		if(!is.null(names)&&length(names)!=1)stop("too many names")
		colnames(tvcv) <- if(!is.null(names)) names
			else paste(deparse(substitute(tvcov)))}
	else if(is.list(tvcov)&&!is.data.frame(tvcov)){
		if(!inherits(tvcov,"tvcov")){
			if(is.null(names)){
				if(!is.null(colnames(tvcov[[1]])))names <- colnames(tvcov[[1]])
				else {
					names <- if(is.matrix(tvcov[[1]]))paste(paste(deparse(substitute(tvcov))),1:ncol(tvcov[[1]]),sep="")
					else paste(deparse(substitute(tvcov)))}}
			for(i in tvcov){
				i <- as.matrix(i)
				nbs <- c(nbs,dim(i)[1])
				tvcv <- rbind(tvcv,i)}
			if(is.character(tvcv)){
				tmp <- tmp2 <- NULL
				for(i in 1:ncol(tvcv)){
					mt <- terms(~as.factor(tvcv[,i]))
					tmp <- cbind(tmp,model.matrix(mt,model.frame(mt))[,-1,drop=F])
					tmp2 <- c(tmp2,paste(names[i],levels(as.factor(tvcv[,i]))[-1],sep=""))}
				tvcv <- tmp
				names <- tmp2
				rm(tmp,tmp2,mt)}
			if(is.null(colnames(tvcov[[1]]))){
				if(is.null(names))names <- paste(deparse(substitute(tvcov)))}
			else if(length(names)!=ncol(tvcv))names <- colnames(tvcov[[1]])
			if(length(names)==1&&ncol(tvcv)>1)
				names <- paste(names,1:ncol(tvcv),sep="")
			colnames(tvcv) <- names}
		else if(inherits(tvcov,"tvcov")){
			nbs <- tvcov$nobs
			tvcv <- tvcov$tvcov}}
	else if(!is.data.frame(tvcov))stop("The time-varying covariates must be a matrix, dataframe, or list")
	if(!is.null(interaction)){
		name <- colnames(tvcv)
		if(is.character(interaction)){
			mat <- match(interaction,name)
			if(any(is.na(mat)))
				stop(paste("Time-varying covariate(s)",ccov[is.na(mat)],"not found"))
			interaction <- mat}
		if(is.vector(interaction,mode="numeric")){
			if(length(interaction)>length(name))stop("too many interactions")
			if(!is.data.frame(tvcv))colnames(tvcv) <- NULL
			if(!is.null(ccov)){
				if(inherits(ccov,"tccov"))ccov <- ccov$ccov
				if(!is.matrix(ccov)&&!is.data.frame(ccov)){
					tmp <- paste(deparse(substitute(ccov)))
					ccov <- matrix(ccov)
					colnames(ccov) <- tmp}
				if(nrow(ccov)!=length(nbs))stop("ccov does not have one observation per individual")
				if(is.null(names))names <- colnames(ccov)
				mat <- match(names,colnames(ccov))
				if(any(is.na(mat)))stop(paste("covariates",names[is.na(mat)],"not found"))
				oldtvcov <- tvcv
				if(!is.data.frame(oldtvcov)&&!is.data.frame(ccov)){
					for(i in 1:length(interaction))for(j in 1:length(mat)){
						oldtvcov <- cbind(oldtvcov,tvcv[,interaction[i]]*rep(ccov[,mat[j]],nbs))
						name <-  c(name,paste(name[interaction[i]],".",names[j],sep=""))}}
				else {
					for(i in 1:length(interaction))for(j in 1:length(mat)){
						mt <- terms(~tvcv[,interaction[i]]:rep(ccov[,mat[j]],nbs))
						tmp <- model.matrix(mt,model.frame(mt))[,-1,drop=F]
				                if(!is.vector(tvcv[,interaction[i]],mode="numeric")){
				                	if(!is.vector(ccov[,mat[j]],mode="numeric")){
				                		nam <- NULL
				                		tmp2 <- paste(name[interaction[i]],levels(tvcv[,interaction[i]])[-1],".",sep="")
				                		for(k in 1:length(levels(ccov[,mat[j]])[-1]))
				                		nam <- c(nam,paste(tmp2,names[j],levels(ccov[,mat[j]])[-1][k],sep=""))
				                		tmp <- tmp[,-c(1:length(levels(tvcv[,interaction[i]])),seq(1,length(levels(tvcv[,interaction[i]]))*length(levels(ccov[,mat[j]])),by=length(levels(tvcv[,interaction[i]])))),drop=F]}
				                	else {
				                		nam <- paste(paste(name[interaction[i]],levels(tvcv[,interaction[i]])[-1],sep=""),".",names[j],sep="")
				                		tmp <- tmp[,-1,drop=F]}}
				                else {
				                	if(!is.vector(ccov[,mat[j]],mode="numeric")){
				                		nam <- paste(name[interaction[i]],".",paste(names[j],levels(ccov[,mat[j]])[-1],sep=""),sep="")
				                		tmp <- tmp[,-1,drop=F]}
				                	else nam <- paste(name[interaction[i]],".",names[j],sep="")}
				                colnames(tmp) <- nam
				                oldtvcov <- cbind(oldtvcov,tmp)}}
				oldtvcov <- list(tvcov=oldtvcov,nobs=nbs)
				if(!is.data.frame(oldtvcov$tvcov))colnames(oldtvcov$tvcov) <- name}
			else if(length(interaction)==2){
				if(is.data.frame(tvcv)){
					mt <- terms(~tvcv[,interaction[1]]:tvcv[,interaction[2]])
					tmp <- model.matrix(mt,model.frame(mt))[,-1,drop=F]
					if(!is.vector(tvcv[,interaction[1]],mode="numeric")){
						if(!is.vector(tvcv[,interaction[2]],mode="numeric")){
							names <- NULL
							tmp2 <- paste(name[interaction[1]],levels(tvcv[,interaction[1]])[-1],".",sep="")
							for(i in 1:length(levels(tvcv[,interaction[2]])[-1]))
							names <- c(names,paste(tmp2,name[interaction[2]],levels(tvcv[,interaction[2]])[-1][i],sep=""))
							tmp <- tmp[,-c(1:length(levels(tvcv[,interaction[1]])),seq(1,length(levels(tvcv[,interaction[1]]))*length(levels(tvcv[,interaction[2]])),by=length(levels(tvcv[,interaction[1]])))),drop=F]}
						else {
							names <- paste(paste(name[interaction[1]],levels(tvcv[,interaction[1]])[-1],sep=""),".",name[interaction[2]],sep="")
							tmp <- tmp[,-1,drop=F]}}
					else {
						if(!is.vector(tvcv[,interaction[2]],mode="numeric")){
							names <- paste(name[interaction[1]],".",paste(name[interaction[2]],levels(tvcv[,interaction[2]])[-1],sep=""),sep="")
							tmp <- tmp[,-1,drop=F]}
						else names <- paste(name[interaction[1]],".",name[interaction[2]],sep="")}
					colnames(tmp) <- names
					oldtvcov <- list(tvcov=cbind(tvcv,tmp),nobs=nbs)}
				else oldtvcov <- list(tvcov=cbind(tvcv,tvcv[,interaction[1]]*tvcv[,interaction[2]]),nobs=nbs)
				if(!is.data.frame(oldtvcov$tvcov))colnames(oldtvcov$tvcov) <- c(name,paste(name[interaction[1]],".",name[interaction[2]],sep=""))}
			else stop("interaction must be a vector containing column numbers or variable names")
			class(oldtvcov) <- "tvcov"}}
	else if(!is.null(oldtvcov)){
		if(!inherits(oldtvcov,"tvcov"))
			stop("oldtvcov must have class, tvcov")
		else if((nrow(oldtvcov$tvcov)==nrow(tvcv))&&
			all(oldtvcov$nobs==nbs)){
			if(dataframe)oldtvcov$tvcov <- data.frame(oldtvcov$tvcov,tvcv)
			else oldtvcov$tvcov <- cbind(oldtvcov$tvcov,tvcv)}
		else stop("old and new covariates do not have the same numbers of observations")}
	else {
		if(dataframe)tvcv <- as.data.frame(tvcv)
		oldtvcov <- list(tvcov=tvcv,nobs=nbs)
		class(oldtvcov) <- "tvcov"}
	fac <- F
	for(i in 1:ncol(oldtvcov$tvcov))if(!is.vector(oldtvcov$tvcov[,i],mode="numeric")){
		fac <- T
		break}
	if(!fac)oldtvcov$tvcov <- as.matrix(oldtvcov$tvcov)
	oldtvcov}

rmna <- function(response, ccov=NULL, tvcov=NULL){
	if(!inherits(response,"response"))response <- restovec(response)
	if(length(response$nobs)==1&&response$nobs==1)
		response$nobs <- rep(1,length(response$y))
	if(!is.null(ccov)){
		if(!inherits(ccov,"tccov"))ccov <- tcctomat(ccov)
		if(length(response$nobs)!=nrow(ccov$ccov))stop("Numbers of individuals for response and for time-constant covariates do not agree.")}
	if(!is.null(tvcov)){
		if(!inherits(tvcov,"tvcov"))tvcov <- tvctomat(tvcov)
		if(any(response$nobs!=tvcov$nobs))stop("Numbers of observations for response and time-varying covariates do not agree.")}
	rna <- !is.na(response$y)
	if(!is.null(response$times))rna <- rna&!is.na(response$times)
	if(!is.null(response$nest))rna <- rna&!is.na(response$nest)
	if(!is.null(response$coordinates))rna <- rna&!is.na(response$coordinates[,1])&!is.na(response$coordinates[,2])
	if(!is.null(response$n))rna <- rna&!is.na(response$n)
	for(i in 1:length(response$nobs))
		if(!is.null(ccov)&&any(is.na(ccov$ccov[i,])))rna[covind(response)==i] <- F
	if(!is.null(tvcov))
		for(i in 1:ncol(tvcov$tvcov))rna <- rna&!is.na(tvcov$tvcov[,i])
	response$y <- response$y[rna]
	if(!is.null(response$times))response$times <- response$times[rna]
	if(!is.null(response$nest))response$nest <- response$nest[rna]
	if(!is.null(response$coordinates))response$coordinates <- response$coordinates[rna,]
	if(!is.null(response$n))response$n <- response$n[rna]
	if(!is.null(response$censor)){
		response$censor <- response$censor[rna]
		if(all(response$censor==1))response$censor <- NULL}
	if(!is.null(response$delta)&&length(response$delta)>1)
		response$delta <- response$delta[rna]
	if(!is.null(tvcov))tvcov$tvcov <- tvcov$tvcov[rna,,drop=F]
	tmp <- NULL
	j <- c(0,cumsum(response$nobs))
	for(i in 1:length(response$nobs)){
		tmp <- c(tmp,sum(rna[(j[i]+1):j[i+1]]))
		if(tmp[i]==0)
			warning(paste("Individual",i,"has no observations"))}
	response$nobs <- tmp[tmp>0]
	if(!is.null(ccov)){
		ccov$ccov <- ccov$ccov[tmp>0,,drop=F]
		for(i in 1: ncol(ccov$ccov))
			if(length(unique(ccov$ccov[,i]))==1)
			warning(paste("covariate",colnames(ccov$ccov)[i],"has only one value\n"))}
	if(!is.null(tvcov)){
		tvcov$nobs <- response$nobs
		for(i in 1: ncol(tvcov$tvcov))
			if(length(unique(tvcov$tvcov[,i]))==1)
			warning(paste("covariate",colnames(tvcov$tvcov)[i],"has only one value\n"))}
	if(all(response$nobs==1))response$nobs <- 1
	z <- list(response=response,tvcov=tvcov,ccov=ccov)
	class(z) <- "repeated"
	z}

print.response <- function(z){
	if(length(z$nobs)>1||z$nobs>1)cat("Number of individuals:                ",length(z$nobs),"\n")
	cat("Number of observations:               ",length(z$y),"\n")
	if(length(z$nobs)>1)cat("Number of observations per individual:\n",z$nobs,"\n")
	if(is.null(z$n))y <- z$y
	else y <- z$y/z$n
	cat("Mean response:                        ",mean(y,na.rm=T),"\n")
	cat("Range of responses:                   ",range(y,na.rm=T),"\n")
	if(any(is.na(z$y)))
		cat("Number of NAs:                        ",sum(is.na(z$y)),"\n")
	if(!is.null(z$wt))
		cat("Number of positive weights:           ",sum(z$wt>0),"\n")
	if(!is.null(z$times)){
		cat("Mean time:                            ",mean(z$times,na.rm=T),"\n")
		cat("Range of times:                       ",range(z$times,na.rm=T),"\n")
		if(length(z$nobs)>1&&is.null(z$nest)){
			mn <- if(any(z$times<0))z$times[cumsum(c(1,z$nobs[1:(length(z$nobs)-1)]))]
				else 0
			cat("Mean total time:                      ",mean(z$times[cumsum(z$nobs)]-mn,na.rm=T),"\n")
			cat("Range of total times:                 ",range(z$times[cumsum(z$nobs)]-mn,na.rm=T),"\n")}}
	if(!is.null(z$nest))
		cat("Maximum number of clusters:           ",max(z$nest),"\n")
	if(!is.null(z$censor)) {
		if(sum(z$censor==0,na.rm=T)>0)cat("Number of right-censored observations:",sum(z$censor==0,na.rm=T),"\n")
		if(sum(z$censor==-1,na.rm=T)>0)cat("Number of left-censored observations: ",sum(z$censor==-1,na.rm=T),"\n")}
	if(!is.null(z$delta)&&length(z$delta)==1)
		cat("Unit of measurement:                  ",z$delta,"\n")}

print.tvcov <- function(z){
	if(is.function(z)){
		print.default(unclass(z))
		return()}
	cat("Number of individuals:            ",length(z$nobs),"\n")
	cat("Number of observations:           ",sum(z$nobs),"\n")
	cat("Number of observations per individual:\n",z$nobs,"\n")
	cat("Number of time-varying covariates:",ncol(z$tvcov),"\n")
	cat("Names of time-varying covariates:\n",colnames(z$tvcov),"\n")}

print.tccov <- function(z){
	if(is.function(z)){
		print.default(unclass(z))
		return()}
	cat("Number of individuals:             ",nrow(z$ccov),"\n")
	cat("Number of time-constant covariates:",ncol(z$ccov),"\n")
	cat("Names of time-constant covariates:\n",colnames(z$ccov),"\n")}

print.repeated <- function(z){
	if(is.function(z)){
		print.default(unclass(z))
		return()}
	cat("Response variable:\n\n")
	print.response(z$response)
	if(!is.null(z$ccov)){
		cat("\nTime-constant covariates:\n\n")
		print.tccov(z$ccov)}
	if(!is.null(z$tvcov)){
		cat("\nTime-varying covariates:\n\n")
		print.tvcov(z$tvcov)}}

plot.response <- function(z, nind=NULL, nest=1, add=F, lty=NULL,
	pch=NULL, main=NULL, ylim=range(z$y,na.rm=T),
	xlim=range(z$times,na.rm=T), xlab=NULL, ylab="Response", ...){
	if(length(z$nobs)==1&&z$nobs==1)z$nobs <- length(z$y)
	if(is.null(z$times)){
		z$times <- sequence(z$nobs)
		if(is.null(xlab)) xlab <- "Index number"}
	else if(is.null(xlab)) xlab <- "Time"
	tnest <- if(!is.null(z$nest)) z$nest
		else 1
	nm <- rep(1:length(z$nobs),z$nobs)
	j <- 1
	lt <- 0
	if(is.null(nind))nind <- 1:length(z$nobs)
	if(is.null(z$n))y <- z$y
	else y <- z$y/z$n
	if(!is.null(lty)){
		if(length(lty)==1)lty <- rep(lty,length(nind))
		else if(length(lty)!=length(nind))stop("lty must have one value for each item in nind")}
	if(!is.null(pch)){
		if(length(pch)==1)pch <- rep(pch,length(nind))
		else if(length(pch)!=length(nind))stop("pch must have one value for each item in nind")}
	for(i in 1:length(z$nobs))if(any(i==nind)){
		if(is.null(lty))lt <- lt%%4+1
		else lt <- lty[j]
		if(!add&&j==1)plot(z$times[nm==i&nest==tnest],
			y[nm==i&nest==tnest],lty=lt,type="l",
			ylim=ylim,xlim=xlim,main=main,
			ylab=ylab,xlab=xlab,...)
		else lines(z$times[nm==i&nest==tnest],
			y[nm==i&nest==tnest],lty=lt)
		if(!is.null(pch))points(z$times[nm==i&nest==tnest],
			y[nm==i&nest==tnest],pch=pch[j])
		j <- j+1}}

plot.repeated <- function(z, variable="response", number=1, nind=NULL,
	add=F, lty=NULL, main=NULL, ylim=range(z$response$y),
	xlim=range(z$response$times), xlab="Time", ylab="Response", ...){
	variable <- match.arg(variable,c("response","time-varying covariate"))
	if(inherits(z,"repeated")){
		if(variable=="response")
			plot.response(z$response, nind=nind,
				add=add, lty=lty, main=main, ylim=ylim,
				xlim=xlim, xlab=xlab, ylab=ylab, ...)
		else if(variable=="time-varying covariate"){
			if(number>ncol(z$tvcov$tvcov))
				stop("Less than",number,"covariates")
			if(missing(ylab))ylab <- colnames(z$tvcov$tvcov)[number]
			zz <- list()
			zz$times <- z$response$times
			zz$y <- z$tvcov$tvcov[,number]
			zz$nobs <- z$tvcov$nobs
			if(missing(ylim))ylim <- range(zz$y)
			class(zz) <- "response"
			plot.response(zz, nind=nind,
				add=add, lty=lty, main=main, ylim=ylim,
				xlim=xlim, xlab=xlab, ylab=ylab, ...)}}}

response <- function(z, ...) UseMethod("response")

response.response <- function(z, nind=NULL){
	if(is.null(nind))nind <- 1:length(z$y)
	else if(length(nind)>length(z$nobs)||any(nind>length(z$nobs)))
		stop("Individual not found")
	else nind <- !is.na(match(covind(z),nind))
	if(all(!nind))stop("No such individuals")
	if(!is.null(z$n)){
		z <- cbind(z$y,z$n-z$y)[nind,]
		colnames(z) <- c("y","n-y")
		return(z)}
	else if(!is.null(z$censor)){
		z <- cbind(z$y,z$censor)[nind,]
		colnames(z) <- c("Response","Censor")
		return(z)}
	else return(z$y[nind])}

response.repeated <- function(z, nind=NULL){
	if(is.null(nind))nind <- 1:length(z$response$y)
	else if(length(nind)>length(z$response$nobs)||any(nind>length(z$response$nobs)))
		stop("Individual not found")
	else nind <- !is.na(match(covind(z),nind))
	if(all(!nind))stop("No such individuals")
	if(!is.null(z$response$n)){
		z <- cbind(z$response$y,z$response$n-z$response$y)[nind,]
		colnames(z) <- c("y","n-y")
		return(z)}
	else if(!is.null(z$response$censor)){
		z <- cbind(z$response$y,z$response$censor)[nind,]
		colnames(z) <- c("Response","Censor")
		return(z)}
	else return(z$response$y[nind])}

covind <- function(z, ...) UseMethod("covind")

covind.response <- function(z) {
	if(length(z$nobs)==1&&z$nobs==1) return(1:length(z$y))
	else return(rep(1:length(z$nobs),z$nobs))}

covind.repeated <- function(z) {
	if(length(z$response$nobs)==1&&z$response$nobs==1)
		return(1:length(z$response$y))
	else return(rep(1:length(z$response$nobs),z$response$nobs))}

nobs <- function(z, ...) UseMethod("nobs")

nobs.response <- function(z) z$nobs

nobs.repeated <- function(z) z$response$nobs

times <- function(z, ...) UseMethod("times")

times.response <- function(z) z$times

times.repeated <- function(z) z$response$times

weights.response <- function(z) z$wt

weights.repeated <- function(z) z$response$wt

nesting <- function(z, ...) UseMethod("nesting")

nesting.response <- function(z){
	if(length(z$nobs)==1)return(NULL)
	else if(is.null(z$nest))return(covind(z))
	else {
		z <- cbind(covind(z),z$nest)
		colnames(z) <- c("Individual","Cluster")
		return(z)}}

nesting.repeated <- function(z){
	if(length(z$response$nobs)==1)return(NULL)
	else if(is.null(z$response$nest))return(covind(z))
	else {
		z <- cbind(covind(z),z$response$nest)
		colnames(z) <- c("Individual","Cluster")
		return(z)}}

covariates <- function(z, ...) UseMethod("covariates")

covariates.tccov <- function(z, names=NULL){
	if(is.null(names))return(z$ccov)
	else {
		mat <- match(names,colnames(z$ccov))
		if(any(is.na(mat)))
			stop(paste("covariate(s)",names[is.na(mat)],"not found"))
		return(z$ccov[,mat,drop=F])}}

covariates.tvcov <- function(z, names=NULL){
	if(is.null(names))return(z$tvcov)
	else {
		mat <- match(names,colnames(z$tvcov))
		if(any(is.na(mat)))
			stop(paste("covariate(s)",names[is.na(mat)],"not found"))
		return(z$tvcov[,mat,drop=F])}}

covariates.repeated <- function(z, names=NULL){
	if(is.null(names)){
		if(is.null(z$tvcov$tvcov))return(z$ccov$ccov)
		else if(is.null(z$ccov$ccov))return(z$tvcov$tvcov)
		else return(list(ccov=z$ccov$ccov,tvcov=z$tvcov$tvcov))}
	else {
		mat1 <- match(names,colnames(z$ccov$ccov))
		mat2 <- match(names,colnames(z$tvcov$tvcov))
		if(any(is.na(mat1)&is.na(mat2)))
			stop(paste("covariate(s)",names[is.na(mat1)&is.na(mat2)],"not found"))
		if(all(is.na(mat1)))return(z$tvcov$tvcov[,mat2,drop=F])
		else if(all(is.na(mat2)))return(z$ccov$ccov[,mat1,drop=F])
		else return(list(ccov=z$ccov$ccov[,mat1,drop=F],tvcov=z$tvcov$tvcov[,mat2,drop=F]))}}

names.tccov <- function(z) colnames(z$ccov)

names.tvcov <- function(z) colnames(z$tvcov)

names.repeated <- function(z)
	list(ccov=colnames(z$ccov$ccov),tvcov=colnames(z$tvcov$tvcov))

formula.tccov <- function(z) z$linear

formula.repeated <- function(z) z$ccov$linear

transform.response <- function(z, y=NULL, times=NULL){
	if(is.call(substitute(y)))y <- substitute(y)
	if(is.call(substitute(times)))times <- substitute(times)
	if(!is.null(y)){
		if(!is.null(z$n))stop("transformations do not make sense with binomial data")
		tran <- eval(deriv(y,"y"),z,NULL)
		z$y <- as.vector(tran)
		z$delta <- if(is.null(z$delta))as.vector(abs(attr(tran,"gradient")))
			else z$delta*as.vector(abs(attr(tran,"gradient")))
		if(all(z$delta==1))z$delta <- NULL}
	if(!is.null(times)){
		z$times <- eval(times,z,NULL)
		for(i in 1:length(z$nobs))if(any(diff(z$times[covind(z)==i])<0))stop("transformation produces negative time steps")}
	z}

transform.repeated <- function(z, y=NULL, times=NULL){
	if(is.call(substitute(y)))y <- substitute(y)
	if(is.call(substitute(times)))times <- substitute(times)
	z$response <- transform.response(z$response,y,times)
	z}

transform.tccov <- function(z, ...){
	isf <- is.data.frame(z$ccov)
	if(!isf)z$ccov <- as.data.frame(z$ccov)
	e <- eval(substitute(list(...)),z$ccov,NULL)
	tags <- names(e)
	for(i in 1:length(e))if(all(is.na(e[[i]])))stop(paste(tags[i],"defines an invalid tranformation\n or attempts to transform a factor variable"))
	inx <- match(tags,colnames(z$ccov))
	matched <- !is.na(inx)
	if(any(matched))z$ccov[inx[matched]] <- e[matched]
	if(!all(matched))z$ccov <- data.frame(z$ccov,e[!matched])
	if(!isf)z$ccov <- as.matrix(z$ccov)
	z}

transform.tvcov <- function(z, ...){
	isf <- is.data.frame(z$tvcov)
	if(!isf)z$tvcov <- as.data.frame(z$tvcov)
	e <- eval(substitute(list(...)),z$tvcov,NULL)
	tags <- names(e)
	for(i in 1:length(e))if(all(is.na(e[[i]])))stop(paste(tags[i],"defines an invalid tranformation\n or attempts to transform a factor variable"))
	inx <- match(tags,colnames(z$tvcov))
	matched <- !is.na(inx)
	if(any(matched))z$tvcov[inx[matched]] <- e[matched]
	if(!all(matched))z$tvcov <- data.frame(z$tvcov,e[!matched])
	if(!isf)z$tvcov <- as.matrix(z$tvcov)
	z}
