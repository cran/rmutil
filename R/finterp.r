#
#  rmutil : A Library of Special Functions for Repeated Measurements
#  Copyright (C) 1999 J.K. Lindsey
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
#   finterp(z, envir=sys.frame(sys.parent()), formula=FALSE,
#	parameters=FALSE, start=1, name=NULL, expand=NULL)
#   fnenvir(z, envir=sys.frame(sys.parent()), name=NULL,
#	expand=TRUE)
#
#  DESCRIPTION
#
#    Function to translate a model formula with unknown parameters
#  into a function.

finterp <- function(z, ...) UseMethod("finterp")

finterp.default <- function(z, envir=sys.frame(sys.parent()), formula=FALSE,
	vector=TRUE, start=1, name=NULL, expand=NULL){
	if(!inherits(z,"formula"))return(NULL)
	if(!is.environment(envir)){
		if(inherits(envir,"repeated"))return(finterp.repeated(z,envir,name,expand))
		if(inherits(envir,"tccov"))return(finterp.tccov(z,envir,name,expand))
		if(inherits(envir,"tvcov"))return(finterp.tvcov(z,envir,name,expand))}
	ch1 <- deparse(z[[length(z)]])
	i <- un <- fcn <- ex <- ch <- par <- NULL
	for(i in 1:length(ch1))ch <- paste(ch,ch1[i],collapse=" ")
	mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(^[0-9]+)"," ",ch)," ")[[1]]
	if(length(mem)>0){
		for(i in 1:length(mem)){
			ex <- c(ex,exists(mem[i],envir=envir))
			fcn <- c(fcn,if(exists(mem[i]))
				is.function(eval(parse(text=mem[i])))
				else F)}
		un <- unique(mem[!ex])
		if(length(unique(mem[ex&!fcn]))==0&&length(un)==0)
			stop("finterp: no variables found")}
	if(is.null(ex)||all(ex)){
		if(formula)return(z)
		else {
			mt <- terms(z)
			if(is.numeric(mt[[2]])){
				dm <- matrix(1)
				colnames(dm) <- "(Intercept)"}
			else dm <- model.matrix(mt,model.frame(mt,envir))
			.fn <- function(.p) as.vector(attr(.fn,"model")%*%
				.p[attr(.fn,"range")[1]:attr(.fn,"range")[2]])
			attributes(.fn) <- list(formula=z,model=dm,
				covariates=if(length(mem)>0)
					unique(mem[ex&!fcn]) else NULL,
				parameters=paste("p[",1:ncol(dm),"]",sep=""),
				range=c(start,start+ncol(dm)-1),
				class="formulafn")
			rm(list=ls())
			return(.fn)}}
	.fn <- function(.p) eval(attr(.fn,"model"))
	if(vector){
		ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch)))
		ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch)))
		ch <- paste(" ",gsub("\\("," ( ",ch)," ",sep="")
		for(i in 1:length(un))ch <- gsub(paste(" ",un[i]," ",sep=""),
			paste(" .p[",start+i-1,"] ",sep=""),ch)}
	else {
		par <- "alist("
		for(i in 1:length(un)){
			if(i>1)par <- paste(par,",",collapse="")
			par <- paste(par,un[i],"=",collapse="")}
# bug in 0.64.0
		if(length(un)==1)par <- paste(par,",...=",collapse="")
		par <- paste(par,")",collapse="")
		formals(.fn) <- eval(parse(text=par))}
	attributes(.fn) <- list(formula=z,model=parse(text=ch),parameters=un,
		covariates=unique(mem[ex&!fcn]),
		range=c(start,start+length(un)-1),class="formulafn")
	rm(list=ls())
	return(.fn)}

finterp.repeated <- function(z, envir=NULL, formula=FALSE, start=1,
	name=NULL, expand=TRUE){
	if(!inherits(z,"formula"))return(NULL)
	if(is.null(envir)||!inherits(envir,"repeated"))stop("envir must be an object of class, repeated")
	ndata <- if(is.null(name))paste(deparse(substitute(envir))) else name
	ch1 <- deparse(z[[length(z)]])
	i <- un <- fcn <- ex1 <- ex2 <- ch <- NULL
	for(i in 1:length(ch1))ch <- paste(ch,ch1[i],collapse=" ")
	mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(^[0-9]+)"," ",ch)," ")[[1]]
	if(length(mem)>0){
		ex1 <- match(mem,colnames(envir$ccov$ccov))
		ex2 <- match(mem,colnames(envir$tvcov$tvcov))
		if(any(!is.na(ex2))&&!expand)stop("time-varying covariates present - time-constant ones must be expanded")
		for(i in 1:length(mem)){
			fcn <- c(fcn,if(exists(mem[i]))
				is.function(eval(parse(text=mem[i])))
				else F)}
		un <- unique(mem[is.na(ex1)&is.na(ex2)&!fcn])
		if(length(unique(mem[(!is.na(ex1)|!is.na(ex2))&!fcn]))==0&&length(un)==0)
			stop("finterp: no variables found")}
	ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch)))
	ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch)))
	ch <- paste(" ",gsub("\\("," ( ",ch)," ",sep="")
	if(expand).i <- covind(envir)
	ex1a <- ex1[!is.na(ex1)]
	if(length(ex1a)>0)for(i in 1:length(ex1a))
		ch <- if(expand)gsub(paste(" ",colnames(envir$ccov$ccov)[ex1a[i]],
			" ",sep=""),paste(" ",ndata,"$ccov$ccov[,",
			ex1a[i],"][.i] ",sep=""),ch)
			else gsub(paste(" ",colnames(envir$ccov$ccov)[ex1a[i]],
			" ",sep=""),paste(" ",ndata,"$ccov$ccov[,",
			ex1a[i],"] ",sep=""),ch)
	ex2a <- ex2[!is.na(ex2)]
	if(length(ex2a)>0)for(i in 1:length(ex2a))
		ch <- gsub(paste(" ",colnames(envir$tvcov$tvcov)[ex2a[i]],
			" ",sep=""),paste(" ",ndata,"$tvcov$tvcov[,",
			ex2a[i],"] ",sep=""),ch)
	if(all(!is.na(ex1)|!is.na(ex2))){
		if(formula)return(z)
		else {
			ch <- as.formula(paste("~",ch))
			mt <- terms(ch)
			if(is.numeric(mt[[2]]))dm <- matrix(1)
			else dm <- model.matrix(mt,model.frame(mt))
			.fn <- function(.p) as.vector(attr(.fn,"model")%*%
				.p[attr(.fn,"range")[1]:attr(.fn,"range")[2]])
			attributes(.fn) <- list(formula=z,model=dm,
				covariates=if(length(mem)>0)
					unique(mem[(!is.na(ex1)|!is.na(ex2))&!fcn]) else NULL,
				parameters=paste("p[",1:ncol(dm),"]",sep=""),
				range=c(start,start+ncol(dm)-1),
				class="formulafn")
			rm(list=ls())
			return(.fn)}}
	.fn <- function(.p) eval(attr(.fn,"model"))
	for(i in 1:length(un))ch <- gsub(paste(" ",un[i]," ",sep=""),
		paste(" .p[",start+i-1,"] ",sep=""),ch)
	attributes(.fn) <- list(formula=z,model=parse(text=ch),parameters=un,
		covariates=unique(mem[(!is.na(ex1)|!is.na(ex2))&!fcn]),
		range=c(start,start+length(un)-1),class="formulafn")
	rm(list=ls())
	return(.fn)}

finterp.tccov <- function(z, envir=NULL, formula=FALSE, start=1,
	name=NULL, expand=NULL){
	if(!inherits(z,"formula"))return(NULL)
	if(is.null(envir)||(!inherits(envir,"repeated")&&!inherits(envir,"tccov")))stop("envir must be an object of class, repeated or tccov")
	ndata <- if(is.null(name))paste(deparse(substitute(envir))) else name
	ch1 <- deparse(z[[length(z)]])
	i <- un <- fcn <- ex1 <- ch <- NULL
	for(i in 1:length(ch1))ch <- paste(ch,ch1[i],collapse=" ")
	mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(^[0-9]+)"," ",ch)," ")[[1]]
	if(inherits(envir,"repeated")){
		cn <- colnames(envir$ccov$ccov)
		cc <- "$ccov"}
	else {
		cn <- colnames(envir$ccov)
		cc <- ""}
	if(length(mem)>0){
		ex1 <- match(mem,cn)
		for(i in 1:length(mem)){
			fcn <- c(fcn,if(exists(mem[i]))
				is.function(eval(parse(text=mem[i])))
				else F)}
		un <- unique(mem[is.na(ex1)&!fcn])
		if(length(unique(mem[!is.na(ex1)&!fcn]))==0&&length(un)==0)
			stop("finterp: no variables found")}
	ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch)))
	ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch)))
	ch <- paste(" ",gsub("\\("," ( ",ch)," ",sep="")
	ex1a <- ex1[!is.na(ex1)]
	if(length(ex1a)>0)for(i in 1:length(ex1a))
		ch <- gsub(paste(" ",cn[ex1a[i]]," ",sep=""),
			paste(" ",ndata,cc,"$ccov[,",
			ex1a[i],"] ",sep=""),ch)
	if(all(!is.na(ex1))){
		if(formula)return(z)
		else {
			ch <- as.formula(paste("~",ch))
			mt <- terms(ch)
			if(is.numeric(mt[[2]]))dm <- matrix(1)
			else dm <- model.matrix(mt,model.frame(mt))
			.fn <- function(.p) as.vector(attr(.fn,"model")%*%
				.p[attr(.fn,"range")[1]:attr(.fn,"range")[2]])
			attributes(.fn) <- list(formula=z,model=dm,
				covariates=if(length(mem)>0)
					unique(mem[!is.na(ex1)&!fcn]) else NULL,
				parameters=paste("p[",1:ncol(dm),"]",sep=""),
				range=c(start,start+ncol(dm)-1),
				class="formulafn")
			rm(list=ls())
			return(.fn)}}
	.fn <- function(.p) eval(attr(.fn,"model"))
	for(i in 1:length(un))ch <- gsub(paste(" ",un[i]," ",sep=""),
		paste(" .p[",start+i-1,"] ",sep=""),ch)
	attributes(.fn) <- list(formula=z,model=parse(text=ch),parameters=un,
		covariates=unique(mem[!is.na(ex1)&!fcn]),
		range=c(start,start+length(un)-1),class="formulafn")
	rm(list=ls())
	return(.fn)}

finterp.tvcov <- function(z, envir=NULL, formula=FALSE, start=1,
	name=NULL, expand=NULL){
	if(!inherits(z,"formula"))return(NULL)
	if(is.null(envir)||(!inherits(envir,"repeated")&&!inherits(envir,"tvcov")))stop("envir must be an object of class, repeated or tvcov")
	ndata <- if(is.null(name))paste(deparse(substitute(envir))) else name
	ch1 <- deparse(z[[length(z)]])
	i <- un <- fcn <- ex1 <- ch <- NULL
	for(i in 1:length(ch1))ch <- paste(ch,ch1[i],collapse=" ")
	mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(^[0-9]+)"," ",ch)," ")[[1]]
	if(inherits(envir,"repeated")){
		cn <- colnames(envir$tvcov$tvcov)
		cc <- "$tvcov"}
	else {
		cn <- colnames(envir$tvcov)
		cc <- ""}
	if(length(mem)>0){
		ex1 <- match(mem,cn)
		for(i in 1:length(mem)){
			fcn <- c(fcn,if(exists(mem[i]))
				is.function(eval(parse(text=mem[i])))
				else F)}
		un <- unique(mem[is.na(ex1)&!fcn])
		if(length(unique(mem[!is.na(ex1)&!fcn]))==0&&length(un)==0)
			stop("finterp: no variables found")}
	ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch)))
	ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch)))
	ch <- paste(" ",gsub("\\("," ( ",ch)," ",sep="")
	ex1a <- ex1[!is.na(ex1)]
	if(length(ex1a)>0)for(i in 1:length(ex1a))
		ch <- gsub(paste(" ",cn[ex1a[i]]," ",sep=""),
			paste(" ",ndata,cc,"$tvcov[,",
			ex1a[i],"] ",sep=""),ch)
	if(all(!is.na(ex1))){
		if(formula)return(z)
		else {
			ch <- as.formula(paste("~",ch))
			mt <- terms(ch)
			if(is.numeric(mt[[2]]))dm <- matrix(1)
			else dm <- model.matrix(mt,model.frame(mt))
			.fn <- function(.p) as.vector(attr(.fn,"model")%*%
				.p[attr(.fn,"range")[1]:attr(.fn,"range")[2]])
			attributes(.fn) <- list(formula=z,model=dm,
				covariates=if(length(mem)>0)
					unique(mem[!is.na(ex1)&!fcn]) else NULL,
				parameters=paste("p[",1:ncol(dm),"]",sep=""),
				range=c(start,start+ncol(dm)-1),
				class="formulafn")
			rm(list=ls())
			return(.fn)}}
	.fn <- function(.p) eval(attr(.fn,"model"))
	for(i in 1:length(un))ch <- gsub(paste(" ",un[i]," ",sep=""),
		paste(" .p[",start+i-1,"] ",sep=""),ch)
	attributes(.fn) <- list(formula=z,model=parse(text=ch),parameters=un,
		covariates=unique(mem[!is.na(ex1)&!fcn]),
		range=c(start,start+length(un)-1),class="formulafn")
	rm(list=ls())
	return(.fn)}

fnenvir <- function(z, ...) UseMethod("fnenvir")

fnenvir.default <- function(z, envir=sys.frame(sys.parent()),
	name=NULL, expand=TRUE){
	if(!is.function(z))return(NULL)
	if(!is.environment(envir)){
		if(inherits(envir,"repeated"))return(fnenvir.repeated(z,envir,name,expand))
		if(inherits(envir,"tccov"))return(fnenvir.tccov(z,envir,name,expand))
		if(inherits(envir,"tvcov"))return(fnenvir.tvcov(z,envir,name,expand))}
	ch1 <- deparse(z,width=500)
	ch2 <- ch1[1]
	ch1 <- ch1[2:length(ch1)]
	mem2 <- strsplit(gsub("[(),]"," ",ch2)," ")[[1]]
	if(length(mem2)>1)mem2 <- mem2[2:length(mem2)]
	else mem2 <- NULL
	i <- un <- fcn <- ex <- ch <- NULL
	for(i in 1:length(ch1))ch <- paste(ch,ch1[i],collapse=" ")
	mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(^[0-9]+)"," ",ch)," ")[[1]]
	if(length(mem)>0){
		for(i in 1:length(mem)){
			ex <- c(ex,exists(mem[i],envir=envir))
			fcn <- c(fcn,if(exists(mem[i])){
				if(mem[i]=="function"||mem[i]=="if"||
					mem[i]=="else"||mem[i]=="for"||
					mem[i]=="while"||mem[i]=="repeat") T
				else is.function(eval(parse(text=mem[i])))}
				else F)}
		un <- unique(mem[!ex])
		if(length(unique(mem[ex&!fcn]))==0&&length(un)==0)
			stop("finterp: no variables found")}
	ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch)))
	ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch)))
	ch <- paste(" ",gsub("\\("," ( ",ch)," ",sep="")
	ch2 <- strsplit(ch," ")[[1]]
	ex2a <- un <- un0 <- un1 <- NULL
	if(length(mem2)>0)for(i in 1:length(mem2)){
		ex1a <- NULL
		for(j in 1:length(ch2))
			if(mem2[i]==ch2[j]){
				if(j<length(ch2)&&length(grep("^\\[",ch2[j+1]))>0){
					ex1a <- c(ex1a,paste(ch2[j],ch2[j+1],sep=""))
					un1 <- c(un1,ch2[j])}
				else un0 <- c(un0,ch2[j])}
		if(!is.null(ex1a)){
			ex1a <- unique(ex1a)
			ex2a <- c(ex2a,length(ex1a))
			o <- gsub("(^[[:alnum:]]\\[)|(\\])","",ex1a)
			un <- if(length(grep("[[:alpha:]]",o))>0)c(un,ex1a)
				else c(un,ex1a[order(as.numeric(o))])}}
	if(length(un0)>0){
		un <- if(length(un1)>0&&length(grep(un1,un0))>0)
			c(un,unique(un0[-grep(un1,un0)]))
			else c(un,unique(un0))}
	.fn <- eval(parse(text=paste("function(",paste(mem2,collapse=","),")",paste("eval(attr(.fn,\"model\"))"))))
	attributes(.fn) <- list(model=parse(text=ch1),parameters=un,
		covariates=unique(mem[ex&!fcn]),class="formulafn")
	rm(list=ls())
	return(.fn)}

fnenvir.repeated <- function(z, envir=NULL, name=NULL, expand=TRUE){
	if(!is.function(z))return(NULL)
	if(is.null(envir)||!inherits(envir,"repeated"))stop("envir must be an object of class, repeated")
	ndata <- if(is.null(name))paste(deparse(substitute(envir)))
	else name
	ch1 <- deparse(z,width=500)
	ch2 <- ch1[1]
	ch1 <- ch1[2:length(ch1)]
	mem2 <- strsplit(gsub("[(),]"," ",ch2)," ")[[1]]
	if(length(mem2)>1)mem2 <- mem2[2:length(mem2)]
	else mem2 <- NULL
	i <- un <- fcn <- ex1 <- ex2 <- ch <- NULL
	for(i in 1:length(ch1))ch <- paste(ch,ch1[i],collapse=" ")
	mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(^[0-9]+)"," ",ch)," ")[[1]]
	if(length(mem)>0){
		ex1 <- match(mem,colnames(envir$ccov$ccov))
		ex2 <- match(mem,colnames(envir$tvcov$tvcov))
		if(any(!is.na(ex2))&&!expand)stop("time-varying covariates present - time-constant ones must be expanded")
		for(i in 1:length(mem)){
			fcn <- c(fcn,if(exists(mem[i])){
				if(mem[i]=="function"||mem[i]=="if"||
					mem[i]=="else"||mem[i]=="for"||
					mem[i]=="while"||mem[i]=="repeat") T
				else is.function(eval(parse(text=mem[i])))}
				else F)}
		un <- unique(mem[is.na(ex1)&is.na(ex2)&!fcn])
		if(length(unique(mem[(!is.na(ex1)|!is.na(ex2))&!fcn]))==0&&length(un)==0)
			stop("finterp: no variables found")}
	for(i in 1:length(ch1)){
		ch1[i] <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch1[i])))
		ch1[i] <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch1[i])))
		ch1[i] <- paste(" ",gsub("\\("," ( ",ch1[i])," ",sep="")}
	if(expand).i <- covind(envir)
	ex1a <- ex1[!is.na(ex1)]
	if(length(ex1a)>0)for(i in 1:length(ex1a))for(j in 1:length(ch1))
		ch1[j] <- if(expand)gsub(paste(" ",colnames(envir$ccov$ccov)[ex1a[i]],
			" ",sep=""),paste(" ",ndata,"$ccov$ccov[,",
			ex1a[i],"][.i] ",sep=""),ch1[j])
			else gsub(paste(" ",colnames(envir$ccov$ccov)[ex1a[i]],
			" ",sep=""),paste(" ",ndata,"$ccov$ccov[,",
			ex1a[i],"] ",sep=""),ch1[j])
	ex2a <- ex2[!is.na(ex2)]
	if(length(ex2a)>0)for(i in 1:length(ex2a))for(j in 1:length(ch1))
		ch1[j] <- gsub(paste(" ",colnames(envir$tvcov$tvcov)[ex2a[i]],
			" ",sep=""),paste(" ",ndata,"$tvcov$tvcov[,",
			ex2a[i],"] ",sep=""),ch1[j])
	ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch)))
	ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch)))
	ch <- paste(" ",gsub("\\("," ( ",ch)," ",sep="")
	ch2 <- strsplit(ch," ")[[1]]
	ex2a <- un <- un0 <- un1 <- NULL
	if(length(mem2)>0)for(i in 1:length(mem2)){
		ex1a <- NULL
		for(j in 1:length(ch2))
			if(mem2[i]==ch2[j]){
				if(j<length(ch2)&&length(grep("^\\[",ch2[j+1]))>0){
					ex1a <- c(ex1a,paste(ch2[j],ch2[j+1],sep=""))
					un1 <- c(un1,ch2[j])}
				else un0 <- c(un0,ch2[j])}
		if(!is.null(ex1a)){
			ex1a <- unique(ex1a)
			ex2a <- c(ex2a,length(ex1a))
			o <- gsub("(^[[:alnum:]]\\[)|(\\])","",ex1a)
			un <- if(length(grep("[[:alpha:]]",o))>0)c(un,ex1a)
				else c(un,ex1a[order(as.numeric(o))])}}
	if(length(un0)>0){
		un <- if(length(un1)>0&&length(grep(un1,un0))>0)
			c(un,unique(un0[-grep(un1,un0)]))
			else c(un,unique(un0))}
	.fn <- eval(parse(text=paste("function(",paste(mem2,collapse=","),")",paste("eval(attr(.fn,\"model\"))"))))
	attributes(.fn) <- list(model=parse(text=ch1),parameters=un,
		covariates=unique(mem[(!is.na(ex1)|!is.na(ex2))&!fcn]),
		class="formulafn")
	rm(list=ls())
	return(.fn)}

fnenvir.tccov <- function(z, envir=NULL, name=NULL, expand=TRUE){
	if(!is.function(z))return(NULL)
	if(is.null(envir)||(!inherits(envir,"repeated")&&!inherits(envir,"tccov")))stop("envir must be an object of class, repeated or tccov")
	ndata <- if(is.null(name))paste(deparse(substitute(envir)))
	else name
	ch1 <- deparse(z,width=500)
	ch2 <- ch1[1]
	ch1 <- ch1[2:length(ch1)]
	mem2 <- strsplit(gsub("[(),]"," ",ch2)," ")[[1]]
	if(length(mem2)>1)mem2 <- mem2[2:length(mem2)]
	else mem2 <- NULL
	i <- un <- fcn <- ex1 <- ch <- NULL
	for(i in 1:length(ch1))ch <- paste(ch,ch1[i],collapse=" ")
	mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(^[0-9]+)"," ",ch)," ")[[1]]
	if(inherits(envir,"repeated")){
		cn <- colnames(envir$ccov$ccov)
		cc <- "$ccov"}
	else {
		cn <- colnames(envir$ccov)
		cc <- ""}
	if(length(mem)>0){
		ex1 <- match(mem,cn)
		for(i in 1:length(mem)){
			fcn <- c(fcn,if(exists(mem[i])){
				if(mem[i]=="function"||mem[i]=="if"||
					mem[i]=="else"||mem[i]=="for"||
					mem[i]=="while"||mem[i]=="repeat") T
				else is.function(eval(parse(text=mem[i])))}
				else F)}
		un <- unique(mem[is.na(ex1)&!fcn])
		if(length(unique(mem[!is.na(ex1)&!fcn]))==0&&length(un)==0)
			stop("finterp: no variables found")}
	for(i in 1:length(ch1)){
		ch1[i] <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch1[i])))
		ch1[i] <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch1[i])))
		ch1[i] <- paste(" ",gsub("\\("," ( ",ch1[i])," ",sep="")}
	ex1a <- ex1[!is.na(ex1)]
	if(length(ex1a)>0)for(i in 1:length(ex1a))for(j in 1:length(ch1))
		ch1[j] <- gsub(paste(" ",cn[ex1a[i]]," ",sep=""),
			paste(" ",ndata,cc,"$ccov[,",
			ex1a[i],"] ",sep=""),ch1[j])
	ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch)))
	ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch)))
	ch <- paste(" ",gsub("\\("," ( ",ch)," ",sep="")
	ch2 <- strsplit(ch," ")[[1]]
	ex2a <- un <- un0 <- un1 <- NULL
	if(length(mem2)>0)for(i in 1:length(mem2)){
		ex1a <- NULL
		for(j in 1:length(ch2))
			if(mem2[i]==ch2[j]){
				if(j<length(ch2)&&length(grep("^\\[",ch2[j+1]))>0){
					ex1a <- c(ex1a,paste(ch2[j],ch2[j+1],sep=""))
					un1 <- c(un1,ch2[j])}
				else un0 <- c(un0,ch2[j])}
		if(!is.null(ex1a)){
			ex1a <- unique(ex1a)
			ex2a <- c(ex2a,length(ex1a))
			o <- gsub("(^[[:alnum:]]\\[)|(\\])","",ex1a)
			un <- if(length(grep("[[:alpha:]]",o))>0)c(un,ex1a)
				else c(un,ex1a[order(as.numeric(o))])}}
	if(length(un0)>0){
		un <- if(length(un1)>0&&length(grep(un1,un0))>0)
			c(un,unique(un0[-grep(un1,un0)]))
			else c(un,unique(un0))}
	.fn <- eval(parse(text=paste("function(",paste(mem2,collapse=","),")",paste("eval(attr(.fn,\"model\"))"))))
	attributes(.fn) <- list(model=parse(text=ch1),parameters=un,
		covariates=unique(mem[!is.na(ex1)&!fcn]),
		class="formulafn")
	rm(list=ls())
	return(.fn)}

fnenvir.tvcov <- function(z, envir=NULL, name=NULL, expand=TRUE){
	if(!is.function(z))return(NULL)
	if(is.null(envir)||(!inherits(envir,"repeated")&&!inherits(envir,"tvcov")))stop("envir must be an object of class, repeated or tvcov")
	ndata <- if(is.null(name))paste(deparse(substitute(envir)))
	else name
	ch1 <- deparse(z,width=500)
	ch2 <- ch1[1]
	ch1 <- ch1[2:length(ch1)]
	mem2 <- strsplit(gsub("[(),]"," ",ch2)," ")[[1]]
	if(length(mem2)>1)mem2 <- mem2[2:length(mem2)]
	else mem2 <- NULL
	i <- un <- fcn <- ex2 <- ch <- NULL
	for(i in 1:length(ch1))ch <- paste(ch,ch1[i],collapse=" ")
	mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(^[0-9]+)"," ",ch)," ")[[1]]
	if(inherits(envir,"repeated")){
		cn <- colnames(envir$tvcov$tvcov)
		cc <- "$tvcov"}
	else {
		cn <- colnames(envir$tvcov)
		cc <- ""}
	if(length(mem)>0){
		ex2 <- match(mem,cn)
		for(i in 1:length(mem)){
			fcn <- c(fcn,if(exists(mem[i])){
				if(mem[i]=="function"||mem[i]=="if"||
					mem[i]=="else"||mem[i]=="for"||
					mem[i]=="while"||mem[i]=="repeat") T
				else is.function(eval(parse(text=mem[i])))}
				else F)}
		un <- unique(mem[is.na(ex2)&!fcn])
		if(length(unique(mem[!is.na(ex2)&!fcn]))==0&&length(un)==0)
			stop("finterp: no variables found")}
	for(i in 1:length(ch1)){
		ch1[i] <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch1[i])))
		ch1[i] <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch1[i])))
		ch1[i] <- paste(" ",gsub("\\("," ( ",ch1[i])," ",sep="")}
	ex2a <- ex2[!is.na(ex2)]
	if(length(ex2a)>0)for(i in 1:length(ex2a))for(j in 1:length(ch1))
		ch1[j] <- gsub(paste(" ",cn[ex2a[i]]," ",sep=""),
			paste(" ",ndata,cc,"$tvcov[,",
			ex2a[i],"] ",sep=""),ch1[j])
	ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",ch)))
	ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",ch)))
	ch <- paste(" ",gsub("\\("," ( ",ch)," ",sep="")
	ch2 <- strsplit(ch," ")[[1]]
	ex2a <- un <- un0 <- un1 <- NULL
	if(length(mem2)>0)for(i in 1:length(mem2)){
		ex1a <- NULL
		for(j in 1:length(ch2))
			if(mem2[i]==ch2[j]){
				if(j<length(ch2)&&length(grep("^\\[",ch2[j+1]))>0){
					ex1a <- c(ex1a,paste(ch2[j],ch2[j+1],sep=""))
					un1 <- c(un1,ch2[j])}
				else un0 <- c(un0,ch2[j])}
		if(!is.null(ex1a)){
			ex1a <- unique(ex1a)
			ex2a <- c(ex2a,length(ex1a))
			o <- gsub("(^[[:alnum:]]\\[)|(\\])","",ex1a)
			un <- if(length(grep("[[:alpha:]]",o))>0)c(un,ex1a)
				else c(un,ex1a[order(as.numeric(o))])}}
	if(length(un0)>0){
		un <- if(length(un1)>0&&length(grep(un1,un0))>0)
			c(un,unique(un0[-grep(un1,un0)]))
			else c(un,unique(un0))}
	.fn <- eval(parse(text=paste("function(",paste(mem2,collapse=","),")",paste("eval(attr(.fn,\"model\"))"))))
	attributes(.fn) <- list(model=parse(text=ch1),parameters=un,
		covariates=unique(mem[!is.na(ex2)&!fcn]),
		class="formulafn")
	rm(list=ls())
	return(.fn)}

print.formulafn <- function(z){
	if(!is.null(attr(z,"formula"))){
		cat("\nformula:\n")
		print.default(unclass(attr(z,"formula")))}
	if(!is.matrix(attr(z,"model"))){
		model <- deparse(attr(z,"model"))
		model[1] <- sub("expression\\(","",model[1])
		model[length(model)] <- sub("\\)$","",model[length(model)])
		cat("\nmodel function:\n")
		cat(model,sep="\n")}
	if(length(attr(z,"covariates"))>0){
		cat(paste("\ncovariates:\n"))
		for(i in 1:length(attr(z,"covariates")))
			cat(attr(z,"covariates")[i]," ")
		cat("\n")}
	if(length(attr(z,"parameters"))>0){
		cat(paste("\nparameters:\n"))
		for(i in 1:length(attr(z,"parameters")))
			cat(attr(z,"parameters")[i]," ")}
	cat("\n\n")}

model <- function(z, ...) UseMethod("model")

model.formulafn <- function(z) attr(z,"model")

formula.formulafn <- function(z) attr(z,"formula")

covariates.formulafn <- function(z) attr(z,"covariates")

parameters <- function(z, ...) UseMethod("parameters")

parameters.formulafn <- function(z) attr(z,"parameters")
