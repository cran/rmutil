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
#   fmobj(z, envir=parent.frame()))
#   finterp(.z, .envir=parent.frame(),
#	.formula=FALSE, .vector=TRUE, .args=NULL, .start=1, .name=NULL,
#	.expand=TRUE, .intercept=TRUE, .old=NULL){
#   fnenvir(.z, .envir=parent.frame(), .name=NULL, .expand=TRUE)
#
#  DESCRIPTION
#
#    Functions to translate a model formula with unknown parameters
#  into a function and to modify a function to read from a data object.

###
### function to find objects specified in a formula, returning
###indicators of which are parameters, covariates, factors, and functions
###
fmobj <- function(z, envir=parent.frame()){
if(!inherits(z,"formula"))return(NULL)
#
# transform formula to one character string
#
local <- fac <- cov <- ch <- NULL
ch1 <- deparse(z[[length(z)]])
for(j in 1:length(ch1))ch <- paste(ch,ch1[j],collapse=" ",sep="\n")
#
# put extra spaces in character string so that substitutions can be made
#
ch <- gsub("\n"," \n ",gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [ ",ch))))
ch <- gsub("/"," / ",gsub(","," ,",gsub("\\("," ( ",gsub(":"," : ",ch))))
ch <- paste(" ",gsub(" -"," - ",gsub(" \\+"," + ",ch))," ",sep="")
#
# find names of variables and parameters
# removing names of functions that are arguments of other functions
# except times and response
#
mem <- all.vars(z)
fcn <- all.names(z)
fcn <- fcn[!(fcn%in%mem)]
if(length(mem)>0){
	tmp <- vector(length=length(mem))
	for(i in 1:length(mem))
		tmp[i] <- exists(mem[i],envir=envir)&&is.function(eval(parse(text=mem[i]),envir=envir))&&(length(grep(paste(mem[i],","),ch))>0||length(grep(paste(",",mem[i]),ch))>0||(length(grep(paste("=",mem[i]),ch))>0&&length(grep(paste(">=",mem[i]),ch))==0&&length(grep(paste("<=",mem[i]),ch))==0)||length(grep(paste("\\(",mem[i],"\\)"),ch))>0)&&mem[i]!="times"&&mem[i]!="response"
	fcn <- unique(c(fcn,mem[tmp]))
	mem <- mem[!tmp]}
#
# Create indicators for numeric covariates (cov) and for factor
# variables (fac) that exist in the given environment.
# Everything else in mem is an unknown parameter.
#
if(length(mem)>0){
	for(j in 1:length(mem)){
		local <- c(local,length(grep(paste(mem[j],"<-"),ch))>0)
		cov <- c(cov,exists(mem[j],envir=envir)&&
			is.numeric(eval(parse(text=mem[j]),envir=envir)))
		fac <- c(fac,exists(mem[j],envir=envir)&&!cov[j]&&
			is.factor(eval(parse(text=mem[j]),envir=envir)))}
	zz <- list(formula=ch,objects=mem,functions=fcn,
		parameters=!cov&!fac&!local,covariates=cov,factors=fac,
		local=local)}
else zz <- list(formula=ch,objects=mem,functions=fcn)
class(zz) <- "fmobj"
zz}

###
### functions to translate model formulae into functions,
### perhaps transforming them to read from a data object
###
finterp <- function(.z, ...) UseMethod("finterp")

### default method
###
finterp.default <- function(.z, .envir=parent.frame(),
	.formula=FALSE, .vector=TRUE, .args=NULL, .start=1, .name=NULL,
	.expand=TRUE, .intercept=TRUE, .old=NULL, .response=FALSE, ...){
if(!inherits(.z,"formula"))return(NULL)
#
# find the appropriate environment if its name is supplied
#
if(is.name(.envir)){
	if(is.null(.name)).name <- as.character(.envir)
	.envir <- eval(.envir)}
#
# call appropriate special method if it exists
#
if(!is.environment(.envir)){
	if(is.null(.name)).name <- paste(deparse(substitute(.envir)))
	if(inherits(.envir,"repeated"))
		return(finterp.repeated(.z,.envir,.formula,.vector,.args,.start,.name,.expand,.intercept,.old,.response))
	if(inherits(.envir,"tccov"))
		return(finterp.tccov(.z,.envir,.formula,.vector,.args,.start,.name,.expand,.intercept,.old))
	if(inherits(.envir,"tvcov"))
		return(finterp.tvcov(.z,.envir,.formula,.vector,.args,.start,.name,.expand,.intercept,.old))
	if(inherits(.envir,"data.frame"))
		return(finterp.data.frame(.z,.envir,.formula,.vector,.args,.start,.name,.expand,.intercept,.old))}
#
# check for common parameters
#
.pars <- .range <- NULL
if(!is.null(.old)){
	if(!is.list(.old)).old <- list(.old)
	for(.j in .old){
		if(!inherits(.j,"formulafn"))
			stop("objects in .old must have class, formulafn")
		.pars <- c(.pars,attr(.j,"parameters"))
		.range <- c(.range,attr(.j,"range")[1]:attr(.j,"range")[2])}
	if(.start<=max(.range))
		warning("possible conflict in vector indexing - check .start")}
if(!is.null(.args)&&!is.character(.args))
	stop(".args must be a character string")
#
# create indicators for covariates that exist in the given environment (.ex)
# and for factor variables (.fac), and a vector of unique parameter names (.un)
#
.zz <- fmobj(.z)
.ch <- .zz$formula
.mem <- .zz$objects
.fcn <- .zz$functions
if("$"%in%.fcn)stop("sublists not allowed (attach dataframes and use variable names)")
.ex <- .zz$covariates
.fac <- .zz$factors
.local <- .zz$local
rm(.zz)
if(length(.mem)>0){
	.un <- unique(.mem[!.ex&!.fac&!.local])
	if(length(unique(.mem[.ex|.fac|.local]))==0&&length(.un)==0)
		warning("finterp.default: no variables found")}
#
# handle W&R formulae
#
if(length(.mem)==0||all(.ex|.fac|.local)){
	if(.formula)return(.z)
	else {
	# create model matrix and return as attribute of function
		if(any("offset"%in%.fcn))stop("offset not allowed")
		.mt <- terms(.z)
		if(is.numeric(.mt[[2]])){
			.dm <- matrix(1)
			colnames(.dm) <- "(Intercept)"}
		else {
			.dm <- model.matrix(.mt,model.frame(.mt,.envir,na.action=NULL))
			if(!.intercept).dm <- .dm[,-1,drop=FALSE]}
		.fna <- function(.p) as.vector(.dm%*%
			.p[attr(.fna,"range")[1]:attr(.fna,"range")[2]])
		attributes(.fna) <- list(formula=.z,model=colnames(.dm),
			covariates=if(length(.mem)>0)
				unique(.mem[.ex|.fac]) else NULL,
			parameters=paste("p[",1:dim(.dm)[2],"]",sep=""),
			range=c(.start,.start+dim(.dm)[2]-1),
			class="formulafn")
		.obj <- ls(all.names=TRUE)
		rm(list=.obj[.obj!=".fna"&.obj!=".dm"])
		rm(.obj)
		return(.fna)}}
#
# create function for formulae with unknowns making sure there are no
# factor variables present
#
if(!is.null(.fac)&&any(.fac))stop(paste("covariates in formulae with unknowns must not be factors\ncheck",.mem[.fac]))
.fna <- function(.p) eval(attr(.fna,"model"))
if(.vector){
#
# return a function with all parameters as one vector
#
	# check if some parameters are not to be put in the vector
	if(!is.null(.args)){
		.tmp <- match(.args,.un)
		if(all(!is.na(.tmp))).un <- .un[-.tmp]
		.par <- "alist(.p="
		for(.j in 1:length(.args)){
			.par <- paste(.par,",",collapse="")
			.par <- paste(.par,.args[.j],"=",collapse="")}
		.par <- paste(.par,")",collapse="")
		formals(.fna) <- eval(parse(text=.par))}
	# check if some parameters are common to another formula
	if(!is.null(.old)){
		.j <- match(.pars,.un)
		.j <- .j[!is.na(.j)]
		.un <- .un[-.j]
		.pars <- .pars[!is.na(.j)]
		.range <- .range[!is.na(.j)]
		for(.j in 1:length(.pars))
			.ch <- gsub(paste(" ",.pars[.j]," ",sep=""),
				paste(" .p[",.range[.j],"] ",sep=""),.ch)}
	if(length(.un)>0)for(.j in 1:length(.un))
		.ch <- gsub(paste(" ",.un[.j]," ",sep=""),
			paste(" .p[",.start+.j-1,"] ",sep=""),.ch)}
else {
#
# return a function with parameters having their original names
#
	.par <- "alist("
	for(.j in 1:length(.un)){
		if(.j>1).par <- paste(.par,",",collapse="")
		.par <- paste(.par,.un[.j],"=",collapse="")}
# bug reported in 0.63.0
#	if(length(.un)==1).par <- paste(.par,",...=",collapse="")
	.par <- paste(.par,")",collapse="")
	formals(.fna) <- eval(parse(text=.par))}
attributes(.fna) <- list(formula=.z,model=parse(text=.ch),parameters=.un,
	common=.pars,covariates=unique(.mem[.ex]),
	range=c(.start,.start+length(.un)-1),class="formulafn")
.obj <- ls(all.names=TRUE)
rm(list=.obj[.obj!=".fna"])
rm(.obj)
return(.fna)}

### method for repeated objects
###
finterp.repeated <- function(.z, .envir=NULL, .formula=FALSE, .vector=TRUE,
	.args=NULL, .start=1, .name=NULL, .expand=TRUE, .intercept=TRUE,
	.old=NULL, .response=FALSE, ...){
if(!inherits(.z,"formula"))return(NULL)
#
# check for common parameters
#
.pars <- .range <- NULL
if(!is.null(.old)){
	if(!is.list(.old)).old <- list(.old)
	for(.j in .old){
		if(!inherits(.j,"formulafn"))
			stop("objects in .old must have class, formulafn")
		.pars <- c(.pars,attr(.j,"parameters"))
		.range <- c(.range,attr(.j,"range")[1]:attr(.j,"range")[2])}
	if(.start<=max(.range))
		warning("possible conflict in vector indexing - check .start")}
if(!is.null(.args)&&!is.character(.args))
	stop(".args must be a character string")
#
# find the appropriate environment if its name is supplied
#
if(is.name(.envir)){
	if(is.null(.name)).name <- as.character(.envir)
	.envir <- eval(.envir)}
if(is.null(.envir)||!inherits(.envir,"repeated"))stop("envir must be an object of class, repeated")
.ndata <- if(is.null(.name))paste(deparse(substitute(.envir))) else .name
#
# dissect formula
#
.ex1 <- .ex2 <- .ex3 <- .ex4 <- .ex5 <- .ex6 <- NULL
.zz <- fmobj(.z)
.ch <- .zz$formula
.mem <- .zz$objects
.fcn <- .zz$functions
.local <- .zz$local
rm(.zz)
#
# create indicators for variables that exist in the repeated object (.ex) and
# a vector of unique parameter names (.un)
#
if(length(.mem)>0){
	.ex1 <- match(.mem,colnames(.envir$ccov$ccov))
	.ex2 <- match(.mem,colnames(.envir$tvcov$tvcov))
	.ex3 <- match(.mem,"times")
	.ex4 <- match(.mem,"individuals")
	if(any(!is.na(.ex4))&&length(nobs(.envir))==1)
		stop("these are not repeated measurements")
	.ex5 <- match(.mem,"nesting")
	if(any(!is.na(.ex5))&&is.null(.envir$response$nest))
		stop("no nesting variable available")
	.ex6 <- match(.mem,colnames(.envir$response$y))
	if(any(!is.na(.ex2))&&!.expand)
		stop("time-varying covariates present - time-constant ones must be expanded")
	.un <- unique(.mem[is.na(.ex1)&is.na(.ex2)&is.na(.ex3)&is.na(.ex4)&is.na(.ex5)&is.na(.ex6)&!.local])
	if(length(unique(.mem[!is.na(.ex1)|!is.na(.ex2)|!is.na(.ex3)|!is.na(.ex4)|!is.na(.ex5)|!is.na(.ex6)]))==0&&length(.un)==0)
		warning("finterp.repeated: no variables found")}
#
# replace time-constant covariate names by location in data object,
# expanding to the proper length, if required
#
if(.expand).i <- covind(.envir)
.ex1a <- if(is.null(.ex1)) NULL else .ex1[!is.na(.ex1)]
if(length(.ex1a)>0)for(.j in 1:length(.ex1a))
	.ch <- if(.expand)gsub(paste(" ",colnames(.envir$ccov$ccov)[.ex1a[.j]],
		" ",sep=""),paste(" ",.ndata,"$ccov$ccov[,",
		.ex1a[.j],"][.i] ",sep=""),.ch)
		else gsub(paste(" ",colnames(.envir$ccov$ccov)[.ex1a[.j]],
		" ",sep=""),paste(" ",.ndata,"$ccov$ccov[,",
		.ex1a[.j],"] ",sep=""),.ch)
#
# replace time-varying covariate names by location in data object
#
.ex2a <- if(is.null(.ex2)) NULL else .ex2[!is.na(.ex2)]
if(length(.ex2a)>0)for(.j in 1:length(.ex2a))
	.ch <- gsub(paste(" ",colnames(.envir$tvcov$tvcov)[.ex2a[.j]],
		" ",sep=""),paste(" ",.ndata,"$tvcov$tvcov[,",
		.ex2a[.j],"] ",sep=""),.ch)
#
# replace special name, times, by location in data object
#
.ex3a <- if(is.null(.ex3)) NULL else .ex3[!is.na(.ex3)]
if(length(.ex3a)>0)
	.ch <- gsub(" times ",paste(" ",.ndata,"$response$times ",sep=""),.ch)
#
# replace special name, individuals, by location in data object
#
.ex4a <- if(is.null(.ex4)) NULL else .ex4[!is.na(.ex4)]
if(length(.ex4a)>0)
	.ch <- gsub(" individuals ",paste(" as.factor(covind(",.ndata,")) ",sep=""),.ch)
#
# replace special name, nesting, by location in data object
#
.ex5a <- if(is.null(.ex5)) NULL else .ex5[!is.na(.ex5)]
if(length(.ex5a)>0)
	.ch <- gsub(" nesting ",paste(" as.factor(",.ndata,"$response$nest) ",sep=""),.ch)
#
# replace response variable names by location in data object
#
.ex6a <- if(is.null(.ex6)) NULL else .ex6[!is.na(.ex6)]
if(length(.ex6a)>0){
	if(dim(.envir$response$y)[2]==1&&!.response)
		stop(paste(colnames(.envir$response$y)[.ex6a]," is the response and cannot be a covariate"))
	for(.j in 1:length(.ex6a)){
		if(!is.null(.envir$response$n)&&!all(is.na(.envir$response$n[.ex6a[.j]])))
			stop(paste(colnames(.envir$response$y)[.ex6a[.j]]," is binomial and cannot be a covariate"))
		if(!is.null(.envir$response$censor)&&!all(is.na(.envir$response$censor[.ex6a[.j]]))&&!all(.envir$response$censor[.ex6a[.j]]==1))
			stop(paste(colnames(.envir$response$y)[.ex6a[.j]]," is censored and cannot be a covariate"))
		.ch <- gsub(paste(" ",colnames(.envir$response$y)[.ex6a[.j]],
			" ",sep=""),paste(" ",.ndata,"$response$y[,",
			.ex6a[.j],"] ",sep=""),.ch)}}
#
# handle W&R formulae
#
if((is.null(.ex1)&&is.null(.ex2)&&is.null(.ex3)&&is.null(.ex4)&&is.null(.ex5)&&is.null(.ex6))||all(!is.na(.ex1)|!is.na(.ex2)|!is.na(.ex3)|!is.na(.ex4)|!is.na(.ex5)|!is.na(.ex6))){
	if(.formula)return(.z)
	else {
	# create model matrix, change column names, and return as
	# attribute of function
		if(any("offset"%in%.fcn))stop("offset not allowed")
		.ch <- as.formula(paste("~",.ch))
		.mt <- terms(.ch)
		if(is.numeric(.mt[[2]])){
			if(!.intercept)return(NULL)
			.n <- if(.expand||is.null(.envir$ccov$ccov))
					dim(.envir$response$y)[1]
				else dim(.envir$ccov$ccov)[1]
			.dm <- matrix(1)
			colnames(.dm) <- "(Intercept)"
			.fna <- function(.p) rep(.p[attr(.fna,"range")[1]],.n)}
		else {
			.dm <- model.matrix(.mt,model.frame(.mt,na.action=NULL))
			if(!.intercept).dm <- .dm[,-1,drop=FALSE]
			if(length(.ex1a)>0)for(.j in 1:length(.ex1a))
				colnames(.dm) <- gsub(paste(.ndata,"\\$ccov\\$ccov\\[, ",.ex1a[.j],"\\]",sep=""),paste(colnames(.envir$ccov$ccov)[.ex1a[.j]],sep=""),colnames(.dm))
			if(length(.ex2a)>0)for(.j in 1:length(.ex2a))
				colnames(.dm) <- gsub(paste(.ndata,"\\$tvcov\\$tvcov\\[, ",.ex2a[.j],"\\]",sep=""),paste(colnames(.envir$tvcov$tvcov)[.ex2a[.j]],sep=""),colnames(.dm))
			if(length(.ex3a)>0)colnames(.dm) <- gsub(paste(.ndata,"\\$response\\$times",sep=""),"times",colnames(.dm))
			if(length(.ex4a)>0)colnames(.dm) <- gsub(paste("as.factor\\(covind\\(",.ndata,"\\)\\)",sep=""),"individuals",colnames(.dm))
			if(length(.ex5a)>0)colnames(.dm) <- gsub(paste("as.factor\\(",.ndata,"\\$response\\$nest\\)",sep=""),"nesting",colnames(.dm))
			if(length(.ex6a)>0)for(.j in 1:length(.ex6a))
				colnames(.dm) <- gsub(paste(.ndata,"\\$response\\$y\\[, ",.ex6a[.j],"\\]",sep=""),paste(colnames(.envir$response$y)[.ex6a[.j]],sep=""),colnames(.dm))
			.fna <- function(.p) as.vector(.dm%*%
				.p[attr(.fna,"range")[1]:attr(.fna,"range")[2]])}
		attributes(.fna) <- list(formula=.z,model=colnames(.dm),
			covariates=if(length(.mem)>0)
				unique(.mem[(!is.na(.ex1)|!is.na(.ex2)|!is.na(.ex3)|!is.na(.ex4)|!is.na(.ex5)|!is.na(.ex6))]) else NULL,
			parameters=paste("p[",1:dim(.dm)[2],"]",sep=""),
			range=c(.start,.start+dim(.dm)[2]-1),
			class="formulafn")
		.obj <- ls(all.names=TRUE)
		rm(list=.obj[.obj!=".i"&.obj!=".fna"&.obj!=".dm"&.obj!=".n"])
		rm(.obj)
		return(.fna)}}
#
# make sure there are no factor variables present
#
if(length(.ex1a)>0)for(.j in 1:length(.ex1a))if(is.factor(.envir$ccov$ccov[,.ex1a[.j]]))stop(paste(colnames(.envir$ccov$ccov)[.ex1a[.j]],"is a factor variable"))
if(length(.ex2a)>0)for(.j in 1:length(.ex2a))if(is.factor(.envir$tvcov$tvcov[,.ex2a[.j]]))stop(paste(colnames(.envir$tvcov$tvcov)[.ex2a[.j]],"is a factor variable"))
if(length(.ex4a)>0)stop("index for individuals cannot be used in formulae with unknowns")
if(length(.ex5a)>0)stop("index for nesting cannot be used in formulae with unknowns")
#
# create and return the function
#
.fna <- function(.p) eval(attr(.fna,"model"))
#
# check if some parameters are not to be put in the vector
#
if(!is.null(.args)){
	.tmp <- match(.args,.un)
	if(all(!is.na(.tmp))).un <- .un[-.tmp]
	.par <- "alist(.p="
	for(.j in 1:length(.args)){
		.par <- paste(.par,",",collapse="")
		.par <- paste(.par,.args[.j],"=",collapse="")}
	.par <- paste(.par,")",collapse="")
	formals(.fna) <- eval(parse(text=.par))}
#
# check if some parameters are common to another formula
#
if(!is.null(.old)){
	.j <- match(.pars,.un)
	.j <- .j[!is.na(.j)]
	.un <- .un[-.j]
	.pars <- .pars[!is.na(.j)]
	.range <- .range[!is.na(.j)]
	for(.j in 1:length(.pars)).ch <- gsub(paste(" ",.pars[.j]," ",sep=""),
		paste(" .p[",.range[.j],"] ",sep=""),.ch)}
if(length(.un)>0)for(.j in 1:length(.un))
	.ch <- gsub(paste(" ",.un[.j]," ",sep=""),
		paste(" .p[",.start+.j-1,"] ",sep=""),.ch)
attributes(.fna) <- list(formula=.z,model=parse(text=.ch),parameters=.un,
	covariates=unique(.mem[(!is.na(.ex1)|!is.na(.ex2)|!is.na(.ex3)|!is.na(.ex4)|!is.na(.ex5)|!is.na(.ex6))]),
	common=.pars,range=c(.start,.start+length(.un)-1),class="formulafn")
.obj <- ls(all.names=TRUE)
rm(list=.obj[.obj!=".i"&.obj!=".fna"])
rm(.obj)
return(.fna)}

### method for tccov objects
###
finterp.tccov <- function(.z, .envir=NULL, .formula=FALSE, .vector=TRUE,
	.args=NULL, .start=1, .name=NULL, .expand=NULL, .intercept=TRUE,
	.old=NULL, ...){
if(!inherits(.z,"formula"))return(NULL)
#
# check for common parameters
#
.pars <- .range <- NULL
if(!is.null(.old)){
	if(!is.list(.old)).old <- list(.old)
	for(.j in .old){
		if(!inherits(.j,"formulafn"))
			stop("objects in .old must have class, formulafn")
		.pars <- c(.pars,attr(.j,"parameters"))
		.range <- c(.range,attr(.j,"range")[1]:attr(.j,"range")[2])}
	if(.start<=max(.range))
		warning("possible conflict in vector indexing - check .start")}
if(!is.null(.args)&&!is.character(.args))
	stop(".args must be a character string")
#
# find the appropriate environment if its name is supplied
#
if(is.name(.envir)){
	if(is.null(.name)).name <- as.character(.envir)
	.envir <- eval(.envir)}
if(is.null(.envir)||(!inherits(.envir,"repeated")&&!inherits(.envir,"tccov")))
	stop("envir must be an object of class, repeated or tccov")
.ndata <- if(is.null(.name))paste(deparse(substitute(.envir))) else .name
#
# create appropriate names for search in data object
#
if(inherits(.envir,"repeated")){
	.cn <- colnames(.envir$ccov$ccov)
	.cc <- "$ccov"
	.cc2 <- "\\$ccov"}
else {
	.cn <- colnames(.envir$ccov)
	.cc2 <- .cc <- ""}
#
# dissect formula
#
.ex1 <- NULL
.zz <- fmobj(.z)
.ch <- .zz$formula
.mem <- .zz$objects
.fcn <- .zz$functions
.local <- .zz$local
rm(.zz)
#
# create indicator for variables that exist in the tccov object (.ex1)
# and a vector of unique parameter names (.un)
#
if(length(.mem)>0){
	.ex1 <- match(.mem,.cn)
	.un <- unique(.mem[is.na(.ex1)&!.local])
	if(length(unique(.mem[!is.na(.ex1)]))==0&&length(.un)==0)
		warning("finterp.tccov: no variables found")}
#
# replace time-constant covariate names by location in data object,
# expanding to the proper length, if required
#
.ex1a <- if(is.null(.ex1)) NULL else .ex1[!is.na(.ex1)]
if(length(.ex1a)>0)for(.j in 1:length(.ex1a))
	.ch <- gsub(paste(" ",.cn[.ex1a[.j]]," ",sep=""),
		paste(" ",.ndata,.cc,"$ccov[,",
		.ex1a[.j],"] ",sep=""),.ch)
#
# handle W&R formulae
#
if(is.null(.ex1)||all(!is.na(.ex1))){
	if(.formula)return(.z)
	else {
	# create model matrix, change column names, and return as
	# attribute of function
		if(any("offset"%in%.fcn))stop("offset not allowed")
		.ch <- as.formula(paste("~",.ch))
		.mt <- terms(.ch)
		if(is.numeric(.mt[[2]])){
			if(!.intercept)return(NULL)
			.n <- dim(.envir$ccov)[1]
			.dm <- matrix(1)
			colnames(.dm) <- "(Intercept)"
			.fna <- function(.p) rep(.p[attr(.fna,"range")[1]],.n)}
		else {
			.dm <- model.matrix(.mt,model.frame(.mt,na.action=NULL))
			if(!.intercept).dm <- .dm[,-1,drop=FALSE]
			if(length(.ex1a)>0)for(.j in 1:length(.ex1a))
				colnames(.dm) <- gsub(paste(.ndata,.cc2,"\\$ccov\\[, ",.ex1a[.j],"\\]",sep=""),paste(.cn[.ex1a[.j]],sep=""),colnames(.dm))
			.fna <- function(.p) as.vector(.dm%*%
				.p[attr(.fna,"range")[1]:attr(.fna,"range")[2]])}
		attributes(.fna) <- list(formula=.z,model=colnames(.dm),
			covariates=if(length(.mem)>0)
				unique(.mem[!is.na(.ex1)]) else NULL,
			parameters=paste("p[",1:dim(.dm)[2],"]",sep=""),
			range=c(.start,.start+dim(.dm)[2]-1),
			class="formulafn")
		.obj <- ls(all.names=TRUE)
		rm(list=.obj[.obj!=".i"&.obj!=".fna"&.obj!=".dm"&.obj!=".n"])
		rm(.obj)
		return(.fna)}}
#
# create function for formulae with unknowns making sure there are no
# factor variables present
#
if(inherits(.envir,"repeated")){
	if(length(.ex1a)>0)for(.j in 1:length(.ex1a))if(is.factor(.envir$ccov$ccov[,.ex1a[.j]]))stop(paste(colnames(.envir$ccov$ccov)[.ex1a[.j]],"is a factor variable"))}
else {
	if(length(.ex1a)>0)for(.j in 1:length(.ex1a))if(is.factor(.envir$ccov[,.ex1a[.j]]))stop(paste(colnames(.envir$ccov)[.ex1a[.j]],"is a factor variable"))}
#
# create and return the function
#
.fna <- function(.p) eval(attr(.fna,"model"))
#
# check if some parameters are not to be put in the vector
#
if(!is.null(.args)){
	.tmp <- match(.args,.un)
	if(all(!is.na(.tmp))).un <- .un[-.tmp]
	.par <- "alist(.p="
	for(.j in 1:length(.args)){
		.par <- paste(.par,",",collapse="")
		.par <- paste(.par,.args[.j],"=",collapse="")}
	.par <- paste(.par,")",collapse="")
	formals(.fna) <- eval(parse(text=.par))}
#
# check if some parameters are common to another formula
#
if(!is.null(.old)){
	.j <- match(.pars,.un)
	.j <- .j[!is.na(.j)]
	.un <- .un[-.j]
	.pars <- .pars[!is.na(.j)]
	.range <- .range[!is.na(.j)]
	for(.j in 1:length(.pars)).ch <- gsub(paste(" ",.pars[.j]," ",sep=""),
		paste(" .p[",.range[.j],"] ",sep=""),.ch)}
if(length(.un)>0)for(.j in 1:length(.un))
	.ch <- gsub(paste(" ",.un[.j]," ",sep=""),
		paste(" .p[",.start+.j-1,"] ",sep=""),.ch)
attributes(.fna) <- list(formula=.z,model=parse(text=.ch),parameters=.un,
	common=.pars,covariates=unique(.mem[!is.na(.ex1)]),
	range=c(.start,.start+length(.un)-1),class="formulafn")
.obj <- ls(all.names=TRUE)
rm(list=.obj[.obj!=".fna"])
rm(.obj)
return(.fna)}

### method for tvcov objects
###
finterp.tvcov <- function(.z, .envir=NULL, .formula=FALSE, .vector=TRUE,
	.args=NULL, .start=1, .name=NULL, .expand=NULL, .intercept=TRUE,
	.old=NULL, ...){
if(!inherits(.z,"formula"))return(NULL)
#
# check for common parameters
#
.pars <- .range <- NULL
if(!is.null(.old)){
	if(!is.list(.old)).old <- list(.old)
	for(.j in .old){
		if(!inherits(.j,"formulafn"))
			stop("objects in .old must have class, formulafn")
		.pars <- c(.pars,attr(.j,"parameters"))
		.range <- c(.range,attr(.j,"range")[1]:attr(.j,"range")[2])}
	if(.start<=max(.range))
		warning("possible conflict in vector indexing - check .start")}
if(!is.null(.args)&&!is.character(.args))
	stop(".args must be a character string")
#
# find the appropriate environment if its name is supplied
#
if(is.name(.envir)){
	if(is.null(.name)).name <- as.character(.envir)
	.envir <- eval(.envir)}
if(is.null(.envir)||(!inherits(.envir,"repeated")&&!inherits(.envir,"tvcov")))
	stop("envir must be an object of class, repeated or tvcov")
.ndata <- if(is.null(.name))paste(deparse(substitute(.envir))) else .name
#
# create appropriate names for search in data object
#
if(inherits(.envir,"repeated")){
	.cn <- colnames(.envir$tvcov$tvcov)
	.cc <- "$tvcov"
	.cc2 <- "\\$tvcov"}
else {
	.cn <- colnames(.envir$tvcov)
	.cc2 <- .cc <- ""}
#
# dissect formula
#
.ex1 <- NULL
.zz <- fmobj(.z)
.ch <- .zz$formula
.mem <- .zz$objects
.fcn <- .zz$functions
.local <- .zz$local
rm(.zz)
#
# create indicator for variables that exist in the tvcov object (.ex1)
# and a vector of unique parameter names (.un)
#
if(length(.mem)>0){
	.ex1 <- match(.mem,.cn)
	.un <- unique(.mem[is.na(.ex1)&!.local])
	if(length(unique(.mem[!is.na(.ex1)]))==0&&length(.un)==0)
		warning("finterp.tvcov: no variables found")}
#
# replace time-varying covariate names by location in data object
#
.ex1a <- if(is.null(.ex1)) NULL else .ex1[!is.na(.ex1)]
if(length(.ex1a)>0)for(.j in 1:length(.ex1a))
	.ch <- gsub(paste(" ",.cn[.ex1a[.j]]," ",sep=""),
		paste(" ",.ndata,.cc,"$tvcov[,",
		.ex1a[.j],"] ",sep=""),.ch)
#
# handle W&R formulae
#
if(is.null(.ex1)||all(!is.na(.ex1))){
	if(.formula)return(.z)
	else {
	# create model matrix, change column names, and return as
	# attribute of function
		if(any("offset"%in%.fcn))stop("offset not allowed")
		.ch <- as.formula(paste("~",.ch))
		.mt <- terms(.ch)
		if(is.numeric(.mt[[2]])){
			if(!.intercept)return(NULL)
			.n <- dim(.envir$tvcov)[1]
			.dm <- matrix(1)
			colnames(.dm) <- "(Intercept)"
			.fna <- function(.p) rep(.p[attr(.fna,"range")[1]],.n)}
		else {
			.dm <- model.matrix(.mt,model.frame(.mt,na.action=NULL))
			if(!.intercept).dm <- .dm[,-1,drop=FALSE]
			if(length(.ex1a)>0)for(.j in 1:length(.ex1a))
				colnames(.dm) <- gsub(paste(.ndata,.cc2,"\\$ccov\\[, ",.ex1a[.j],"\\]",sep=""),paste(.cn[.ex1a[.j]],sep=""),colnames(.dm))
			.fna <- function(.p) as.vector(.dm%*%
				.p[attr(.fna,"range")[1]:attr(.fna,"range")[2]])}
		attributes(.fna) <- list(formula=.z,model=colnames(.dm),
			covariates=if(length(.mem)>0)
				unique(.mem[!is.na(.ex1)]) else NULL,
			parameters=paste("p[",1:dim(.dm)[2],"]",sep=""),
			range=c(.start,.start+dim(.dm)[2]-1),
			class="formulafn")
		.obj <- ls(all.names=TRUE)
		rm(list=.obj[.obj!=".fna"&.obj!=".dm"&.obj!=".n"])
		rm(.obj)
		return(.fna)}}
#
# create function for formulae with unknowns making sure there are no
# factor variables present
#
if(inherits(.envir,"repeated")){
	if(length(.ex1a)>0)for(.j in 1:length(.ex1a))if(is.factor(.envir$tvcov$tvcov[,.ex1a[.j]]))stop(paste(colnames(.envir$tvcov$tvcov)[.ex1a[.j]],"is a factor variable"))}
else {
	if(length(.ex1a)>0)for(.j in 1:length(.ex1a))if(is.factor(.envir$tvcov[,.ex1a[.j]]))stop(paste(colnames(.envir$tvcov)[.ex1a[.j]],"is a factor variable"))}
#
# create and return the function
#
.fna <- function(.p) eval(attr(.fna,"model"))
#
# check if some parameters are not to be put in the vector
#
if(!is.null(.args)){
	.tmp <- match(.args,.un)
	if(all(!is.na(.tmp))).un <- .un[-.tmp]
	.par <- "alist(.p="
	for(.j in 1:length(.args)){
		.par <- paste(.par,",",collapse="")
		.par <- paste(.par,.args[.j],"=",collapse="")}
	.par <- paste(.par,")",collapse="")
	formals(.fna) <- eval(parse(text=.par))}
#
# check if some parameters are common to another formula
#
if(!is.null(.old)){
	.j <- match(.pars,.un)
	.j <- .j[!is.na(.j)]
	.un <- .un[-.j]
	.pars <- .pars[!is.na(.j)]
	.range <- .range[!is.na(.j)]
	for(.j in 1:length(.pars)).ch <- gsub(paste(" ",.pars[.j]," ",sep=""),
		paste(" .p[",.range[.j],"] ",sep=""),.ch)}
if(length(.un)>0)for(.j in 1:length(.un))
	.ch <- gsub(paste(" ",.un[.j]," ",sep=""),
		paste(" .p[",.start+.j-1,"] ",sep=""),.ch)
attributes(.fna) <- list(formula=.z,model=parse(text=.ch),parameters=.un,
	common=.pars,covariates=unique(.mem[!is.na(.ex1)]),
	range=c(.start,.start+length(.un)-1),class="formulafn")
.obj <- ls(all.names=TRUE)
rm(list=.obj[.obj!=".fna"])
rm(.obj)
return(.fna)}

### method for dataframes
###
finterp.data.frame <- function(.z, .envir=NULL, .formula=FALSE, .vector=TRUE,
	.args=NULL, .start=1, .name=NULL, .expand=NULL, .intercept=TRUE,
	.old=NULL, ...){
if(!inherits(.z,"formula"))return(NULL)
#
# check for common parameters
#
.pars <- .range <- NULL
if(!is.null(.old)){
	if(!is.list(.old)).old <- list(.old)
	for(.j in .old){
		if(!inherits(.j,"formulafn"))
			stop("objects in .old must have class, formulafn")
		.pars <- c(.pars,attr(.j,"parameters"))
		.range <- c(.range,attr(.j,"range")[1]:attr(.j,"range")[2])}
	if(.start<=max(.range))
		warning("possible conflict in vector indexing - check .start")}
if(!is.null(.args)&&!is.character(.args))
	stop(".args must be a character string")
#
# find the appropriate environment if its name is supplied
#
if(is.name(.envir)){
	if(is.null(.name)).name <- as.character(.envir)
	.envir <- eval(.envir)}
.ndata <- if(is.null(.name))paste(deparse(substitute(.envir))) else .name
#
# create appropriate names for search in data object
#
.cn <- colnames(.envir)
#
# dissect formula
#
.ex1 <- NULL
.zz <- fmobj(.z)
.ch <- .zz$formula
.mem <- .zz$objects
.fcn <- .zz$functions
.local <- .zz$local
rm(.zz)
#
# create indicator for variables that exist in the dataframe (.ex1)
# and a vector of unique parameter names (.un)
#
if(length(.mem)>0){
	.ex1 <- match(.mem,.cn)
	.un <- unique(.mem[is.na(.ex1)&!.local])
	if(length(unique(.mem[!is.na(.ex1)]))==0&&length(.un)==0)
		warning("finterp.data.frame: no variables found")}
#
# replace covariate names by location in data object
#
.ex1a <- if(is.null(.ex1)) NULL else .ex1[!is.na(.ex1)]
if(length(.ex1a)>0)for(.j in 1:length(.ex1a))
	.ch <- gsub(paste(" ",.cn[.ex1a[.j]]," ",sep=""),
		paste(" ",.ndata,"$",.cn[.ex1a[.j]],sep=""),.ch)
#
# handle W&R formulae
#
if(is.null(.ex1)||all(!is.na(.ex1))){
	if(.formula)return(.z)
	else {
	# create model matrix, change column names, and return as
	# attribute of function
		if(any("offset"%in%.fcn))stop("offset not allowed")
		.ch <- as.formula(paste("~",.ch))
		.mt <- terms(.ch)
		if(is.numeric(.mt[[2]])){
			if(!.intercept)return(NULL)
			.n <- dim(.envir)[1]
			.dm <- matrix(1)
			colnames(.dm) <- "(Intercept)"
			.fna <- function(.p) rep(.p[attr(.fna,"range")[1]],.n)}
		else {
			.dm <- model.matrix(.mt,model.frame(.mt,data=.envir,na.action=NULL))
			if(!.intercept).dm <- .dm[,-1,drop=FALSE]
			.fna <- function(.p) as.vector(.dm%*%
				.p[attr(.fna,"range")[1]:attr(.fna,"range")[2]])}
		attributes(.fna) <- list(formula=.z,model=colnames(.dm),
			covariates=if(length(.mem)>0)
				unique(.mem[!is.na(.ex1)]) else NULL,
			parameters=paste("p[",1:dim(.dm)[2],"]",sep=""),
			range=c(.start,.start+dim(.dm)[2]-1),
			class="formulafn")
		.obj <- ls(all.names=TRUE)
		rm(list=.obj[.obj!=".i"&.obj!=".fna"&.obj!=".dm"&.obj!=".n"])
		rm(.obj)
		return(.fna)}}
#
# create function for formulae with unknowns making sure there are no
# factor variables present
#
if(length(.ex1a)>0)for(.j in 1:length(.ex1a))if(is.factor(.envir[,.ex1a[.j]]))
	stop(paste(colnames(.envir)[.ex1a[.j]],"is a factor variable"))
#
# create and return the function
#
.fna <- function(.p) eval(attr(.fna,"model"))
#
# check if some parameters are not to be put in the vector
#
if(!is.null(.args)){
	.tmp <- match(.args,.un)
	if(all(!is.na(.tmp))).un <- .un[-.tmp]
	.par <- "alist(.p="
	for(.j in 1:length(.args)){
		.par <- paste(.par,",",collapse="")
		.par <- paste(.par,.args[.j],"=",collapse="")}
	.par <- paste(.par,")",collapse="")
	formals(.fna) <- eval(parse(text=.par))}
#
# check if some parameters are common to another formula
#
if(!is.null(.old)){
	.j <- match(.pars,.un)
	.j <- .j[!is.na(.j)]
	.un <- .un[-.j]
	.pars <- .pars[!is.na(.j)]
	.range <- .range[!is.na(.j)]
	for(.j in 1:length(.pars)).ch <- gsub(paste(" ",.pars[.j]," ",sep=""),
		paste(" .p[",.range[.j],"] ",sep=""),.ch)}
if(length(.un)>0)for(.j in 1:length(.un))
	.ch <- gsub(paste(" ",.un[.j]," ",sep=""),
		paste(" .p[",.start+.j-1,"] ",sep=""),.ch)
attributes(.fna) <- list(formula=.z,model=parse(text=.ch),parameters=.un,
	common=.pars,covariates=unique(.mem[!is.na(.ex1)]),
	range=c(.start,.start+length(.un)-1),class="formulafn")
.obj <- ls(all.names=TRUE)
rm(list=.obj[.obj!=".fna"])
rm(.obj)
return(.fna)}

### functions to find the variables and parameters in functions,
### perhaps transforming them to read from a data object
###
fnenvir <- function(.z, ...) UseMethod("fnenvir")

### default method
###
fnenvir.default <- function(.z, .envir=parent.frame(), .name=NULL,
	.expand=TRUE, .response=FALSE, ...){
if(!is.function(.z))return(NULL)
#
# find the appropriate environment if its name is supplied
#
if(is.name(.envir)){
	if(is.null(.name)).name <- as.character(.envir)
	.envir <- eval(.envir)}
#
# call appropriate special method if it exists
#
if(!is.environment(.envir)){
	if(is.null(.name)).name <- paste(deparse(substitute(.envir)))
	if(inherits(.envir,"repeated"))return(fnenvir.repeated(.z,.envir,.name=.name,.expand,.response))
	if(inherits(.envir,"tccov"))return(fnenvir.tccov(.z,.envir,.name=.name,.expand))
	if(inherits(.envir,"tvcov"))return(fnenvir.tvcov(.z,.envir,.name=.name,.expand))
	if(inherits(.envir,"data.frame"))return(fnenvir.data.frame(.z,.envir,.name=.name,.expand))}
#
# transform function to a character string
#
.ch1 <- deparse(.z,width.cutoff=500)
.ch2 <- .ch1[1]
.ch1 <- .ch1[-1]
#
# find arguments of function
#
.mem2 <- strsplit(gsub("[(),]"," ",.ch2)," ")[[1]]
if(length(.mem2)>0).mem2 <- .mem2[.mem2!=""]
if(length(.mem2)>1).mem2 <- .mem2[2:length(.mem2)]
else .mem2 <- NULL
#
# find body of function
#
.fcn <- .ex <- .ch <- NULL
for(.j in 1:length(.ch1)).ch <- paste(.ch,.ch1[.j],collapse=" ")
#
# remove punctuation from function body
#
#.mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)"," ",.ch)," ")[[1]]
.mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"=\\-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)"," ",.ch)," ")[[1]]
if(length(.mem)>0).mem <- .mem[.mem!=""]
#
# create indicators for variables that exist in the given environment (.ex)
# and for functions (.fcn) and a vector of unique parameter names (.un)
#
if(length(.mem)>0){
	for(.j in 1:length(.mem)){
		.ex <- c(.ex,exists(.mem[.j],envir=.envir))
		.fcn <- c(.fcn,if(exists(.mem[.j])){
			if(.mem[.j]=="function"||.mem[.j]=="if"||
				.mem[.j]=="else"||.mem[.j]=="for"||
				.mem[.j]=="while"||.mem[.j]=="repeat") TRUE
			else is.function(eval(parse(text=.mem[.j])))}
			else FALSE)}
	for(.j in 1:length(.mem)){
		if(!.fcn[.j]&&.ex[.j]&&is.factor(eval(parse(text=.mem[.j]),envir=.envir)))stop(paste(.mem[.j],"is a factor variable"))}
	.un <- unique(.mem[!.ex])
	if(length(unique(.mem[.ex&!.fcn]))==0&&length(.un)==0)
		warning("fnenvir.default: no variables found")}
#
# put extra spaces in character string so that substitutions can be made
#
.ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",.ch)))
.ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",gsub("\\+","+ ",.ch))))
.ch <- paste(" ",gsub("\\("," ( ",gsub(":"," : ",.ch))," ",sep="")
.ch2 <- strsplit(.ch," ")[[1]]
#
# find arguments in body with (.un) and without (.un0) subscripts
#
.un <- .un0 <- .un1 <- NULL
if(length(.mem2)>0)for(.j in 1:length(.mem2)){
	.ex1a <- NULL
	for(.k in 1:length(.ch2))if(.mem2[.j]==.ch2[.k]){
		if(.k<length(.ch2)&&length(grep("^\\[",.ch2[.k+1]))>0){
			.ex1a <- c(.ex1a,paste(.ch2[.k],.ch2[.k+1],sep=""))
			.un1 <- c(.un1,.ch2[.k])}
		else .un0 <- c(.un0,.ch2[.k])}
	if(!is.null(.ex1a)){
		.ex1a <- unique(.ex1a)
		.o <- gsub("(^[[:alnum:]]\\[)|(\\])","",.ex1a)
		.un <- if(length(grep("[[:alpha:]]",.o))>0)c(.un,.ex1a)
			else c(.un,.ex1a[order(as.numeric(.o))])}}
#
# find unique arguments without subscripts and add to .un
#
if(length(.un0)>0){
	if(length(.un1)>0){
		.tmp <- NULL
		for(.k in 1:length(.un1))
			if(length(grep(.un1[.k],.un0))>0)
				.tmp <- c(.tmp,grep(.un1[.k],.un0))
		.un <- c(.un,unique(if(!is.null(.tmp)).un0[-.tmp]else .un0))}
	else .un <- c(.un,unique(.un0))}
#
# create the new function that evaluates its model attribute
#
.fnb <- eval(parse(text=paste("function(",paste(.mem2,collapse=","),")",paste("eval(attr(.fnb,\"model\"))"))))
.ex <- if(length(.fcn)>0&&!is.null(.ex)).ex&!.fcn else NULL
attributes(.fnb) <- list(model=parse(text=.ch1),parameters=.un,
	covariates=unique(.mem[.ex]),class="formulafn")
.obj <- ls(all.names=TRUE)
rm(list=.obj[.obj!=".fnb"])
rm(.obj)
return(.fnb)}

### method for repeated objects
###
fnenvir.repeated <- function(.z, .envir=NULL, .name=NULL, .expand=TRUE,
	.response=FALSE, ...){
if(!is.function(.z))return(NULL)
#
# find the appropriate environment if its name is supplied
#
if(is.name(.envir)){
	if(is.null(.name)).name <- as.character(.envir)
	.envir <- eval(.envir)}
if(is.null(.envir)||!inherits(.envir,"repeated"))stop("envir must be an object of class, repeated")
.ndata <- if(is.null(.name))paste(deparse(substitute(.envir)))
else .name
#
# transform function to a character string
#
.ch1 <- deparse(.z,width.cutoff=500)
.ch2 <- .ch1[1]
.ch1 <- .ch1[-1]
#
# find arguments of function
#
.mem2 <- strsplit(gsub("[(),]"," ",.ch2)," ")[[1]]
if(length(.mem2)>0).mem2 <- .mem2[.mem2!=""]
if(length(.mem2)>1).mem2 <- .mem2[2:length(.mem2)]
else .mem2 <- NULL
#
# find body of function
#
.fcn <- .ex1 <- .ex2 <- .ex3 <- .ex4 <- .ch <- NULL
for(.j in 1:length(.ch1)).ch <- paste(.ch,.ch1[.j],collapse=" ")
#
# remove punctuation from function body
#
#.mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)"," ",.ch)," ")[[1]]
.mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"=\\-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)"," ",.ch)," ")[[1]]
if(length(.mem)>0).mem <- .mem[.mem!=""]
#
# create indicators for variables that exist in the repeated environment (.ex)
# and for functions (.fcn) and a vector of unique parameter names (.un)
#
if(length(.mem)>0){
	.ex1 <- match(.mem,colnames(.envir$ccov$ccov))
	.ex2 <- match(.mem,colnames(.envir$tvcov$tvcov))
	.ex3 <- match(.mem,"times")
	.ex4 <- match(.mem,colnames(.envir$response$y))
	if(any(!is.na(.ex2))&&!.expand)stop("time-varying covariates present - time-constant ones must be expanded")
	for(.j in 1:length(.mem)){
		.fcn <- c(.fcn,if(exists(.mem[.j])){
			if(.mem[.j]=="function"||.mem[.j]=="if"||
				.mem[.j]=="else"||.mem[.j]=="for"||
				.mem[.j]=="while"||.mem[.j]=="repeat") TRUE
			else is.function(eval(parse(text=.mem[.j])))&&is.na(.ex1[.j])&&is.na(.ex2[.j])&&is.na(.ex3[.j])&&is.na(.ex4[.j])}
			else FALSE)}
	.un <- unique(.mem[is.na(.ex1)&is.na(.ex2)&is.na(.ex3)&is.na(.ex4)&!.fcn])
	if(length(unique(.mem[(!is.na(.ex1)|!is.na(.ex2)|!is.na(.ex3)|!is.na(.ex4))&!.fcn]))==0&&length(.un)==0)
		warning("fnenvir.repeated: no variables found")}
#
# put extra spaces in character string so that substitutions can be made
#
for(.j in 1:length(.ch1)){
	.ch1[.j] <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",.ch1[.j])))
	.ch1[.j] <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",.ch1[.j])))
	.ch1[.j] <- paste(" ",gsub("\\("," ( ",.ch1[.j])," ",sep="")}
#
# replace time-constant covariate names by location in data object,
# expanding to the proper length, if required
#
if(.expand).i <- covind(.envir)
.ex1a <- .ex1[!is.na(.ex1)]
if(length(.ex1a)>0)for(.j in 1:length(.ex1a)){
	if(is.factor(.envir$ccov$ccov[,.ex1a[.j]]))stop(paste(colnames(.envir$ccov$ccov)[.ex1a[.j]],"is a factor variable"))
	for(.k in 1:length(.ch1))
	.ch1[.k] <- if(.expand)gsub(paste(" ",colnames(.envir$ccov$ccov)[.ex1a[.j]],
		" ",sep=""),paste(" ",.ndata,"$ccov$ccov[,",
		.ex1a[.j],"][.i] ",sep=""),.ch1[.k])
		else gsub(paste(" ",colnames(.envir$ccov$ccov)[.ex1a[.j]],
		" ",sep=""),paste(" ",.ndata,"$ccov$ccov[,",
		.ex1a[.j],"] ",sep=""),.ch1[.k])}
#
# replace time-varying covariate names by location in data object
#
.ex2a <- .ex2[!is.na(.ex2)]
if(length(.ex2a)>0)for(.j in 1:length(.ex2a)){
	if(is.factor(.envir$tvcov$tvcov[,.ex2a[.j]]))stop(paste(colnames(.envir$tvcov$tvcov)[.ex2a[.j]],"is a factor variable"))
	for(.k in 1:length(.ch1))
		.ch1[.k] <- gsub(paste(" ",colnames(.envir$tvcov$tvcov)[.ex2a[.j]],
			" ",sep=""),paste(" ",.ndata,"$tvcov$tvcov[,",
			.ex2a[.j],"] ",sep=""),.ch1[.k])}
#
# replace special name, times, by location in data object
#
.ex3a <- if(is.null(.ex3)) NULL else .ex3[!is.na(.ex3)]
if(length(.ex3a)>0)for(.k in 1:length(.ch1))
	.ch1[.k] <- gsub(" times ",paste(" ",.ndata,"$response$times ",sep=""),.ch1[.k])
#
# replace response by location in data object
#
.ex4a <- .ex4[!is.na(.ex4)]
if(length(.ex4a)>0){
	if(dim(.envir$response$y)[2]==1&&!.response)
		stop(paste(colnames(.envir$response$y)[.ex4a]," is the response and cannot be a covariate"))
	for(.k in 1:length(.ch1))
		.ch1[.k] <- gsub(paste(" ",colnames(.envir$response$y)[.ex4a[.j]],
			" ",sep=""),paste(" ",.ndata,"$response$y[,",
			.ex4a[.j],"] ",sep=""),.ch1[.k])}
#
# put extra spaces in character string so that substitutions can be made
#
.ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",.ch)))
.ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",gsub("\\+","+ ",.ch))))
.ch <- paste(" ",gsub("\\("," ( ",gsub(":"," : ",.ch))," ",sep="")
.ch2 <- strsplit(.ch," ")[[1]]
#
# find arguments in body with (.un) and without (.un0) subscripts
#
.un <- .un0 <- .un1 <- NULL
if(length(.mem2)>0)for(.j in 1:length(.mem2)){
	.ex1a <- NULL
	for(.k in 1:length(.ch2))if(.mem2[.j]==.ch2[.k]){
		if(.k<length(.ch2)&&length(grep("^\\[",.ch2[.k+1]))>0){
			.ex1a <- c(.ex1a,paste(.ch2[.k],.ch2[.k+1],sep=""))
			.un1 <- c(.un1,.ch2[.k])}
		else .un0 <- c(.un0,.ch2[.k])}
	if(!is.null(.ex1a)){
		.ex1a <- unique(.ex1a)
		.o <- gsub("(^[[:alnum:]]\\[)|(\\])","",.ex1a)
		.un <- if(length(grep("[[:alpha:]]",.o))>0)c(.un,.ex1a)
			else c(.un,.ex1a[order(as.numeric(.o))])}}
#
# find unique arguments without subscripts and add to .un
#
if(length(.un0)>0){
	if(length(.un1)>0){
		.tmp <- NULL
		for(.k in 1:length(.un1))
			if(length(grep(.un1[.k],.un0))>0)
				.tmp <- c(.tmp,grep(.un1[.k],.un0))
		.un <- c(.un,unique(if(!is.null(.tmp)).un0[-.tmp]else .un0))}
	else .un <- c(.un,unique(.un0))}
#
# create the new function that evaluates its model attribute
#
.fnb <- eval(parse(text=paste("function(",paste(.mem2,collapse=","),")",paste("eval(attr(.fnb,\"model\"))"))))
.ex1 <- if(!is.null(.ex1)&&!is.null(.ex2)&&!is.null(.ex3)&&length(.fcn)>0)
	(!is.na(.ex1)|!is.na(.ex2)|!is.na(.ex3))&!.fcn else NULL
attributes(.fnb) <- list(model=parse(text=.ch1),parameters=.un,
	covariates=unique(.mem[.ex1]),class="formulafn")
.obj <- ls(all.names=TRUE)
rm(list=.obj[.obj!=".i"&.obj!=".fnb"])
rm(.obj)
return(.fnb)}

### method for tccov objects
###
fnenvir.tccov <- function(.z, .envir=NULL, .name=NULL, .expand=TRUE, ...){
if(!is.function(.z))return(NULL)
#
# find the appropriate environment if its name is supplied
#
if(is.null(.envir)||(!inherits(.envir,"repeated")&&!inherits(.envir,"tccov")))stop("envir must be an object of class, repeated or tccov")
.ndata <- if(is.null(.name))paste(deparse(substitute(.envir)))
else .name
#
# transform function to a character string
#
.ch1 <- deparse(.z,width.cutoff=500)
.ch2 <- .ch1[1]
.ch1 <- .ch1[-1]
#
# find arguments of function
#
.mem2 <- strsplit(gsub("[(),]"," ",.ch2)," ")[[1]]
if(length(.mem2)>0).mem2 <- .mem2[.mem2!=""]
if(length(.mem2)>1).mem2 <- .mem2[2:length(.mem2)]
else .mem2 <- NULL
#
# find body of function
#
.fcn <- .ex1 <- .ch <- NULL
for(.j in 1:length(.ch1)).ch <- paste(.ch,.ch1[.j],collapse=" ")
#
# remove punctuation from function body
#
#.mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)"," ",.ch)," ")[[1]]
.mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"=\\-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)"," ",.ch)," ")[[1]]
if(length(.mem)>0).mem <- .mem[.mem!=""]
#
# create appropriate names for search in data object
#
if(inherits(.envir,"repeated")){
	.cn <- colnames(.envir$ccov$ccov)
	.cc <- "$ccov"}
else {
	.cn <- colnames(.envir$ccov)
	.cc <- ""}
#
# create indicators for variables that exist in the tccov environment (.ex)
# and for functions (.fcn) and a vector of unique parameter names (.un)
#
if(length(.mem)>0){
	.ex1 <- match(.mem,.cn)
	for(.j in 1:length(.mem)){
		.fcn <- c(.fcn,if(exists(.mem[.j])){
			if(.mem[.j]=="function"||.mem[.j]=="if"||
				.mem[.j]=="else"||.mem[.j]=="for"||
				.mem[.j]=="while"||.mem[.j]=="repeat") TRUE
			else is.function(eval(parse(text=.mem[.j])))&&is.na(.ex1[.j])}
			else FALSE)}
	.un <- unique(.mem[is.na(.ex1)&!.fcn])
	if(length(unique(.mem[!is.na(.ex1)&!.fcn]))==0&&length(.un)==0)
		warning("fnenvir.tccov: no variables found")}
#
# put extra spaces in character string so that substitutions can be made
#
for(.j in 1:length(.ch1)){
	.ch1[.j] <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",.ch1[.j])))
	.ch1[.j] <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",.ch1[.j])))
	.ch1[.j] <- paste(" ",gsub("\\("," ( ",.ch1[.j])," ",sep="")}
#
# replace time-constant covariate names by location in data object,
# expanding to the proper length, if required
#
.ex1a <- .ex1[!is.na(.ex1)]
if(length(.ex1a)>0)for(.j in 1:length(.ex1a)){
	if(inherits(.envir,"repeated")){
		if(is.factor(.envir$ccov$ccov[,.ex1a[.j]]))stop(paste(colnames(.envir$ccov$ccov)[.ex1a[.j]],"is a factor variable"))}
	else if(is.factor(.envir$ccov[,.ex1a[.j]]))stop(paste(colnames(.envir$ccov)[.ex1a[.j]],"is a factor variable"))
	for(.k in 1:length(.ch1)).ch1[.k] <- gsub(paste(" ",.cn[.ex1a[.j]],
		" ",sep=""),paste(" ",.ndata,.cc,"$ccov[,",
		.ex1a[.j],"] ",sep=""),.ch1[.k])}
#
# put extra spaces in character string so that substitutions can be made
#
.ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",.ch)))
.ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",gsub("\\+","+ ",.ch))))
.ch <- paste(" ",gsub("\\("," ( ",gsub(":"," : ",.ch))," ",sep="")
.ch2 <- strsplit(.ch," ")[[1]]
#
# find arguments in body with (.un) and without (.un0) subscripts
#
.un <- .un0 <- .un1 <- NULL
if(length(.mem2)>0)for(.j in 1:length(.mem2)){
	.ex1a <- NULL
	for(.k in 1:length(.ch2))if(.mem2[.j]==.ch2[.k]){
		if(.k<length(.ch2)&&length(grep("^\\[",.ch2[.k+1]))>0){
			.ex1a <- c(.ex1a,paste(.ch2[.k],.ch2[.k+1],sep=""))
			.un1 <- c(.un1,.ch2[.k])}
		else .un0 <- c(.un0,.ch2[.k])}
	if(!is.null(.ex1a)){
		.ex1a <- unique(.ex1a)
		.o <- gsub("(^[[:alnum:]]\\[)|(\\])","",.ex1a)
		.un <- if(length(grep("[[:alpha:]]",.o))>0)c(.un,.ex1a)
			else c(.un,.ex1a[order(as.numeric(.o))])}}
#
# find unique arguments without subscripts and add to .un
#
if(length(.un0)>0){
	if(length(.un1)>0){
		.tmp <- NULL
		for(.k in 1:length(.un1))
			if(length(grep(.un1[.k],.un0))>0)
				.tmp <- c(.tmp,grep(.un1[.k],.un0))
		.un <- c(.un,unique(if(!is.null(.tmp)).un0[-.tmp]else .un0))}
	else .un <- c(.un,unique(.un0))}
#
# create the new function that evaluates its model attribute
#
.fnb <- eval(parse(text=paste("function(",paste(.mem2,collapse=","),")",paste("eval(attr(.fnb,\"model\"))"))))
.ex1 <- if(!is.null(.ex1)&&length(.fcn)>0)!is.na(.ex1)&!.fcn else NULL
attributes(.fnb) <- list(model=parse(text=.ch1),parameters=.un,
	covariates=unique(.mem[.ex1]),class="formulafn")
.obj <- ls(all.names=TRUE)
rm(list=.obj[.obj!=".fnb"])
rm(.obj)
return(.fnb)}

### method for tvcov objects
###
fnenvir.tvcov <- function(.z, .envir=NULL, .name=NULL, .expand=TRUE, ...){
if(!is.function(.z))return(NULL)
#
# find the appropriate environment if its name is supplied
#
if(is.null(.envir)||(!inherits(.envir,"repeated")&&!inherits(.envir,"tvcov")))stop("envir must be an object of class, repeated or tvcov")
.ndata <- if(is.null(.name))paste(deparse(substitute(.envir)))
else .name
#
# transform function to a character string
#
.ch1 <- deparse(.z,width.cutoff=500)
.ch2 <- .ch1[1]
.ch1 <- .ch1[-1]
#
# find arguments of function
#
.mem2 <- strsplit(gsub("[(),]"," ",.ch2)," ")[[1]]
if(length(.mem2)>0).mem2 <- .mem2[.mem2!=""]
if(length(.mem2)>1).mem2 <- .mem2[2:length(.mem2)]
else .mem2 <- NULL
#
# find body of function
#
.fcn <- .ex2 <- .ch <- NULL
for(.j in 1:length(.ch1)).ch <- paste(.ch,.ch1[.j],collapse=" ")
#
# remove punctuation from function body
#
#.mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)"," ",.ch)," ")[[1]]
.mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"=\\-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)"," ",.ch)," ")[[1]]
if(length(.mem)>0).mem <- .mem[.mem!=""]
#
# create appropriate names for search in data object
#
if(inherits(.envir,"repeated")){
	.cn <- colnames(.envir$tvcov$tvcov)
	.cc <- "$tvcov"}
else {
	.cn <- colnames(.envir$tvcov)
	.cc <- ""}
#
# create indicators for variables that exist in the tccov environment (.ex)
# and for functions (.fcn) and a vector of unique parameter names (.un)
#
if(length(.mem)>0){
	.ex2 <- match(.mem,.cn)
	for(.j in 1:length(.mem)){
		.fcn <- c(.fcn,if(exists(.mem[.j])){
			if(.mem[.j]=="function"||.mem[.j]=="if"||
				.mem[.j]=="else"||.mem[.j]=="for"||
				.mem[.j]=="while"||.mem[.j]=="repeat") TRUE
		  #else is.function(eval(parse(text=.mem[.j])))&&is.na(.ex1[.j])} # bruce edit
			else is.function(eval(parse(text=.mem[.j])))&&is.na(.ex2[.j])}
			else FALSE)}
	.un <- unique(.mem[is.na(.ex2)&!.fcn])
	if(length(unique(.mem[!is.na(.ex2)&!.fcn]))==0&&length(.un)==0)
		warning("fnenvir.tvcov: no variables found")}
#
# put extra spaces in character string so that substitutions can be made
#
for(.j in 1:length(.ch1)){
	.ch1[.j] <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",.ch1[.j])))
	.ch1[.j] <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",.ch1[.j])))
	.ch1[.j] <- paste(" ",gsub("\\("," ( ",.ch1[.j])," ",sep="")}
#
# replace time-varying covariate names by location in data object,
# expanding to the proper length, if required
#
.ex2a <- .ex2[!is.na(.ex2)]
if(length(.ex2a)>0)for(.j in 1:length(.ex2a)){
	if(inherits(.envir,"repeated")){
		if(is.factor(.envir$tvcov$tvcov[,.ex2a[.j]]))stop(paste(colnames(.envir$tvcov$tvcov)[.ex2a[.j]],"is a factor variable"))}
	else if(is.factor(.envir$tvcov[,.ex2a[.j]]))stop(paste(colnames(.envir$tvcov)[.ex2a[.j]],"is a factor variable"))
	for(.k in 1:length(.ch1)).ch1[.k] <- gsub(paste(" ",.cn[.ex2a[.j]],
		" ",sep=""),paste(" ",.ndata,.cc,"$tvcov[,",
		.ex2a[.j],"] ",sep=""),.ch1[.k])}
#
# put extra spaces in character string so that substitutions can be made
#
.ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",.ch)))
.ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",gsub("\\+","+ ",.ch))))
.ch <- paste(" ",gsub("\\("," ( ",gsub(":"," : ",.ch))," ",sep="")
.ch2 <- strsplit(.ch," ")[[1]]
#
# find arguments in body with (.un) and without (.un0) subscripts
#
.un <- .un0 <- .un1 <- NULL
if(length(.mem2)>0)for(.j in 1:length(.mem2)){
	.ex1a <- NULL
	for(.k in 1:length(.ch2))if(.mem2[.j]==.ch2[.k]){
		if(.k<length(.ch2)&&length(grep("^\\[",.ch2[.k+1]))>0){
			.ex1a <- c(.ex1a,paste(.ch2[.k],.ch2[.k+1],sep=""))
			.un1 <- c(.un1,.ch2[.k])}
		else .un0 <- c(.un0,.ch2[.k])}
	if(!is.null(.ex1a)){
		.ex1a <- unique(.ex1a)
		.o <- gsub("(^[[:alnum:]]\\[)|(\\])","",.ex1a)
		.un <- if(length(grep("[[:alpha:]]",.o))>0)c(.un,.ex1a)
			else c(.un,.ex1a[order(as.numeric(.o))])}}
#
# find unique arguments without subscripts and add to .un
#
if(length(.un0)>0){
	if(length(.un1)>0){
		.tmp <- NULL
		for(.k in 1:length(.un1))
			if(length(grep(.un1[.k],.un0))>0)
				.tmp <- c(.tmp,grep(.un1[.k],.un0))
		.un <- c(.un,unique(if(!is.null(.tmp)).un0[-.tmp]else .un0))}
	else .un <- c(.un,unique(.un0))}
#
# create the new function that evaluates its model attribute
#
.fnb <- eval(parse(text=paste("function(",paste(.mem2,collapse=","),")",paste("eval(attr(.fnb,\"model\"))"))))
.ex2 <- if(!is.null(.ex2)&&length(.fcn)>0)!is.na(.ex2)&!.fcn else NULL
attributes(.fnb) <- list(model=parse(text=.ch1),parameters=.un,
	covariates=unique(.mem[.ex2]),class="formulafn")
.obj <- ls(all.names=TRUE)
rm(list=.obj[.obj!=".fnb"])
rm(.obj)
return(.fnb)}

### method for dataframes
###
fnenvir.data.frame <- function(.z, .envir=NULL, .name=NULL, .expand=TRUE, ...){
if(!is.function(.z))return(NULL)
#
# find the appropriate environment if its name is supplied
#
.ndata <- if(is.null(.name))paste(deparse(substitute(.envir)))
else .name
#
# transform function to a character string
#
.ch1 <- deparse(.z,width.cutoff=500)
.ch2 <- .ch1[1]
.ch1 <- .ch1[-1]
#
# find arguments of function
#
.mem2 <- strsplit(gsub("[(),]"," ",.ch2)," ")[[1]]
if(length(.mem2)>0).mem2 <- .mem2[.mem2!=""]
if(length(.mem2)>1).mem2 <- .mem2[2:length(.mem2)]
else .mem2 <- NULL
#
# find body of function
#
.fcn <- .ex1 <- .ch <- NULL
for(.j in 1:length(.ch1)).ch <- paste(.ch,.ch1[.j],collapse=" ")
#
# remove punctuation from function body
#
#.mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)"," ",.ch)," ")[[1]]
.mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"=\\-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)"," ",.ch)," ")[[1]]
if(length(.mem)>0).mem <- .mem[.mem!=""]
#
# create appropriate names for search in data object
#
.cn <- colnames(.envir)
#
# create indicators for variables that exist in the dataframe (.ex)
# and for functions (.fcn) and a vector of unique parameter names (.un)
#
if(length(.mem)>0){
	.ex1 <- match(.mem,.cn)
	for(.j in 1:length(.mem)){
		.fcn <- c(.fcn,if(exists(.mem[.j])){
			if(.mem[.j]=="function"||.mem[.j]=="if"||
				.mem[.j]=="else"||.mem[.j]=="for"||
				.mem[.j]=="while"||.mem[.j]=="repeat") TRUE
			else is.function(eval(parse(text=.mem[.j])))&&is.na(.ex1[.j])}
			else FALSE)}
	.un <- unique(.mem[is.na(.ex1)&!.fcn])
	if(length(unique(.mem[!is.na(.ex1)&!.fcn]))==0&&length(.un)==0)
		warning("fnenvir.data.frame: no variables found")}
#
# put extra spaces in character string so that substitutions can be made
#
for(.j in 1:length(.ch1)){
	.ch1[.j] <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",.ch1[.j])))
	.ch1[.j] <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",.ch1[.j])))
	.ch1[.j] <- paste(" ",gsub("\\("," ( ",.ch1[.j])," ",sep="")}
#
# replace covariate names by location in data object
#
.ex1a <- .ex1[!is.na(.ex1)]
if(length(.ex1a)>0)for(.j in 1:length(.ex1a)){
	if(is.factor(.envir[,.ex1a[.j]]))stop(paste(colnames(.envir)[.ex1a[.j]],"is a factor variable"))
	for(.k in 1:length(.ch1)).ch1[.k] <- gsub(paste(" ",.cn[.ex1a[.j]],
		" ",sep=""),paste(" ",.ndata,"$",.cn[.ex1a[.j]],
		sep=""),.ch1[.k])}
#
# put extra spaces in character string so that substitutions can be made
#
.ch <- gsub("\\^"," ^ ",gsub("\\)"," )",gsub("\\["," [",.ch)))
.ch <- gsub("-","- ",gsub("/"," / ",gsub(","," ,",gsub("\\+","+ ",.ch))))
.ch <- paste(" ",gsub("\\("," ( ",gsub(":"," : ",.ch))," ",sep="")
.ch2 <- strsplit(.ch," ")[[1]]
#
# find arguments in body with (.un) and without (.un0) subscripts
#
.un <- .un0 <- .un1 <- NULL
if(length(.mem2)>0)for(.j in 1:length(.mem2)){
	.ex1a <- NULL
	for(.k in 1:length(.ch2))if(.mem2[.j]==.ch2[.k]){
		if(.k<length(.ch2)&&length(grep("^\\[",.ch2[.k+1]))>0){
			.ex1a <- c(.ex1a,paste(.ch2[.k],.ch2[.k+1],sep=""))
			.un1 <- c(.un1,.ch2[.k])}
		else .un0 <- c(.un0,.ch2[.k])}
	if(!is.null(.ex1a)){
		.ex1a <- unique(.ex1a)
		.o <- gsub("(^[[:alnum:]]\\[)|(\\])","",.ex1a)
		.un <- if(length(grep("[[:alpha:]]",.o))>0)c(.un,.ex1a)
			else c(.un,.ex1a[order(as.numeric(.o))])}}
#
# find unique arguments without subscripts and add to .un
#
if(length(.un0)>0){
	if(length(.un1)>0){
		.tmp <- NULL
		for(.k in 1:length(.un1))
			if(length(grep(.un1[.k],.un0))>0)
				.tmp <- c(.tmp,grep(.un1[.k],.un0))
		.un <- c(.un,unique(if(!is.null(.tmp)).un0[-.tmp]else .un0))}
	else .un <- c(.un,unique(.un0))}
#
# create the new function that evaluates its model attribute
#
.fnb <- eval(parse(text=paste("function(",paste(.mem2,collapse=","),")",paste("eval(attr(.fnb,\"model\"))"))))
.ex1 <- if(!is.null(.ex1)&&length(.fcn)>0)!is.na(.ex1)&!.fcn else NULL
attributes(.fnb) <- list(model=parse(text=.ch1),parameters=.un,
	covariates=unique(.mem[.ex1]),class="formulafn")
.obj <- ls(all.names=TRUE)
rm(list=.obj[.obj!=".fnb"])
rm(.obj)
return(.fnb)}

### print methods
###
print.formulafn <- function(x, ...){
  z <- x; rm(x)
if(!is.null(attr(z,"formula"))){
	cat("\nformula:\n")
	print.default(unclass(attr(z,"formula")))}
if(!is.character(attr(z,"model"))){
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
	if(length(attr(z,"common"))>0)cat(paste("\nnew parameters:\n"))
	else cat(paste("\nparameters:\n"))
	for(i in 1:length(attr(z,"parameters")))
		cat(attr(z,"parameters")[i]," ")
	cat("\n")}
if(length(attr(z,"common"))>0){
	cat(paste("\ncommon parameters:\n"))
	for(i in 1:length(attr(z,"common")))
		cat(attr(z,"common")[i]," ")
	cat("\n")}
cat("\n")}

print.fmobj <- function(x, ...){
  z <- x; rm(x)
if(any(z$parameters)){
	tmp <- unique(z$objects[z$parameters])
	cat(paste("\nparameters:\n"))
	for(i in 1:length(tmp))
		cat(tmp[i]," ")
	cat("\n")}
if(any(z$covariates)){
	tmp <- unique(z$objects[z$covariates])
	cat(paste("\ncovariates:\n"))
	for(i in 1:length(tmp))
		cat(tmp[i]," ")
	cat("\n")}
if(any(z$factors)){
	tmp <- unique(z$objects[z$factors])
	cat(paste("\nfactors:\n"))
	for(i in 1:length(tmp))
		cat(tmp[i]," ")
	cat("\n")}
if(!is.null(z$functions)){
	cat(paste("\nfunctions:\n"))
	for(i in 1:length(z$functions))
		cat(z$functions[i]," ")
	cat("\n")}
cat("\n")
}

### miscellaneous methods
###
### extract the model
###
model <- function(z, ...) UseMethod("model")

model.formulafn <- function(z, ...) attr(z,"model")

### extract the original formula
###
formula.formulafn <- function(x, ...) attr(x,"formula")

### extract the covariate names
###
covariates.formulafn <- function(z, ...) attr(z,"covariates")

### extract the parameter names
###
parameters <- function(z, ...) UseMethod("parameters")

parameters.formulafn <- function(z, ...) attr(z,"parameters")
