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
#     plot.profile(z, ...)
#     plot.iprofile(z, ...)
#     plot.residuals(z, ...)
#
#  DESCRIPTION
#
#    Utility functions for plotting repeated measurements profiles
# and residuals

plot.profile <- function(z, ...)
	UseMethod("plot.profile")

plot.profile.default <- function(z, times=NULL, nind=1, mu=NULL, add=F,
	ylim=NULL, lty=NULL, ylab="Fitted value", xlab="Time", ...){
	if(max(nind)>length(z$response$nobs))stop("no such individual")
	if(is.null(mu)||!is.function(mu)){
		if(is.null(z$pred))stop("Fitted values not available")
		ns <- length(nind)
		ii <- covind(z$response)
		if(!is.null(lty)){
			if(length(lty)!=1&&length(lty)!=ns)
				stop("lty must have a value for each observation")
			else if(length(lty)==1)lty <- rep(lty,ns)}
		first <- !add
		j <- 0
		lt <- 0
		for(i in nind){
			j <- j+1
			if(is.null(lty)) lt <- (lt%%4)+1
			else lt <- lty[j]
			if(first){
				plot(z$response$times[ii==i], z$pred[ii==i], type="l", ylim=ylim, lty=lt, xlab=xlab, ylab=ylab, ...)
				first <- FALSE}
			else lines(z$response$times[ii==i], z$pred[ii==i], lty=lt)}}
	else {
		yy <- mu(z$coef,times)
		if(is.null(ylim))ylim <- c(min(yy),max(yy))
		if(is.null(lty))lty <- 1
		if(is.null(times))times <- seq(min(z$response$times),
			max(z$response$times),length.out=25)
		if(!add)plot(times, yy, type="l", ylim=ylim, lty=lty,
			xlab=xlab, ylab=ylab, ...)
		else lines(times, yy, lty=lty)}}

plot.iprofile <- function(z, ...)
	UseMethod("plot.iprofile")

plot.iprofile.default <- function(z, nind=1, obs=T, add=F, plotsd=F, lty=NULL,
	pch=NULL, ylab="Recursive fitted value", xlab="Time", main=NULL,
	ylim=NULL, xlim=NULL, ...){
	if(!inherits(z,"recursive"))
		stop("The object must have class, recursive")
	else if(is.null(z$rpred))stop("Individual profiles not available")
	if(max(nind)>length(z$response$nobs))stop("no such individual")
	ns <- length(nind)
	ii <- covind(z$response)
	pc <- -1
	lt <- 0
	first <- !add
	if(is.null(main))
		main <- ifelse(ns==1,paste("Individual ",nind),"")
	if(is.null(ylim))ylim <- c(min(c(z$rpred,z$response$y),na.rm=T),
		max(c(z$rpred,z$response$y),na.rm=T))
	if(is.null(xlim))xlim <- c(min(z$resp$times),max(z$resp$times))
	if(!is.null(lty)){
		if(length(lty)!=1&&length(lty)!=ns)
			stop("lty must have a value for each observation")
		else if(length(lty)==1)lty <- rep(lty,ns)}
	if(!is.null(pch)){
		if(length(pch)!=1&&length(pch)!=ns)
			stop("pch must have a value for each observation")
		else if(length(pch)==1)pch <- rep(pch,ns)}
	na <- !is.na(z$rpred)
	j <- 0
	for(i in nind){
	        j <- j+1
		if(is.null(pch))pc <- (pc+1)%%4
		else pc <- pch[j]
		if(is.null(lty)) lt <- (lt%%4)+1
		else lt <- lty[j]
		if(first){
			plot(z$resp$times[ii==i&na],z$rpred[ii==i&na],type="l",
				lty=lt, ylab=ylab, xlab=xlab, main=main, 
				ylim=ylim, xlim=xlim, ...)
			first <- F}
		else lines(z$resp$times[ii==i&na],z$rpred[ii==i&na],lty=lt)
		if(obs)points(z$resp$times[ii==i],z$resp$y[ii==i],pch=pc)
		if(plotsd){
			lines(z$resp$times[ii==i&na],z$rpred[ii==i&na]+
				2*z$sdr[ii==i],lty=3)
			lines(z$resp$times[ii==i&na],z$rpred[ii==i&na]-
				2*z$sdr[ii==i],lty=3)}}}

plot.residuals <- function(z, ...)
	UseMethod("plot.residuals")

plot.residuals.default <- function(z, x=NULL, subset=NULL, ccov=NULL,
	nind=NULL, recursive=TRUE, pch=20, ylab="Residual",
	xlab=NULL, main=NULL, ...){
	if(!inherits(z,"recursive"))
		stop("The object must have class, recursive")
	na <- TRUE
	if(!is.null(subset)){
		tmp <- rep(FALSE,length(z$response$y))
		tmp[subset] <- TRUE
		subset <- tmp}
	if(!is.null(nind)){
		if(is.null(subset))subset <- rep(FALSE,length(z$response$y))
		for(i in nind)subset <- subset|covind(z)==i}
	if(is.null(subset))subset <- rep(TRUE,length(z$response$y))
	res <- residuals(z, recursive=recursive)
	if(is.character(x))x <- match.arg(x,c("response","fitted"))
	if(is.null(x)){
		x <- z$response$times
		if(is.null(xlab))xlab <- "Time"}
	else if(is.numeric(x)){
		if(length(x)!=length(z$response$y))
			stop("x variable must have same length as residuals")
		if(is.null(xlab))xlab <- paste(deparse(substitute(x)))}
	else if(x=="response"){
		x <- z$response$y
		if(is.null(xlab))xlab <- "Response"}
	else if(x=="fitted"){
		x <- fitted(z, recursive=recursive)
		na <- !is.na(x)
		if(is.null(xlab))xlab <- "Fitted values"}
	if(is.null(ccov))
		plot(x[subset&na], res[subset&na], pch=pch, ylab=ylab,
			xlab=xlab, main=main, ...)
	else if(length(ccov)>1)stop("Only one covariate can be specified")
	else {
		un <- unique(z$response$ccov[,ccov])
		tmp <- par()$mfg[3:4]
		ss <- ceiling(sqrt(length(un)))
		if(length(un)==ss*(ss-1))ss1 <- ss-1
		else ss1 <- ss
		par(mfrow=c(ss,ss1))
		for(i in un){
			ind <- (1:sum(z$response$nobs))[rep(z$response$ccov[,ccov],z$resp$nobs)==i]
			main <- paste("Covariate ",cov,"=",i)
			plot(x[subset&ind&na],res[subset&ind&na], pch=pch,
				ylab=ylab, xlab=xlab, main=main, ...)}
		par(mfrow=tmp)}}
