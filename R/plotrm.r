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
#     plot(profile(z, ...))
#     plot(iprofile(z, ...))
#     plot(residuals(z, ...))
#
#  DESCRIPTION
#
#    Utility functions for plotting repeated measurements profiles
# and residuals

profile <- function(z, ...) UseMethod("profile")

profile.default <- function(z, times=NULL, mu=NULL){
	if(is.null(mu)||!is.function(mu)){
		if(is.null(z$pred))stop("Fitted values not available")
		if(!is.null(z$transform)){
			if(z$transform=="exp")z$pred <- log(z$pred)
			else if(z$transform=="square")z$pred  <- sqrt(z$pred)
			else if(z$transform=="sqrt")z$pred <- z$pred^2
			else if(z$transform=="log")z$pred <- exp(z$pred)}
		z$ptimes <- z$response$times}
	else {
		z$ptimes <- if(is.null(times))seq(min(z$response$times),
			max(z$response$times),length.out=25) else times
		z$pred <- mu(z$coef,z$ptimes)}
	class(z) <- c("profile",class(z))
	invisible(z)}

plot.profile <- function(z, nind=1, intensity=F, add=F,
	ylim=range(z$pred,na.rm=T), lty=NULL, ylab=NULL, xlab=NULL, ...){
	if(max(nind)>length(z$response$nobs))stop("no such individual")
	if(inherits(z,"kalsurv")){
		for(i in 1:length(z$response$y))if(z$response$y[i]==0)
			z$response$y[i] <- z$response$y[i-1]
		if(is.null(xlab))xlab <- "Chronological time"
		if(intensity){
			z$pred <- 1/z$pred
			if(is.null(ylab))ylab <- "Mean intensity"}
		else if(is.null(ylab))ylab <- "Time between events"}
	else {
		if(is.null(xlab))xlab <- "Time"
		if(is.null(ylab))ylab <- "Fitted value"}
	if(length(z$ptimes)==length(z$response$times)){
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
			if(is.null(z$response$nest)) kk <- nest <- 1
			else {
				nest <- unique(z$response$nest)
				kk <- z$response$nest}
			for(k in nest){
				j <- j+1
				if(is.null(lty)) lt <- (lt%%4)+1
				else lt <- lty[j]
				if(first){
					plot(z$ptimes[ii==i&kk==k],z$pred[ii==i&kk==k],type="l",ylim=ylim,lty=lt,xlab=xlab,ylab=ylab,...)
					first <- FALSE}
				else lines(z$ptimes[ii==i&kk==k],z$pred[ii==i&kk==k],lty=lt)}}}
	else {
		if(is.null(ylim))ylim <- c(min(z$pred),max(z$pred))
		if(is.null(lty))lty <- 1
		if(!add)plot(z$ptimes, z$pred, type="l", ylim=ylim, lty=lty,
			xlab=xlab, ylab=ylab, ...)
		else lines(z$ptimes, z$pred, lty=lty)}
	if(!is.null(z$pse)){
		lines(z$ptimes, z$pse[,1], lty=3)
		lines(z$ptimes, z$pse[,2], lty=3)}}

iprofile <- function(z, ...) UseMethod("iprofile")

iprofile.default <- function(z, plotsd=F){
	if(!inherits(z,"recursive"))
		stop("The object must have class, recursive")
	else if(is.null(z$rpred))stop("Individual profiles not available")
	if(!is.null(z$transform)){
                if(z$transform=="exp"){
                	if(plotsd){
                		sd1 <- log(z$rpred+2*z$sdr)
                		sd2 <- log(z$rpred-2*z$sdr)}
                	z$rpred <- log(z$rpred)}
                else if(z$transform=="square"){
                	if(plotsd){
                		sd1 <- sqrt(z$rpred+2*z$sdr)
                		sd2 <- sqrt(z$rpred-2*z$sdr)}
                	z$rpred  <- sqrt(z$rpred)}
                else if(z$transform=="sqrt"){
                	if(plotsd){
                		sd1 <- (z$rpred+2*z$sdr)^2
                		sd2 <- (z$rpred-2*z$sdr)^2}
                	z$rpred <- z$rpred^2}
                else if(z$transform=="log"){
                	if(plotsd){
                		sd1 <- exp(z$rpred+2*z$sdr)
                		sd2 <- exp(z$rpred-2*z$sdr)}
                	z$rpred <- exp(z$rpred)}
		else {
			sd1 <- z$rpred+2*z$sdr
			sd2 <- z$rpred-2*z$sdr}
		if(plotsd)z$psd <- cbind(sd1,sd2)}
	class(z) <- c("iprofile",class(z))
	invisible(z)}

plot.iprofile <- function(z, nind=1, observed=T, intensity=F, add=F, lty=NULL,
	pch=NULL, ylab=NULL, xlab=NULL, main=NULL, ylim=NULL, xlim=NULL, ...){
	if(max(nind)>length(z$response$nobs))stop("no such individual")
	if(inherits(z,"kalsurv")){
		for(i in 1:length(z$response$y))if(z$response$y[i]==0)
			z$response$y[i] <- z$response$y[i-1]
		if(is.null(xlab))xlab <- "Chronological time"
		if(intensity){
			z$rpred <- 1/z$rpred
			z$response$y <- 1/z$response$y
			if(is.null(ylab))ylab <- "Mean intensity"}
		else if(is.null(ylab))ylab <- "Time between events"}
	else {
		if(is.null(xlab))xlab <- "Time"
		if(is.null(ylab))ylab <- "Recursive fitted value"}
	if(is.null(ylim)&&!is.null(z$sdr)&&z$transform=="identity")
		ylim <- c(min(z$rpred-3*z$sdr,na.rm=T),max(z$rpred+3*z$sdr,na.rm=T))
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
		if(is.null(z$response$nest)) kk <- nest <- 1
		else {
			nest <- unique(z$response$nest)
			kk <- z$response$nest}
		for(k in nest){
			j <- j+1
			if(is.null(pch))pc <- (pc+1)%%4
			else pc <- pch[j]
			if(is.null(lty)) lt <- (lt%%4)+1
			else lt <- lty[j]
			if(first){
				plot(z$resp$times[ii==i&kk==k&na],z$rpred[ii==i&kk==k&na],type="l",lty=lt,ylab=ylab,xlab=xlab,main=main,ylim=ylim,xlim=xlim,...)
				first <- F}
			else {
				if(k==1)lines(z$resp$times[ii==i&kk==k&na],z$rpred[ii==i&kk==k&na],lty=lt)
				else lines(z$resp$times[ii==i&kk==k&na],c(z$pred[ii==i&kk==k&na][1],z$rpred[ii==i&kk==k&na][-1]),lty=lt)}
			if(observed)points(z$resp$times[ii==i&kk==k],z$resp$y[ii==i&kk==k],pch=pc)
			if(!is.null(z$psd)){
				lines(z$resp$times[ii==i&kk==k&na],z$psd[ii==i&kk==k&na,1],lty=3)
				lines(z$resp$times[ii==i&kk==k&na],z$psd[ii==i&kk==k&na,2],lty=3)}}}}

plot.residuals <- function(z, x=NULL, subset=NULL, ccov=NULL,
	nind=NULL, recursive=TRUE, pch=20, ylab="Residual",
	xlab=NULL, main=NULL, ...){
	na <- TRUE
	reps <- !is.null(z$response$y)
	if(!reps){
		nind <- ccov <- NULL
		recursive <- FALSE}
	if(is.character(x))x <- match.arg(x,c("response","fitted"))
	if(reps){
		n <- length(z$response$y)
		res <- if(inherits(z,"recursive"))
			residuals(z, recursive=recursive) else residuals(z)
		if(inherits(z,"kalsurv"))for(i in 1:length(z$response$y))
			if(z$response$y[i]==0)
				z$response$y[i] <- z$response$y[i-1]}
	else {
		res <- residuals(z)
		n <- length(res)}
	if(is.null(subset))subset <- rep(TRUE,n)
	else {
		tmp <- rep(FALSE,n)
		tmp[subset] <- TRUE
		subset <- tmp}
	if(reps&&!is.null(nind)){
		if(is.null(subset))subset <- rep(FALSE,n)
		for(i in nind)subset <- subset|covind(z)==i}
	if(is.null(x)){
		x <- z$response$times
		if(is.null(x))stop("x must be specified")
		if(is.null(xlab))xlab <- "Time"}
	else if(is.numeric(x)){
		if(length(x)!=n)
			stop("x variable must have same length as residuals")
		if(is.null(xlab))xlab <- paste(deparse(substitute(x)))}
	else if(x=="response"){
		x <- if(reps)z$response$y else z$y
		if(is.null(x))stop("response variable not found")
		if(is.null(xlab))xlab <- "Response"}
	else if(x=="fitted"){
		x <- if(reps) fitted(z, recursive=recursive) else fitted(z)
		na <- !is.na(x)
		if(is.null(xlab))xlab <- "Fitted values"}
	if(is.null(ccov))
		plot(x[subset&na], res[subset&na], pch=pch, ylab=ylab,
			xlab=xlab, main=main, ...)
	else if(length(ccov)>1)stop("Only one covariate can be specified")
	else {
		mat <- match(ccov,colnames(z$response$ccov))
		if(is.na(mat))stop("covariate not found")
		un <- unique(z$response$ccov[,mat])
		tmp <- par()$mfg[3:4]
		ss <- ceiling(sqrt(length(un)))
		if(length(un)==ss*(ss-1))ss1 <- ss-1
		else ss1 <- ss
		par(mfrow=c(ss,ss1))
		for(i in un){
			ind <- (1:sum(z$response$nobs))[rep(z$response$ccov[,ccov],z$resp$nobs)==i]
			main <- paste("Covariate ",cov,"=",i)
			plot(x[subset&ind&na],res[subset&ind&na],
				pch=pch,ylab=ylab,xlab=xlab,main=main,...)}
		par(mfrow=tmp)}}
