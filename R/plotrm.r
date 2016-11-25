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
#     plot(mprofile(z, ...))
#     plot(iprofile(z, ...))
#     plot.residuals(z, ...)
#
#  DESCRIPTION
#
#    Utility functions for plotting repeated measurements profiles
# and residuals

### functions to create information for plotting marginal profiles
###   from dynamic models
mprofile <- function(z, ...) UseMethod("mprofile")

mprofile.default <- function(z, times=NULL, mu=NULL, ...){
if(is.null(times)&&is.null(z$response$times)&&!inherits(z,"kalsurv"))
	stop("No times available")
if(is.null(mu)||!is.function(mu)){
#
# use fitted values from model
#
	if(is.null(z$pred))stop("Fitted values not available")
	# if response is transformed, transform predicted values
	if(!is.null(z$transform)){
		if(z$transform=="exp")z$pred <- log(z$pred)
		else if(z$transform=="square")z$pred  <- sqrt(z$pred)
		else if(z$transform=="sqrt")z$pred <- z$pred^2
		else if(z$transform=="log")z$pred <- exp(z$pred)}
	# store times corresponding to predicted values
	if(inherits(z,"kalsurv")){
		z$ptimes <- vector("double",length(z$response$y))
		tmp1 <- 1
		for(i in 1:length(nobs(z))){
			tmp2 <- sum(nobs(z)[1:i])
			z$ptimes[tmp1:tmp2] <- cumsum(z$response$y[tmp1:tmp2])
			tmp1 <- tmp2+1}
		z$response$times <- z$ptimes}
	else z$ptimes <- z$response$times}
else {
#
# if a function is supplied, calculate predicted values
#
	z$ptimes <- if(is.null(times))seq(min(z$response$times),
		max(z$response$times),length.out=25) else times
	z$pred <- mu(z$coef,z$ptimes)}
class(z) <- c("mprofile",class(z))
invisible(z)}

### function to plot marginal profiles
###
plot.mprofile <- function(x, nind=1, intensity=FALSE, add=FALSE,
	ylim=range(z$pred,na.rm=TRUE), lty=NULL, ylab=NULL, xlab=NULL, ...){
  z <- x; rm(x)
if(max(nind)>length(nobs(z)))stop("no such individual")
if(inherits(z,"kalsurv")){
#
# event history data may have ties
#
	for(i in 1:length(z$response$y))if(z$response$y[i]==0)
		z$response$y[i] <- z$response$y[i-1]
	if(is.null(xlab))xlab <- "Chronological time"
	if(intensity){
	# modify predicted values, if intensity wanted
		z$pred <- 1/z$pred
		if(is.null(ylab))ylab <- "Mean intensity"}
	else if(is.null(ylab))ylab <- "Time between events"}
else {
	if(is.null(xlab))xlab <- "Time"
	if(is.null(ylab))ylab <- "Fitted value"}
if(length(z$ptimes)==length(z$response$times)){
#
# use fitted values from model
#
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
	# plot separately for each individual chosen
		if(is.null(z$response$nest)) kk <- nest <- 1
		else {
			nest <- unique(z$response$nest)
			kk <- z$response$nest}
		for(k in nest){
		# if nesting, start lines over in each cluster
			j <- j+1
			if(is.null(lty)) lt <- (lt%%6)+1
			else lt <- lty[j]
			if(first){
				plot(z$ptimes[ii==i&kk==k],z$pred[ii==i&kk==k],type="l",ylim=ylim,lty=lt,xlab=xlab,ylab=ylab,...)
				first <- FALSE}
			else lines(z$ptimes[ii==i&kk==k],z$pred[ii==i&kk==k],lty=lt)}}}
else {
#
# a function was supplied to calculate predicted values
#
	if(is.null(ylim))ylim <- c(min(z$pred),max(z$pred))
	if(is.null(lty))lty <- 1
	if(!add)plot(z$ptimes, z$pred, type="l", ylim=ylim, lty=lty,
		xlab=xlab, ylab=ylab, ...)
	else lines(z$ptimes, z$pred, lty=lty)}
if(!is.null(z$pse)){
#
# plot standard errors if available (carma only)
#
	lines(z$ptimes, z$pse[,1], lty=3)
	lines(z$ptimes, z$pse[,2], lty=3)}}

### functions to create information for plotting individual profiles
### from dynamic models
###
iprofile <- function(z, ...) UseMethod("iprofile")

iprofile.default <- function(z, plotsd=FALSE, ...){
if(is.null(z$response$times)&&!inherits(z,"kalsurv"))stop("No times available")
if(!inherits(z,"recursive"))
	stop("The object must have class, recursive")
else if(is.null(z$rpred))stop("Individual profiles not available")
#
# if response is transformed, transform predicted values
# and standard deviations
#
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
if(inherits(z,"kalsurv")){
	z$response$times <- vector("double",length(z$response$y))
	tmp1 <- 1
	for(i in 1:length(nobs(z))){
		tmp2 <- sum(nobs(z)[1:i])
		z$response$times[tmp1:tmp2] <- cumsum(z$response$y[tmp1:tmp2])
		tmp1 <- tmp2+1}}
class(z) <- c("iprofile",class(z))
invisible(z)}

### function to plot individual profiles
###
plot.iprofile <- function(x, nind=1, observed=TRUE, intensity=FALSE, add=FALSE,
	lty=NULL, pch=NULL, ylab=NULL, xlab=NULL, main=NULL, ylim=NULL,
	xlim=NULL, ...){
  z <- x; rm(x)
if(max(nind)>length(nobs(z)))stop("no such individual")
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
	ylim <- c(min(z$rpred-3*z$sdr,na.rm=TRUE),max(z$rpred+3*z$sdr,na.rm=TRUE))
ns <- length(nind)
ii <- covind(z$response)
pc <- -1
lt <- 0
first <- !add
#
# set up plotting controls
#
if(is.null(main))
	main <- ifelse(ns==1,paste("Individual ",nind),"")
if(is.null(ylim))ylim <- c(min(c(z$rpred,z$response$y),na.rm=TRUE),
	max(c(z$rpred,z$response$y),na.rm=TRUE))
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
#
# plot separately for each individual chosen
#
if(is.null(z$response$nest)) kk <- nest <- 1
else {
	nest <- unique(z$response$nest)
	kk <- z$response$nest}
for(i in nind){
	for(k in nest){
	# if nesting, start lines over in each cluster
		j <- j+1
		if(is.null(pch))pc <- (pc+1)%%6
		else pc <- pch[j]
		if(is.null(lty)) lt <- (lt%%6)+1
		else lt <- lty[j]
		if(first){
			plot(z$resp$times[ii==i&kk==k&na],z$rpred[ii==i&kk==k&na],type="l",lty=lt,ylab=ylab,xlab=xlab,main=main,ylim=ylim,xlim=xlim,...)
			first <- FALSE}
		else lines(z$resp$times[ii==i&kk==k&na],z$rpred[ii==i&kk==k&na],lty=lt)
		# if required, plot data points
		if(observed)points(z$resp$times[ii==i&kk==k],z$resp$y[ii==i&kk==k],pch=pc)
		# if required, plot standard deviations
		if(!is.null(z$psd)){
			lines(z$resp$times[ii==i&kk==k&na],z$psd[ii==i&kk==k&na,1],lty=3)
			lines(z$resp$times[ii==i&kk==k&na],z$psd[ii==i&kk==k&na,2],lty=3)}}}}

### function to plot residuals
###
plot.residuals <- function(x, X=NULL, subset=NULL, ccov=NULL,
	nind=NULL, recursive=TRUE, pch=20, ylab="Residual",
	xlab=NULL, main=NULL, ...){
  z <- x; rm(x)
  x <- X; rm(X)
na <- TRUE
#
# check if model produced by a function in one of my libraries
#
resp <- !is.null(z$response$y)
if(!resp){
	nind <- ccov <- NULL
	recursive <- FALSE}
if(is.character(x))x <- match.arg(x,c("response","fitted"))
#
# find residuals to be plotted
#
if(resp){
	n <- length(z$response$y)
	res <- if(inherits(z,"recursive"))
		residuals(z, recursive=recursive) else residuals(z)
	# handle ties in event histories
	if(inherits(z,"kalsurv"))for(i in 1:length(z$response$y))
		if(z$response$y[i]==0)z$response$y[i] <- z$response$y[i-1]}
else {
	res <- residuals(z)
	n <- length(res)}
#
# indicator of subset to be plotted
#
if(is.null(subset))subset <- rep(TRUE,n)
else {
	tmp <- rep(FALSE,n)
	tmp[subset] <- TRUE
	subset <- tmp}
#
# find individuals and subset to be plotted
#
if(resp&&!is.null(nind)){
	if(is.null(subset))subset <- rep(FALSE,n)
	for(i in nind)subset <- subset&covind(z)==i}
#
# find values for x axis
#
if(is.null(x)){
	# plot by times if possible
	x <- z$response$times
	if(is.null(x))stop("x must be specified")
	if(is.null(xlab))xlab <- "Time"}
else if(is.numeric(x)){
	# plot against vector supplied if possible
	if(length(x)!=n)
		stop("x variable must have same length as residuals")
	if(is.null(xlab))xlab <- paste(deparse(substitute(x)))}
else if(x=="response"){
	# plot against response if possible
	x <- if(resp)z$response$y else z$y
	if(is.null(x))stop("response variable not found")
	if(is.null(xlab))xlab <- "Response"}
else if(x=="fitted"){
	# plot against fitted values if possible
	x <- if(resp)fitted(z, recursive=recursive) else fitted(z)
	na <- !is.na(x)
	if(is.null(xlab))xlab <- "Fitted values"}
#
# plot
#
if(is.null(ccov))
	plot(x[subset&na], res[subset&na], pch=pch, ylab=ylab,
		xlab=xlab, main=main, ...)
else if(length(ccov)>1)stop("Only one covariate can be specified")
else {
	mat <- match(ccov,colnames(z$ccov$ccov))
	if(is.na(mat))stop("covariate not found")
	un <- unique(z$ccov$ccov[,mat])
	tmp <- par()$mfg[3:4]
	ss <- ceiling(sqrt(length(un)))
	if(length(un)==ss*(ss-1))ss1 <- ss-1
	else ss1 <- ss
	par(mfrow=c(ss,ss1))
	for(i in un){
		ind <- 1:sum(nobs(z))
		ind[rep(z$ccov$ccov[,mat],nobs(z))!=i] <- NA
		main <- paste("Covariate ",ccov,"=",i)
		plot(x[subset&ind&na],res[subset&ind&na],
			pch=pch,ylab=ylab,xlab=xlab,main=main,...)}
	par(mfrow=tmp)}}
