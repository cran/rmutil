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
#  DESCRIPTION
#
#    Utility functions for printing repeated measurements results

residuals.gnlr <- function(z) z$residuals
fitted.values.gnlr <- function(z) z$fitted.values
coefficients.gnlr <- function(z) z$coefficients
weights.gnlr <- function(z) z$prior.weights
df.residual.gnlr <- function(z) z$df
deviance.gnlr <- function(z) 2*z$maxlike

print.gnlr <- function(z) {
	sht <- z$nps>0
	mix <- z$npm>0
	gnlmm <- !is.null(z$points)
	censor <- if(mix)!is.null(z$censor) else z$censor
	npl <- z$npl-gnlmm
	np1 <- z$npl+1
	np1a <- z$npl+z$npm+1
	np2 <- z$npl+z$npm+z$nps
	np3 <- np2+1
	np <- z$npl+z$npm+z$nps+z$npf
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	if(mix&&censor)cat(z$censor,"")
	if(censor)cat("censored ")
	if(z$npf>0&(z$dist=="inverse Gauss"|z$dist=="logistic"|z$dist=="gamma"|z$dist=="Weibull"|z$dist=="extreme value"))
		cat("generalized ")
	cat(z$dist,"distribution\n\n")
	if(gnlmm){
		cat(" with normal mixing distribution on",z$scale,"scale\n")
		cat(" (",z$points," point Gauss-Hermite integration)\n\n",sep="")}
	t <- deparse(z$likefn)
	cat("Log likelihood function:",t[2:length(t)],"",sep="\n")
	if(z$npl>0){
		cat("Location function:\n")
		if(!is.null(attr(z$mu,"formula")))cat(deparse(attr(z$mu,"formula")),sep="\n")
		else if(!is.null(attr(z$mu,"model"))){
			t <- deparse(attr(z$mu,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		if(inherits(z$linear[[1]],"formulafn"))
			cat("Linear part: ",deparse(attr(z$linear[[1]],"formula")),sep="\n")}
	if(mix){
		cat("\nMixture function:\n")
		if(!is.null(attr(z$mix,"formula")))cat(deparse(attr(z$mix,"formula")),sep="\n")
		else if(!is.null(attr(z$mix,"model"))){
			t <- deparse(attr(z$mix,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		if(inherits(z$linear[[2]],"formulafn"))
			cat("Linear part: ",deparse(attr(z$linear[[2]],"formula")),sep="\n")}
	if(sht){
		if(!mix){
			cat("\nLog shape function:\n")
			if(!is.null(attr(z$shape,"formula")))cat(deparse(attr(z$shape,"formula")),sep="\n")
			else if(!is.null(attr(z$shape,"model"))){
				t <- deparse(attr(z$shape,"model"))
				t[1] <- sub("expression\\(","",t[1])
				t[length(t)] <- sub("\\)$","",t[length(t)])
				cat(t,sep="\n")}
			if(inherits(z$linear[[2]],"formulafn"))
				cat("Linear part: ",deparse(attr(z$linear[[2]],"formula")),sep="\n")}
		if(!is.null(z$family)){
			cat("\n(Log) family function:\n")
			if(!is.null(attr(z$family,"formula")))cat(deparse(attr(z$family,"formula")),sep="\n")
			else if(!is.null(attr(z$family,"model"))){
				t <- deparse(attr(z$family,"model"))
				t[1] <- sub("expression\\(","",t[1])
				t[length(t)] <- sub("\\)$","",t[length(t)])
				cat(t,sep="\n")}
			if(inherits(z$linear[[3]],"formulafn"))
				cat("Linear part: ",deparse(attr(z$linear[[3]],"formula")),sep="\n")}}
	cat("\n-Log likelihood   ",z$maxlike,"\n")
	cat("Degrees of freedom",z$df,"\n")
	cat("AIC               ",z$aic,"\n")
	cat("Iterations        ",z$iterations,"\n\n")
	if(npl>0){
		cat("Location parameters:\n")
		cname <- if(is.matrix(attr(z$mu,"model")))colnames(attr(z$mu,"model"))
			else if(length(grep("linear",attr(z$mu,"parameters")))>0)
			attr(z$mu,"parameters")[grep("\\[",attr(z$mu,"parameters"))]
			else attr(z$mu,"parameters")
		if(!is.null(z$linear[[1]])&&!is.null(attr(z$linear[[1]],"parameters")))cname <- c(colnames(attr(z$linear[[1]],"model")),cname)
		coef.table <- cbind(z$coefficients[1:npl],z$se[1:npl])
		dimnames(coef.table) <- list(cname, c("estimate", "se"))
		print.default(coef.table,digits=4,print.gap=2)}
	if(z$npm>0){
		cat("\nMixture parameters:\n")
		cname <- if(is.matrix(attr(z$mix,"model")))colnames(attr(z$mix,"model"))
			else if(length(grep("linear",attr(z$mix,"parameters")))>0)
			attr(z$mix,"parameters")[grep("\\[",attr(z$mix,"parameters"))]
			else attr(z$mix,"parameters")
		if(!is.null(z$linear[[2]])&&!is.null(attr(z$linear[[2]],"parameters")))cname <- c(colnames(attr(z$linear[[2]],"model")),cname)
		coef.table <- cbind(z$coefficients[np1:(np-sht)],z$se[np1:(np-sht)])
		dimnames(coef.table) <- list(cname, c("estimate", "se"))
		print.default(coef.table,digits=4,print.gap=2)}
	if(gnlmm){
		cat("\nMixing standard deviation:\n")
		coef.table <- cbind(z$coefficients[z$npl],z$se[z$npl])
		dimnames(coef.table) <- list(" ", c("estimate", "se"))
		print.default(coef.table, digits=4, print.gap=2)}
	if(z$nps>0){
		cat("\nShape parameters:\n")
		cname <- if(z$npm>0)" "
			else if(is.matrix(attr(z$shape,"model")))colnames(attr(z$shape,"model"))
			else if(length(grep("linear",attr(z$shape,"parameters")))>0||length(grep("mu",attr(z$shape,"parameters")))>0)
			attr(z$shape,"parameters")[grep("\\[",attr(z$shape,"parameters"))]
			else attr(z$shape,"parameters")
		if(!is.null(z$linear[[2]])&&!is.null(attr(z$linear[[2]],"parameters")))cname <- c(colnames(attr(z$linear[[2]],"model")),cname)
		coef.table <- cbind(z$coefficients[np1a:np2],z$se[np1a:np2])
		dimnames(coef.table) <- list(cname, c("estimate", "se"))
		print.default(coef.table,digits=4,print.gap=2)}
	if(z$npf>0){
		cat("\nFamily parameters:\n")
		cname <- if(is.matrix(attr(z$family,"model")))colnames(attr(z$family,"model"))
			else if(length(grep("linear",attr(z$family,"parameters")))>0)
			attr(z$family,"parameters")[grep("\\[",attr(z$family,"parameters"))]
			else attr(z$family,"parameters")
		if(!is.null(z$linear[[3]])&&!is.null(attr(z$linear[[3]],"parameters")))cname <- c(colnames(attr(z$linear[[3]],"model")),cname)
		coef.table <- cbind(z$coefficients[np3:np],z$se[np3:np])
		dimnames(coef.table) <- list(cname, c("estimate", "se"))
		print.default(coef.table, digits=4, print.gap=2)}
	if(np>1){
		cat("\nCorrelations:\n")
		dimnames(z$corr) <- list(seq(1,np),seq(1,np))
		print.default(z$corr, digits=4)}
	invisible(z)}
