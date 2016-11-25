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
#  DESCRIPTION
#
#    Utility functions for printing nonlinear model results

### standard methods for gnlm models
###

weights.gnlm <- function(object, ...) object$prior.weights

df.residual.gnlm <- function(object, ...) object$df

deviance.gnlm <- function(object, ...) 2*object$maxlike

coef.gnlm <- function(object, ...){
  z <- object; rm(object)
gnlmm <- !is.null(z$points)
mix <- z$npm>0||!is.null(z$mix)
npl <- z$npl-gnlmm
gnlmix <- z$npm>0&&is.null(z$mix) ## bruce added this line
if(npl>0){
	if(z$distribution=="own"){
		cname <- if(is.character(attr(z$likefn,"model")))
			attr(z$likefn,"model")
		else if(length(grep("linear",attr(z$likefn,"parameters")))>0)
		attr(z$likefn,"parameters")[grep("\\[",attr(z$likefn,"parameters"))]
		else attr(z$likefn,"parameters")}
	else {
		cname <- if(is.character(attr(z$mu,"model")))attr(z$mu,"model")
		else if(length(grep("linear",attr(z$mu,"parameters")))>0)
		attr(z$mu,"parameters")[grep("\\[",attr(z$mu,"parameters"))]
		else attr(z$mu,"parameters")}
	if(!is.null(z$linmodel[[1]]))cname <-
	c(cname,z$linmodel[[1]])}
if(gnlmm)cname <- c(cname," ")
if(z$npm>0&&!is.null(z$mix)){
	cname <- c(cname,if(is.character(attr(z$mix,"model")))
		attr(z$mix,"model")
		else if(length(grep("linear",attr(z$mix,"parameters")))>0)
		attr(z$mix,"parameters")[grep("\\[",attr(z$mix,"parameters"))]
		else attr(z$mix,"parameters"))
	if(!is.null(z$linmodel[[2]]))cname <- c(cname,z$linmodel[[2]])}
if(z$common||z$nps>0){
	if(!is.null(z$shape))cname <- c(cname,if(mix&&!gnlmix)" "
		else if(is.character(attr(z$shape,"model")))
			attr(z$shape,"model")
		else if(length(grep("linear",attr(z$shape,"parameters")))>0||
		length(grep("mu",attr(z$shape,"parameters")))>0)
		attr(z$shape,"parameters")[grep("\\[",attr(z$shape,"parameters"))]
		else attr(z$shape,"parameters"))
	if(!is.null(z$linmodel[[2]]))cname <- c(cname,z$linmodel[[2]])}
if(z$npf>0||!is.null(z$family)){
	cname <- c(cname,if(is.character(attr(z$family,"model")))
		attr(z$family,"model")
		else if(length(grep("linear",attr(z$family,"parameters")))>0)
		attr(z$family,"parameters")[grep("\\[",attr(z$family,"parameters"))]
		else attr(z$family,"parameters"))
	if(!is.null(z$linmodel[[3]]))cname <- c(cname,z$linmodel[[3]])}
if(z$common)cname <- unique(cname)
if(mix)cname <- c(cname," ")
coef <- z$coef
names(coef) <- cname
coef}

vcov.gnlm <- function(object, ...){
  z <- object; rm(object)
gnlmm <- !is.null(z$points)
mix <- z$npm>0||!is.null(z$mix)
npl <- z$npl-gnlmm
gnlmix <- z$npm>0&&is.null(z$mix) ## bruce added this line
if(npl>0){
	if(z$distribution=="own"){
		cname <- if(is.character(attr(z$likefn,"model")))
			attr(z$likefn,"model")
		else if(length(grep("linear",attr(z$likefn,"parameters")))>0)
		attr(z$likefn,"parameters")[grep("\\[",attr(z$likefn,"parameters"))]
		else attr(z$likefn,"parameters")}
	else {
		cname <- if(is.character(attr(z$mu,"model")))attr(z$mu,"model")
		else if(length(grep("linear",attr(z$mu,"parameters")))>0)
		attr(z$mu,"parameters")[grep("\\[",attr(z$mu,"parameters"))]
		else attr(z$mu,"parameters")}
	if(!is.null(z$linmodel[[1]]))cname <-
	c(cname,z$linmodel[[1]])}
if(gnlmm)cname <- c(cname," ")
if(z$npm>0&&!is.null(z$mix)){
	cname <- c(cname,if(is.character(attr(z$mix,"model")))
		attr(z$mix,"model")
		else if(length(grep("linear",attr(z$mix,"parameters")))>0)
		attr(z$mix,"parameters")[grep("\\[",attr(z$mix,"parameters"))]
		else attr(z$mix,"parameters"))
	if(!is.null(z$linmodel[[2]]))cname <- c(cname,z$linmodel[[2]])}
if(z$common||z$nps>0){
	if(!is.null(z$shape))cname <- c(cname,if(mix&&!gnlmix)" "
		else if(is.character(attr(z$shape,"model")))
			attr(z$shape,"model")
		else if(length(grep("linear",attr(z$shape,"parameters")))>0||
		length(grep("mu",attr(z$shape,"parameters")))>0)
		attr(z$shape,"parameters")[grep("\\[",attr(z$shape,"parameters"))]
		else attr(z$shape,"parameters"))
	if(!is.null(z$linmodel[[2]]))cname <- c(cname,z$linmodel[[2]])}
if(z$npf>0||!is.null(z$family)){
	cname <- c(cname,if(is.character(attr(z$family,"model")))
		attr(z$family,"model")
		else if(length(grep("linear",attr(z$family,"parameters")))>0)
		attr(z$family,"parameters")[grep("\\[",attr(z$family,"parameters"))]
		else attr(z$family,"parameters"))
	if(!is.null(z$linmodel[[3]]))cname <- c(cname,z$linmodel[[3]])}
if(z$common)cname <- unique(cname)
if(mix)cname <- c(cname," ")
vcov <- z$cov
rownames(vcov) <- cname
colnames(vcov) <- cname
vcov}

### print function for bnlr, gnlr, gnlr3, gnlmm, gnlmm3, and fmr
###
print.gnlm <- function(x,digits=max(4,.Options$digits-3),correlation=TRUE, ...) {
  z <- x; rm(x)
sht <- z$nps>0||!is.null(z$shape)
mix <- z$npm>0||!is.null(z$mix)
gnlmm <- !is.null(z$points)
gnlmix <- z$npm>0&&is.null(z$mix)
censor <- if(mix&&!gnlmix)!is.null(z$censor) else z$censor
npl <- z$npl-gnlmm
np1 <- z$npl+1
np1a <- z$npl+z$npm*(!gnlmix)+1
np2 <- z$npl+z$npm*(!gnlmix)+z$nps
np3 <- np2+1
np <- z$npl+z$npm+z$nps+z$npf
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
if(mix&&censor&&!is.logical(z$censor))cat(z$censor,"")
if(censor)cat("censored ")
if(z$npf>0&&(z$dist=="inverse Gauss"||z$dist=="logistic"||z$dist=="gamma"||
	z$dist=="Weibull"||z$dist=="extreme value"))cat("generalized ")
if(!is.null(z$dist))cat(z$dist,"distribution\n\n")
else if(!is.null(z$link))cat(z$link,"link\n\n")
if(!is.null(z$mixture))cat(z$mixture,"mixing distribution\n\n")
if(gnlmm){
	cat(" with normal mixing distribution on",z$scale,"scale\n")
	cat(" (",z$points," point Gauss-Hermite integration)\n\n",sep="")}
cat("Response:",z$respname,"\n\n")
cat("Log likelihood function:\n")
if(z$distribution=="own"){
	if(!is.null(attr(z$likefn,"formula")))
		cat(deparse(attr(z$likefn,"formula")),sep="\n")
	else if(!is.null(attr(z$likefn,"model"))){
		t <- deparse(attr(z$likefn,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}}
else {
	t <- deparse(z$likefn)
	cat(t[2:length(t)],"",sep="\n")}
if((z$npl>0||!is.null(z$mu))&&z$distribution!="own"){
	cat("Location function:\n")
	if(!is.null(attr(z$mu,"formula")))
		cat(deparse(attr(z$mu,"formula")),sep="\n")
	else if(!is.null(attr(z$mu,"model"))){
		t <- deparse(attr(z$mu,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	if(!is.null(z$linear[[1]])){
		cat("Linear part:\n")
		print(z$linear[[1]])}}
if(mix){
	if(!is.null(z$mix))cat("\nMixture function:\n")
	if(!is.null(attr(z$mix,"formula")))
		cat(deparse(attr(z$mix,"formula")),sep="\n")
	else if(!is.null(attr(z$mix,"model"))){
		t <- deparse(attr(z$mix,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	if(!is.null(z$linear[[2]])){
		cat("Linear part:\n")
		print(z$linear[[2]])}}
if(sht){
	if(!mix){
		cat("\nLog shape function:\n")
		if(!is.null(attr(z$shape,"formula")))
			cat(deparse(attr(z$shape,"formula")),sep="\n")
		else if(!is.null(attr(z$shape,"model"))){
			t <- deparse(attr(z$shape,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		if(!is.null(z$linear[[2]])){
			cat("Linear part:\n")
			print(z$linear[[2]])}}
	if(!is.null(z$family)){
		cat("\n(Log) family function:\n")
		if(!is.null(attr(z$family,"formula")))
			cat(deparse(attr(z$family,"formula")),sep="\n")
		else if(!is.null(attr(z$family,"model"))){
			t <- deparse(attr(z$family,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		if(!is.null(z$linear[[3]])){
			cat("Linear part:\n")
			print(z$linear[[3]])}}}
cat("\n-Log likelihood   ",z$maxlike,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n\n")
if(npl>0){
	if(z$common)cat("Common parameters:\n")
	else if(z$distribution=="own")cat("Model parameters:\n")
	else cat("Location parameters:\n")
	if(z$distribution=="own"){
		cname <- if(is.character(attr(z$likefn,"model")))
			attr(z$likefn,"model")
		else if(length(grep("linear",attr(z$likefn,"parameters")))>0)
		attr(z$likefn,"parameters")[grep("\\[",attr(z$likefn,"parameters"))]
		else attr(z$likefn,"parameters")}
	else {
		cname <- if(is.character(attr(z$mu,"model")))attr(z$mu,"model")
		else if(length(grep("linear",attr(z$mu,"parameters")))>0)
		attr(z$mu,"parameters")[grep("\\[",attr(z$mu,"parameters"))]
		else attr(z$mu,"parameters")}
	if(!is.null(z$linmodel[[1]]))cname <- c(cname,z$linmodel[[1]])
	coef.table <- cbind(z$coefficients[1:npl],z$se[1:npl])
	if(!z$common){
		dimnames(coef.table) <- list(cname, c("estimate", "se"))
		print.default(coef.table,digits=digits,print.gap=2)
		cname <- coef.table <- NULL}}
if(z$npm>0&&!is.null(z$mix)){
	if(!z$common)cat("\nMixture parameters:\n")
	cname <- c(cname,if(is.character(attr(z$mix,"model")))
		attr(z$mix,"model")
		else if(length(grep("linear",attr(z$mix,"parameters")))>0)
		attr(z$mix,"parameters")[grep("\\[",attr(z$mix,"parameters"))]
		else attr(z$mix,"parameters"))
	if(!is.null(z$linmodel[[2]]))cname <- c(cname,z$linmodel[[2]])
	if(!z$common)coef.table <- cbind(z$coefficients[np1:(np-sht)],
		z$se[np1:(np-sht)])
	dimnames(coef.table) <- list(cname, c("estimate", "se"))
	print.default(coef.table,digits=digits,print.gap=2)
	cname <- coef.table <- NULL}
if(z$common||z$nps>0){
	if(!is.null(z$shape))cname <- c(cname,if(mix&&!gnlmix)" "
		else if(is.character(attr(z$shape,"model")))
			attr(z$shape,"model")
		else if(length(grep("linear",attr(z$shape,"parameters")))>0||
		length(grep("mu",attr(z$shape,"parameters")))>0)
		attr(z$shape,"parameters")[grep("\\[",attr(z$shape,"parameters"))]
		else attr(z$shape,"parameters"))
	if(!is.null(z$linmodel[[2]]))cname <- c(cname,z$linmodel[[2]])
	if(!z$common)coef.table <- cbind(z$coefficients[np1a:np2],
		z$se[np1a:np2])
	if(z$common&&is.null(z$family)){
		dimnames(coef.table) <- list(unique(cname),c("estimate","se"))
		print.default(coef.table,digits=digits,print.gap=2)}
	if(is.null(z$shape)&&z$nps==1){
		coef.table <- cbind(z$coefficients[np2],z$se[np2])
		cname <- " "}}
if(gnlmm){
	cat("\nMixing standard deviation:\n")
	coef.table2 <- cbind(z$coefficients[z$npl],z$se[z$npl])
	dimnames(coef.table2) <- list(" ", c("estimate", "se"))
	print.default(coef.table2,digits=digits,print.gap=2)}
if(z$nps>0&&(!z$common||mix)){
	cat("\nShape parameters:\n")
	dimnames(coef.table) <- list(cname,c("estimate","se"))
	print.default(coef.table,digits=digits,print.gap=2)
	cname <- coef.table <- NULL}
if(z$npf>0||!is.null(z$family)){
	if(!z$common)cat("\nFamily parameters:\n")
	cname <- c(cname,if(is.character(attr(z$family,"model")))
		attr(z$family,"model")
		else if(length(grep("linear",attr(z$family,"parameters")))>0)
		attr(z$family,"parameters")[grep("\\[",attr(z$family,"parameters"))]
		else attr(z$family,"parameters"))
	if(!is.null(z$linmodel[[3]]))cname <- c(cname,z$linmodel[[3]])
	if(z$common){
		dimnames(coef.table) <- list(unique(cname),c("estimate","se"))
		print.default(coef.table,digits=digits,print.gap=2)}
	else {
		coef.table <- cbind(z$coefficients[np3:np],z$se[np3:np])
		dimnames(coef.table) <- list(cname,c("estimate","se"))
		print.default(coef.table,digits=digits,print.gap=2)}}
if(z$npm>0&&is.null(z$mix)){
	cat("\nMixing dispersion parameter:\n")
	coef.table <- cbind(z$coefficients[np],z$se[np])
	dimnames(coef.table) <- list(" ", c("estimate", "se"))
	print.default(coef.table,digits=digits,print.gap=2)}
if(np>1&&correlation){
	cat("\nCorrelations:\n")
	dimnames(z$corr) <- list(seq(1,np),seq(1,np))
	print.default(z$corr,digits=digits)}
invisible(z)}
