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
#     read.list(file="", skip=0, nlines=2, order=NULL)
#     read.surv(file="", skip=0, nlines=1, cumulative=T, all=T)
#
#  DESCRIPTION
#
#    Utility functions for reading repeated measurements data

read.list <- function(file="", skip=0, nlines=2, order=NULL){
	if(!is.null(order)){
		if(length(order)!=nlines)
			stop("order must have length",nlines,"\n")
		else if(range(order)!=c(1,nlines))
			stop("order must have values in",1:nlines,"\n")}
	continue <- T
	result <- list()
	while(continue){
		x <- scan(file,skip=skip,nlines=nlines,quiet=T)
		skip <- skip+nlines
		if(length(x)==0)continue <- F
		else {
			tmp <- matrix(x,ncol=nlines)
			if(!is.null(order))tmp <- tmp[,order]
			result <- c(result,list(tmp))}}
	invisible(result)}

read.surv <- function(file="", skip=0, nlines=1, cumulative=T, all=T){
	continue <- T
	result <- list()
	censor <- NULL
	while(continue){
		x <- scan(file,skip=skip,nlines=nlines,quiet=T)
		skip <- skip+nlines
		if(length(x)==0)continue <- F
		else {
			if(all)mm <- matrix(x,ncol=2,byrow=T)[,1]
			else mm <- x[1:(length(x)-1)]
			if(cumulative)mm <- c(mm[1],diff(mm))
			result <- c(result,list(mm))
			censor <- c(censor,x[length(x)])}}
	invisible(list(result,censor))}
