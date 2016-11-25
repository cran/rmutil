#
#  rmutil : A Library of Special Functions for Repeated Measurements
#  Copyright (C) 2001 J.K. Lindsey
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
#     contr.mean(n, contrasts=TRUE)
#
#  DESCRIPTION
#
#    A function to provide correct constraints about the mean
#  (correcting contr.sum)

contr.mean <- function(n, contrasts=TRUE){
	if(length(n) <= 1){
		if(is.numeric(n)&&length(n)==1&&n>1)levels <- 1:n
		else stop("Not enough degrees of freedom to define contrasts")}
	else levels <- n
	lenglev <- length(levels)
	if(contrasts){
		cont <- array(0,c(lenglev,lenglev-1),
			list(levels,levels[1:(lenglev-1)]))
		cont[col(cont)==row(cont)] <- 1
		cont[lenglev,] <- -1}
	else {
		cont <- array(0,c(lenglev,lenglev),list(levels,levels))
		cont[col(cont) == row(cont)] <- 1}
	cont}
