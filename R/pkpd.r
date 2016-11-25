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
#     mu1.0o1c(p, times, dose=1, end=0.5)
#     mu1.1o1c(p, times, dose=1)
#     mu1.1o2c(p, times, dose=1)
#     mu1.1o2cl(p, times, dose=1)
#     mu1.1o2cc(p, times, dose=1)
#     mu2.0o1c(p, times, dose=1, ind, end=0.5)
#     mu2.0o2c1(p, times, dose=1, ind, end=0.5)
#     mu2.0o2c2(p, times, dose=1, ind, end=0.5)
#     mu2.1o1c(p, times, dose=1, ind)
#     mu2.0o1cfp(p, times, dose=1, ind, end=0.5)
#     mu2.0o2c1fp(p, times, dose=1, ind, end=0.5)
#     mu2.0o2c2fp(p, times, dose=1, ind, end=0.5)
#     mu2.1o1cfp(p, times, dose=1, ind)
#
#  DESCRIPTION
#
#    Functions giving nonlinear regressions models for various PKPD
#  compartment models

### standard pharmacokinetic compartment models
###
### open zero-order one-compartment model
# p[1]: log volume (V)
# p[2]: log elimination rate (ke)
# end:  time when infusion stops
mu1.0o1c <- function(p, times, dose=1, end=0.5) {
	ke <- exp(p[2])
	dose/(exp(p[1])*ke)*((1-exp(-ke*times))*(times<=end)+
		(1-exp(-ke*end))*exp(-ke*(times-end))*(times>end))}
###
### open first-order one-compartment model
# p[1]: log volume (V)
# p[2]: log absorption rate (ka)
# p[3]: log elimination rate (ke)
mu1.1o1c <- function(p, times, dose=1) {
	ka <- exp(p[2])
	ke <- exp(p[3])
	exp(p[2]-p[1])*dose/(ka-ke)*(exp(-ke*times)-exp(-ka*times))}
###
### open first-order two-compartment model (ordered)
# p[1]: log volume (V)
# p[2]: log absorption rate (ka)
# p[3]: log elimination rate (ke)
# p[4]: log transfer rate between compartments (k12)
mu1.1o2c <- function(p, times, dose=1) {
	ka <- exp(p[2])
	ke <- exp(p[3])
	k12 <- exp(p[4])
	ka*k12*exp(-p[1])*dose/(k12-ka)*((exp(-ka*times)-exp(-ke*times))/
		(ke-ka)-(exp(-k12*times)-exp(-ke*times))/(ke-k12))}
###
### open first-order two-compartment model (ordered, absorption and transfer equal)
# p[1]: log volume (V)
# p[2]: log absorption rate (ka)
# p[3]: log elimination rate (ke)
mu1.1o2cl <- function(p, times, dose=1) {
	ka <- exp(p[2])
	ke <- exp(p[3])
	ka^2*exp(-p[1])*dose/(ka-ke)*((exp(-ka*times)-exp(-ke*times))/(ke-ka)
		-times*exp(-ka*times))}
###
### open first-order two-compartment model (circular)
# p[1]: log volume (V)
# p[2]: log absorption rate (ka)
# p[3]: log elimination rate (ke)
# p[4]: log rate to second compartment (k12)
# p[5]: log rate from second compartment (k21)
mu1.1o2cc <- function(p, times, dose=1) {
	ka <- exp(p[2])
	ke <- exp(p[3])
	k12 <- exp(p[4])
	k21 <- exp(p[5])
	beta <- 0.5*(k12+k21+ke-sqrt((k12+k21+ke)^2-4*k21*ke))
	alpha <- (k21*ke)/beta
	exp(p[2]-p[1])*dose*((k21-alpha)*exp(-alpha*times)/
		((ka-alpha)*(beta-alpha))+(k21-beta)*exp(-beta*times)/
		((ka-beta)*(alpha-beta))+(k21-ka)*exp(-ka*times)/
		((beta-ka)*(beta-ka)))}
###
### simultaneous models for parent drug and metabolite
###
### zero-order one-compartment model
# p[1]: log parent drug volume (Vp)
# p[2]: log parent drug direct elimination rate (kep)
# p[3]: log transformation rate from parent to metabolite (kpm)
# p[4]: log metabolite elimination rate (kem)
# p[5]: log metabolite volume (Vm)
# ind:  indicator vector: 1 for parent, 0 for metabolite
# end:  time when infusion stops
mu2.0o1c <- function(p, times, dose=1, ind, end=0.5) {
	Vp <- exp(p[1])
	kpm <- exp(p[3])
	kp <- exp(p[2])+kpm
	kem <- exp(p[4])
	Vm <- exp(p[5])
	kemp <- kem-kp
	tmp1 <- exp(-kp*times)
	tmp2 <- kpm/(kp*kem*Vm)
	g1 <- exp(-kp*end)
	g2 <- exp(-kem*end)
	cend <- (1-g1)/(Vp*kp)
	cexp <- exp(-kp*(times-end))*(times>end)
	cmend <- tmp2*(1+kp/kemp*g2-kem/kemp*g1)
	tmp3 <- cend*kpm*Vp/(kemp*Vm)
	dose*(ind*((1-tmp1)/(Vp*kp)*(times<=end)+cend*cexp)+
		(1-ind)*(tmp2*(1+kp/kemp*exp(-kem*times)-kem/kemp*
		tmp1)*(times<=end)+(g2/g1*cexp*tmp3+
		(cmend-tmp3)*exp(-kem*(times-end))/g2*(times>end))))}
###
### zero-order two-compartment for parent, one-compartment for
###  metabolite, model
# p[1]: log parent drug volume (Vp)
# p[2]: log parent drug direct elimination rate (kep)
# p[3]: log parent drug rate to second compartment (k12)
# p[4]: log parent drug rate from second compartment (k21)
# p[5]: log transformation rate from parent to metabolite (kpm)
# p[6]: log metabolite elimination rate (kem)
# ind:  indicator vector: 1 for parent, 0 for metabolite
# end:  time when infusion stops
mu2.0o2c1 <- function(p, times, dose=1, ind, end=0.5) {
	Vp <- exp(p[1])
	kp12 <- exp(p[3])
	kp21 <- exp(p[4])
	kpm <- exp(p[5])
	kp <- exp(p[2])+kpm
	kem <- exp(p[6])
	tmp1 <- sqrt((kp+kp12+kp21)^2-4*kp21*kp)
	lamp1 <- 0.5*(kp+kp12+kp21+tmp1)
	lamp2 <- 0.5*(kp+kp12+kp21-tmp1)
	tmp10 <- exp(-kem*times)
	tmp13 <- exp(-kem*end)
	tmp2 <- (1-tmp10)/kem
	tmp3 <- (1-tmp13)/kem
	tmp4 <- lamp1-kp21
	tmp5 <- lamp2-kp21
	tmp6 <- exp(-lamp1*times)
	tmp7 <- exp(-lamp2*times)
	tmp8 <- exp(-lamp1*end)
 	tmp9 <- exp(-lamp2*end)
	tmp11 <- tmp6/tmp8
	tmp12 <- tmp7/tmp9
	tmp14 <- tmp10/tmp13
	dose/(Vp*tmp1)*(ind*((tmp4*(1-tmp6)/lamp1-tmp5*(1-tmp7)/lamp2)*
		(times<=end)+(tmp4*(1-tmp8)*tmp11/lamp1-tmp5*(1-tmp9)*
		tmp12/lamp2)*(times>end))+(1-ind)*kpm*((tmp4*(tmp2-
		(tmp6-tmp10)/(kem-lamp1))/lamp1-tmp5*(tmp2-(tmp7-tmp10)/
		(kem-lamp2))/lamp2)*(times<=end)+((tmp4*(tmp3-(tmp8-
		tmp13)/(kem-lamp1))/lamp1-tmp5*(tmp3-(tmp9-tmp13)/
		(kem-lamp2))/lamp2)*tmp14+tmp4*(1-tmp8)*
		(tmp14-tmp11)/(lamp1*(lamp1-kem))-tmp5*(1-tmp9)*
		(tmp14-tmp12)/(lamp2*(lamp2-kem)))*(times>end)))}
###
### zero-order two-compartment model for both parent and metabolite
# p[1]: log parent drug volume (Vp)
# p[2]: log parent drug direct elimination rate (kep)
# p[3]: log parent drug rate to second compartment (kp12)
# p[4]: log parent drug rate from second compartment (kp21)
# p[5]: log transformation rate from parent to metabolite (kpm)
# p[6]: log metabolite elimination rate (kem)
# p[7]: log metabolite drug rate to second compartment (km12)
# p[8]: log metabolite drug rate from second compartment (km21)
# ind:  indicator vector: 1 for parent, 0 for metabolite
# end:  time when infusion stops
mu2.0o2c2 <- function(p, times, dose=1, ind, end=0.5) {
        Vp <- exp(p[1])
	kp12 <- exp(p[3])
	kp21 <- exp(p[4])
	kpm <- exp(p[5])
	kp <- exp(p[2])+kpm
	kem <- exp(p[6])
	km12 <- exp(p[7])
	km21 <- exp(p[8])
	tmp1 <- sqrt((kp+kp12+kp21)^2-4*kp21*kp)
	lamp1 <- 0.5*(kp+kp12+kp21+tmp1)
	lamp2 <- 0.5*(kp+kp12+kp21-tmp1)
	tmp2 <- lamp1-kp21
	tmp3 <- lamp2-kp21
	tmp4 <- exp(-lamp1*times)
	tmp5 <- exp(-lamp2*times)
	tmp6 <- exp(-lamp1*end)
	tmp7 <- exp(-lamp2*end)
	tmp8 <- tmp4/tmp6
	tmp9 <- tmp5/tmp7
	tmp10 <- sqrt((kem+km12+km21)^2-4*km21*kem)
	lamm1 <- 0.5*(kem+km12+km21+tmp10)
	lamm2 <- 0.5*(kem+km12+km21-tmp10)
	tmp11 <- kem-lamm1
	tmp12 <- kem-lamm2
        tmp13 <- lamp1-lamm2
        tmp14 <- lamp1-lamm1
	tmp15 <- exp(-lamm1*times)
	tmp16 <- exp(-lamm2*times)
        tmp17 <- lamp1-km21
        tmp18 <- lamp2-km21
        tmp19 <- lamp2-lamm1
        tmp20 <- lamp2-lamm2
	tmp21 <- exp(-lamm1*end)
	tmp22 <- exp(-lamm2*end)
        tmp23 <- km21-lamp1
        tmp24 <- km21-lamp2
        tmp25 <- tmp15/tmp21
        tmp26 <- tmp16/tmp22
	dose/Vp*(ind*((tmp2*(1-tmp4)/(lamp1*tmp1)-tmp3*(1-tmp5)/
		(lamp2*tmp1))*(times<=end)+(tmp2*(1-tmp6)*tmp8/(lamp1*tmp1)-
		tmp3*(1-tmp7)* tmp9/(lamp2*tmp1))*(times>end))+(1-ind)*kpm* 
		((tmp2/tmp1*((tmp11*tmp16/tmp13-tmp12*tmp15/tmp14)/
		(kem*tmp10)+tmp17*tmp4/(lamp1*tmp14*tmp13)+1/(kem*lamp1))
		-tmp3/tmp1*((tmp12*tmp15/tmp19-tmp11*tmp16/tmp20)/
		(kem*(-tmp10))+tmp18*tmp5/(lamp2*tmp19*tmp20)+1/(kem*lamp2)))*
		(times<=end)+((-tmp2/(tmp1*tmp10)*(tmp12*tmp21/(kem*tmp14)+
		lamm2/(kem*lamp1)+tmp23/(lamp1*tmp14))+tmp3/(tmp1*tmp10)*
		(tmp12*tmp21/(kem*tmp19)+lamm2/(kem*lamp2)+tmp24/
		(lamp2*tmp19)))*tmp25+(tmp2/(tmp1*tmp10)*(tmp11*tmp22/
		(kem*tmp13)+lamm1/(kem*lamp1)+tmp23/(lamp1*tmp13))-
		tmp3/(tmp1*tmp10)*( tmp11*tmp22/(kem*tmp20)+lamm1/(kem*lamp2)+
		tmp24/(lamp2*tmp20)))*tmp26+tmp2*tmp23*(1-tmp6)*tmp8/
		(lamp1*tmp1*tmp14*tmp13)-tmp3*tmp24*(1-tmp7)*tmp9/
		(lamp2*tmp1*tmp19*tmp20))*(times>end)))}
###
### first-order one-compartment model
# p[1]: log volume (V)
# p[2]: log parent drug absorption rate (kap)
# p[3]: log parent drug direct elimination rate (kep)
# p[4]: log transformation rate from parent to metabolite (kpm)
# p[5]: log metabolite elimination rate (kem)
# ind:  indicator vector: 1 for parent, 0 for metabolite
mu2.1o1c <- function(p, times, dose=1, ind) {
	kap <- exp(p[2])
	kep <- exp(p[3])
	kem <- exp(p[5])
	kap*exp(p[1])*dose/(kap-kep)*(ind*(exp(-kep*times)-exp(-kap*times))+
	(1-ind)*exp(p[4])*(exp(-kap*times)/(kap-kem)-exp(-kep*times)/(kep-kem)+
	(1/(kep-kem)-1/(kap-kem))*exp(-kem*times)))}
###
### zero-order one-compartment first-pass model
# p[1]: log parent drug volume (Vp)
# p[2]: log parent drug direct elimination rate (kep)
# p[3]: log transformation rate from parent to metabolite (kpm)
# p[4]: log metabolite elimination rate (kem)
# p[5]: log metabolite volume (Vm)
# p[7]: logit of proportion going to first pass (lpfp)
# ind:  indicator vector: 1 for parent, 0 for metabolite
# end:  time when infusion stops
mu2.0o1cfp <- function(p, times, dose=1, ind, end=0.5) {
	Vp <- exp(p[1])
	kpm <- exp(p[3])
	kp <- exp(p[2])+kpm
	kem <- exp(p[4])
	Vm <- exp(p[5])
	kemp <- kem-kp
	tmp1 <- exp(-kp*times)
	tmp2 <- exp(-kem*times)
	tmp3 <- kpm/(kp*kem*Vm)
	g1 <- exp(-kp*end)
	g2 <- exp(-kem*end)
	cend <- (1-g1)/(Vp*kp)
	cexp <- exp(-kp*(times-end))*(times>end)
	cmend <- tmp3*(1+kp/kemp*g2-kem/kemp*g1)
	tmp4 <- cend*kpm*Vp/(kemp*Vm)
	lpfp <- exp(p[6])
	lpfp <- lpfp/(1+lpfp)
	dose*(ind*((1-tmp1)/(Vp*kp)*(times<=end)+cend*cexp)*lpfp+
		(1-ind)*(((1-tmp2)/(Vm*kem)*(times<=0.5)+(1-g2)/(Vm*kem)*
		exp(-kem*(times-0.5))*(times>0.5))*(1-lpfp)+(tmp3*
		(1+kp/kemp*tmp2-kem/kemp*tmp1)*(times<=end)+(g2/g1*cexp*tmp4+
		(cmend-tmp4)*exp(-kem*(times-end))/g2*(times>end)))*lpfp))}
###
### zero-order two-compartment for parent, one-compartment for metabolite,
###   first-pass model
# p[1]: log parent drug volume (Vp)
# p[2]: log parent drug direct elimination rate (kep)
# p[3]: log parent drug rate to second compartment (k12)
# p[4]: log parent drug rate from second compartment (k21)
# p[5]: log transformation rate from parent to metabolite (kpm)
# p[6]: log metabolite elimination rate (kem)
# p[7]: logit of proportion going to first pass (lpfp)
# ind:  indicator vector: 1 for parent, 0 for metabolite
# end:  time when infusion stops
mu2.0o2c1fp <- function(p, times, dose=1, ind, end=0.5) {
	Vp <- exp(p[1])
	kp12 <- exp(p[3])
	kp21 <- exp(p[4])
	kpm <- exp(p[5])
	kp <- exp(p[2])+kpm
	kem <- exp(p[6])
	tmp1 <- sqrt((kp+kp12+kp21)^2-4*kp21*kp)
	lamp1 <- 0.5*(kp+kp12+kp21+tmp1)
	lamp2 <- 0.5*(kp+kp12+kp21-tmp1)
	tmp10 <- exp(-kem*times)
	tmp13 <- exp(-kem*end)
	tmp2 <- (1-tmp10)/kem
	tmp3 <- (1-tmp13)/kem
	tmp4 <- lamp1-kp21
	tmp5 <- lamp2-kp21
	tmp6 <- exp(-lamp1*times)
	tmp7 <- exp(-lamp2*times)
	tmp8 <- exp(-lamp1*end)
 	tmp9 <- exp(-lamp2*end)
	tmp11 <- tmp6/tmp8
	tmp12 <- tmp7/tmp9
	tmp14 <- tmp10/tmp13
	lpfp <- exp(p[7])
	lpfp <- lpfp/(1+lpfp)
	dose/(Vp*tmp1)*(ind*((tmp4*(1-tmp6)/lamp1-tmp5*(1-tmp7)/lamp2)*
		(times<=end)+(tmp4*(1-tmp8)*tmp11/lamp1-tmp5*(1-tmp9)*
		tmp12/lamp2)*(times>end))*lpfp+(1-ind)*(((1-exp(-kem*times))/
		kem*(times<=end)+(1-exp(-kem*end))/kem*exp(-kem*(times-end))*
		(times>end))*tmp1*(1-lpfp)+(kpm*((tmp4*(tmp2-(tmp6-tmp10)/
		(kem-lamp1))/lamp1-tmp5*(tmp2-(tmp7-tmp10)/(kem-lamp2))/
		lamp2)*(times<=end)+((tmp4*(tmp3-(tmp8-
		tmp13)/(kem-lamp1))/lamp1-tmp5*(tmp3-(tmp9-tmp13)/
		(kem-lamp2))/lamp2)*tmp14+tmp4*(1-tmp8)*
		(tmp14-tmp11)/(lamp1*(lamp1-kem))-tmp5*(1-tmp9)*
		(tmp14-tmp12)/(lamp2*(lamp2-kem)))*(times>end)))*lpfp))}
###
### zero-order two-compartment model for both parent and metabolite
###   first-pass model
# p[1]: log parent drug volume (Vp)
# p[2]: log parent drug direct elimination rate (kep)
# p[3]: log parent drug rate to second compartment (kp12)
# p[4]: log parent drug rate from second compartment (kp21)
# p[5]: log transformation rate from parent to metabolite (kpm)
# p[6]: log metabolite elimination rate (kem)
# p[7]: log metabolite drug rate to second compartment (km12)
# p[8]: log metabolite drug rate from second compartment (km21)
# p[9]: logit of proportion going to first pass (lpfp)
# ind:  indicator vector: 1 for parent, 0 for metabolite
# end:  time when infusion stops
mu2.0o2c2fp <- function(p, times, dose=1, ind, end=0.5) {
        Vp <- exp(p[1])
	kp12 <- exp(p[3])
	kp21 <- exp(p[4])
	kpm <- exp(p[5])
	kp <- exp(p[2])+kpm
	kem <- exp(p[6])
	km12 <- exp(p[7])
	km21 <- exp(p[8])
	tmp1 <- sqrt((kp+kp12+kp21)^2-4*kp21*kp)
	lamp1 <- 0.5*(kp+kp12+kp21+tmp1)
	lamp2 <- 0.5*(kp+kp12+kp21-tmp1)
	tmp2 <- lamp1-kp21
	tmp3 <- lamp2-kp21
	tmp4 <- exp(-lamp1*times)
	tmp5 <- exp(-lamp2*times)
	tmp6 <- exp(-lamp1*end)
	tmp7 <- exp(-lamp2*end)
	tmp8 <- tmp4/tmp6
	tmp9 <- tmp5/tmp7
	tmp10 <- sqrt((kem+km12+km21)^2-4*km21*kem)
	lamm1 <- 0.5*(kem+km12+km21+tmp10)
	lamm2 <- 0.5*(kem+km12+km21-tmp10)
	tmp11 <- kem-lamm1
	tmp12 <- kem-lamm2
        tmp13 <- lamp1-lamm2
        tmp14 <- lamp1-lamm1
	tmp15 <- exp(-lamm1*times)
	tmp16 <- exp(-lamm2*times)
        tmp17 <- lamp1-km21
        tmp18 <- lamp2-km21
        tmp19 <- lamp2-lamm1
        tmp20 <- lamp2-lamm2
	tmp21 <- exp(-lamm1*end)
	tmp22 <- exp(-lamm2*end)
        tmp23 <- km21-lamp1
        tmp24 <- km21-lamp2
        tmp25 <- tmp15/tmp21
        tmp26 <- tmp16/tmp22
	lpfp <- exp(p[9])
	lpfp <- lpfp/(1+lpfp)
	dose/Vp*(ind*((tmp2*(1-tmp4)/(lamp1*tmp1)-tmp3*(1-tmp5)/(lamp2*tmp1))*
		(times<=end)+(tmp2*(1-tmp6)*tmp8/(lamp1*tmp1)-tmp3*(1-tmp7)*
		tmp9/(lamp2*tmp1))*(times>end))*lpfp+(1-ind)*
		(((1-exp(-kem*times))/kem*(times<=end)+(1-exp(-kem*end))/
		kem*exp(-kem*(times-end))*(times>end))*(1-lpfp)+
		(kpm*((tmp2/tmp1*((tmp11*tmp16/tmp13-tmp12*tmp15/tmp14)/
		(kem*tmp10)+tmp17*tmp4/(lamp1*tmp14*tmp13)+1/(kem*lamp1))-
		tmp3/tmp1*((tmp12*tmp15/tmp19-tmp11*tmp16/tmp20)/
		(kem*(-tmp10))+tmp18*tmp5/(lamp2*tmp19*tmp20)+1/(kem*lamp2)))*
		(times<=end)+((-tmp2/(tmp1*tmp10)*(tmp12*tmp21/(kem*tmp14)+
		lamm2/(kem*lamp1)+tmp23/(lamp1*tmp14))+tmp3/(tmp1*tmp10)*
		(tmp12*tmp21/(kem*tmp19)+lamm2/(kem*lamp2)+tmp24/
		(lamp2*tmp19)))*tmp25+(+tmp2/(tmp1*tmp10)*(tmp11*tmp22/
		(kem*tmp13)+lamm1/(kem*lamp1)+tmp23/(lamp1*tmp13))-tmp3/
		(tmp1*tmp10)*(tmp11*tmp22/(kem*tmp20)+lamm1/(kem*lamp2)+
		tmp24/(lamp2*tmp20)))*tmp26+tmp2*tmp23*(1-tmp6)*tmp8/
		(lamp1*tmp1*tmp14*tmp13)-tmp3*tmp24*(1-tmp7)*tmp9/
		(lamp2*tmp1*tmp19*tmp20))*(times>end)))*lpfp))}
###
### first-order one-compartment first-pass model
# p[1]: log volume (V)
# p[2]: log parent drug absorption rate (kap)
# p[3]: log parent drug direct elimination rate (ked)
#	(kep: total parent drug elimination)
# p[4]: log transformation rate from parent to metabolite (kpm)
# p[5]: log metabolite first-pass formation rate (kmfp)
# p[6]: log metabolite elimination rate (kem)
# p[7]: logit of proportion going to first pass (lpfp)
# ind:  indicator vector: 1 for parent, 0 for metabolite
mu2.1o1cfp <- function(p, times, dose=1, ind) {
	kap <- exp(p[2])
	kf <- exp(p[4])
	kep <- exp(p[3])+kf
	kmfp <- exp(p[5])
	kem <- exp(p[6])
	lpfp <- exp(p[7])
	lpfp <- lpfp/(1+lpfp)
	exp(-p[1])*dose*(ind*(exp(-kep*times)-exp(-kap*times))/
		(kap-kep)*lpfp+(1-ind)*((exp(-kmfp*times)-exp(-kap*times))/
		(kap-kmfp)*(1-lpfp)+kf*(exp(-kap*times)/(kap-kem)-
		exp(-kep*times)/(kep-kem)+(1/(kep-kem)-1/(kap-kem))*
		exp(-kem*times))/(kap-kep)*lpfp))}
