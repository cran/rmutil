/*
 *  rmutil : A Library of Special Functions for Repeated Measurements
 *  Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  SYNOPSIS
 *
 *    void romberg(void *fcn, double *a, double *b, int *len, double *eps,
 *	      int *pts, int *max, int *err, double sumlen[])
 *
 *  DESCRIPTION
 *
 *    This function performs vectorized Romberg integration for the
 *    R function supplied.
 *
 */

#include <math.h>
#include <stddef.h>
#include "R.h"

/* polynomial interpolation algorithm for function with values fx
   evaluated at x */
void interp(double x[], double fx[], int pts, double tab1[],
		   double tab2[], double *f, double *df, int *err)
{
  int i,j,ni=0;
  double diff1,diff2,tmp1,tmp2,lim1,lim2;
 
  *err=0;
  /* create an initial table of values from those supplied */
  tmp1=fabs(x[0]);
  for(i=0;i<pts;i++){
    tmp2=fabs(x[i]);
    if(tmp2<tmp1){
      ni=i;
      tmp1=tmp2;}
    tab1[i]=tab2[i]=fx[i];}
  /* find the best interpolated value by modifying the table */
  *f=fx[ni--];
  for(j=0;j<pts-1;j++){
    for(i=0;i<pts-j-1;i++){
      lim1=x[i];
      lim2=x[i+j+1];
      diff1=tab1[i+1]-tab2[i];
      diff2=lim1-lim2;
      if(diff2==0.0){
	*err=2;
	return;}
      diff2=diff1/diff2;
      tab2[i]=lim2*diff2;
      tab1[i]=lim1*diff2;}
    /* calculate estimated error and final interpolated value */
    *df=2*ni<(pts-j-3)?tab1[ni+1]:tab2[ni--];
    *f+=*df;}}

/* evaluate the R function to be integrated using a trapezoid rule
   between limits a and b.
   At each step, 2^(n-1) interior points are added. */
static void evalRfn(void *fcn, double a[], double b[], int n, int len,
		double sum[], double tmpsum[], double zz[],
		double pnt1[], double pnt2[])
{
  double *x,nn;
  char *mode[1],*ss[1];
  long length[1];
  void *args[1];
  int i,j,k;

  mode[0] = "double";
  length[0] = len;
  args[0] = (void *)(zz);

  if(n==1){
    /* one trapezoid */
    for(k=0;k<len;k++)zz[k]=0.5*(a[k]+b[k]);
    call_R(fcn, 1L, args, mode, length, 0L, 1L, ss);
    x=(double *)ss[0];
    for(k=0;k<len;k++)sum[k]=(b[k]-a[k])*x[k];
    return;}
  else {
    /* several trapezoids */
    for(i=1,j=1;j<n-1;j++) i*=3;
    nn=i;
    /* calculate points for evaluation */
    for(k=0;k<len;k++){
      pnt1[k]=(b[k]-a[k])/(3.0*nn);
      pnt2[k]=2.0*pnt1[k];
      zz[k]=a[k]+0.5*pnt1[k];
      tmpsum[k]=0.0;}
    /* evaluate R function at these points */
    for(j=1;j<=i;j++){
      call_R(fcn, 1L, args, mode, length, 0L, 1L, ss);
      x=(double *)ss[0];
      for(k=0;k<len;k++){
	tmpsum[k]+=x[k];
	zz[k]+=pnt2[k];}
      call_R(fcn, 1L, args, mode, length, 0L, 1L, ss);
      x=(double *)ss[0];
      for(k=0;k<len;k++){
	tmpsum[k]+=x[k];
	zz[k]+=pnt1[k];}}
    /* calculate total area */
    for(k=0;k<len;k++)sum[k]=(sum[k]+(b[k]-a[k])*tmpsum[k]/nn)/3.0;
    return;}}

void romberg(void *fcn, double *a, double *b, int *len, double *eps,
	     int *pts, int *max, int *err, double sumlen[])
{
  int i,j,j1,finish;
  double errsum,*tab1,*tab2,*x,*fx,*sum,*tmpsum,*zz,*pnt1,*pnt2;

  /* allocate vectors acccording to number of points and length of vector */
  tab1=(double*)R_alloc((size_t)(*pts),sizeof(double));
  tab2=(double*)R_alloc((size_t)(*pts),sizeof(double));
  x=(double*)R_alloc((size_t)(*max*(*len+1)),sizeof(double));
  fx=(double*)R_alloc((size_t)(*max*(*len+1)),sizeof(double));
  sum=(double*)R_alloc((size_t)(*len),sizeof(double));
  tmpsum=(double*)R_alloc((size_t)(*len),sizeof(double));
  zz=(double*)R_alloc((size_t)(*len),sizeof(double));
  pnt1=(double*)R_alloc((size_t)(*len),sizeof(double));
  pnt2=(double*)R_alloc((size_t)(*len),sizeof(double));
  if(!tab1||!tab2||!x||!fx||!sum||!tmpsum||!zz||!pnt1||!pnt2){
    *err=1;
    return;}
  *err=0;
  for(i=0;i<*len;i++)x[i**max]=1.0;
  /* iterate, decreasing step size, until convergence or max number of steps */
  for(j=0;j<*max;j++){
    j1=j+1;
    /* evaluate R function and calculate sum of trapezoids */
    evalRfn(fcn,a,b,j1,*len,sum,tmpsum,zz,pnt1,pnt2);
    finish=(j1>=*pts?1:0);
    /* repeatedly call polynomial interpolation routine */
    for(i=0;i<*len;i++){
      fx[j+i**max]=sum[i];
      if(j1>=*pts){
	interp(&x[j1-*pts+i**max],&fx[j1-*pts+i**max],*pts,tab1,tab2,&sumlen[i],&errsum,err);
	if(*err)return;
	/* check convergence */
	if(fabs(errsum)>*eps*fabs(sumlen[i]))finish=0;}
      /* decrease step size */
      x[j1+i**max]=x[j+i**max]/9.0;
      fx[j1+i**max]=fx[j+i**max];}
    if(finish)return;}
  *err=3;
  return;}
