/*
 *  rmutil : A Library of Special Functions for Repeated Measurements
 *  Copyright (C) 1998 J.K. Lindsey
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

void call_R(char *fcn, long nargs, void **arguments, char **modes,
	    long *lengths, char **names, long nres, char **results);

void interp(double x[], double fx[], int pts, double tab1[],
		   double tab2[], double *f, double *df, int *err)
{
  int i,j,ni=0;
  double diff1,diff2,tmp1,tmp2,lim1,lim2;
 
  *err=0;
  tmp1=fabs(x[0]);
  for(i=0;i<pts;i++){
    tmp2=fabs(x[i]);
    if(tmp2<tmp1){
      ni=i;
      tmp1=tmp2;}
    tab1[i]=fx[i];
    tab2[i]=fx[i];}
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
    *df=2*ni<(pts-j-3)?tab1[ni+1]:tab2[ni--];
    *f+=*df;}}

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
    for(k=0;k<len;k++)zz[k]=0.5*(a[k]+b[k]);
    call_R(fcn, 1L, args, mode, length, 0L, 1L, ss);
    x=(double *)ss[0];
    for(k=0;k<len;k++)sum[k]=(b[k]-a[k])*x[k];
    return;}
  else {
    for(i=1,j=1;j<n-1;j++) i*=3;
    nn=i;
    for(k=0;k<len;k++){
      pnt1[k]=(b[k]-a[k])/(3.0*nn);
      pnt2[k]=2.0*pnt1[k];
      zz[k]=a[k]+0.5*pnt1[k];
      tmpsum[k]=0.0;}
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
    for(k=0;k<len;k++)sum[k]=(sum[k]+(b[k]-a[k])*tmpsum[k]/nn)/3.0;
    return;}}

void romberg(void *fcn, double *a, double *b, int *len, double *eps,
	     int *pts, int *max, int *err, double sumlen[])
{
  int i,j,j1,finish;
  double errsum,*tab1,*tab2,*x,*fx,*sum,*tmpsum,*zz,*pnt1,*pnt2;

  tab1=(double*)malloc((size_t)(*pts*sizeof(double)));
  tab2=(double*)malloc((size_t)(*pts*sizeof(double)));
  x=(double*)malloc((size_t)((*max**len)*sizeof(double)));
  fx=(double*)malloc((size_t)((*max**len)*sizeof(double)));
  sum=(double*)malloc((size_t)(*len*sizeof(double)));
  tmpsum=(double*)malloc((size_t)(*len*sizeof(double)));
  zz=(double*)malloc((size_t)(*len*sizeof(double)));
  pnt1=(double*)malloc((size_t)(*len*sizeof(double)));
  pnt2=(double*)malloc((size_t)(*len*sizeof(double)));
  if(!tab1||!tab2||!x||!fx||!sum||!tmpsum||!zz||!pnt1||!pnt2){
    *err=1;
    return;}
  *err=0;
  for(i=0;i<*len;i++)x[i**max]=1.0;
  for(j=0;j<*max;j++){
    j1=j+1;
    evalRfn(fcn,a,b,j1,*len,sum,tmpsum,zz,pnt1,pnt2);
    finish=(j1>=*pts?1:0);
    for(i=0;i<*len;i++){
      fx[j+i**max]=sum[i];
      if(j1>=*pts){
	interp(&x[j1-*pts+i**max],&fx[j1-*pts+i**max],*pts,tab1,tab2,&sumlen[i],&errsum,err);
	if(*err)goto end;
	if(fabs(errsum)>*eps*fabs(sumlen[i]))finish=0;}
      x[j1+i**max]=x[j+i**max]/9.0;
      fx[j1+i**max]=fx[j+i**max];}
    if(finish)goto end;}
  *err=3;
 end: free((char *)x);
  free((char *)fx);
  free((char *)sum);
  free((char *)tmpsum);
  free((char *)zz);
  free((char *)pnt1);
  free((char *)pnt2);
  free((char *)tab2);
  free((char *)tab1);
  return;}
