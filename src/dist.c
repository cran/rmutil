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
 * void pdp(int q[], int *my, double m[], double s[], int *nn, double res[])
 * void ddp(int y[], int *my, double m[], double s[], int *nn,
 *	 double wt[], double res[])
 * void pmp(int q[], int *my, double m[], double s[], int *nn, double res[])
 * void dmp(int y[], int *my, double m[], double s[], int *nn,
 *	 double wt[], double res[])
 * void ppvfp(int q[], double m[], double s[], double f[], int *nn,
 *       double res[])
 * void dpvfp(int y[], double m[], double s[], double f[], int *nn,
 *	 double wt[], double res[])
 * void pdb(int q[], int n[], double m[], double s[], int *nn, double res[])
 * void ddb(int y[], int n[], double m[], double s[], int *nn,
 *	 double wt[], double res[])
 * void pmb(int q[], int n[], double m[], double s[], int *nn, double res[])
 * void dmb(int y[], int n[], double m[], double s[], int *nn,
 *	 double wt[], double res[])
 *
 * void psimplex_c(double y[], double m[], double s[], double f[], int *len,
 *	   double *eps, int *pts, int *max, int *err, double res[])
 * void pginvgauss_c(double y[], double m[], double s[], double f[], int *len,
 *	   double *eps, int *pts, int *max, int *err, double res[])
 * void ppowexp_c(double y[], double m[], double s[], double f[], int *len,
 *	   double *eps, int *pts, int *max, int *err, double res[])
 *
 *  DESCRIPTION
 *
 *    Functions to compute the probability and cumulative probability
 * functions of the following overdispersed discrete distributions:
 * double Poisson, multiplicative Poisson, double binomial, and
 * multiplicative binomial,
 * and the cumulative probability functions for the following
 * continuous distributions:
 * Levy, generalized inverse Gaussian, and power exponential.
 *
 */

#include <math.h>
#include <stddef.h>
#include "dist.h"
#include "R.h"
#include "Rmath.h"

/* double Poisson */
static double dpnc(int my, double m, double s){
  int i;
  double r;
  r=exp(-s*m);
  for(i=1;i<=my;i++)
    r+=exp(i*(1-s)*log((double)i)+i*s*log(m)+i*(s-1)-s*m-lgammafn(i+1.));
  return(r);}

void pdp(int q[], int *my, double m[], double s[], int *nn, double res[]){
  int i;
  for(i=0;i<*nn;i++)
    res[i]=dpnc(q[i],m[i],s[i])/dpnc(*my,m[i],s[i]);}

void ddp(int y[], int *my, double m[], double s[], int *nn,
	 double wt[], double res[]){
  int i,y1;
  for(i=0;i<*nn;i++){
    if(wt[i]>0){
      y1=y[i]>0?y[i]:1;
      res[i]=wt[i]*(-s[i]*m[i]+y[i]*s[i]*(1+log(m[i]/y1))+y[i]*log((double)y1)-y[i]-lgammafn(y[i]+1.)-log(dpnc(*my,m[i],s[i])));}
    else res[i]=0;}}

/* multiplicative Poisson */
static double mpnc(int my, double m, double s){
  int i;
  double r;
  r=0.;
  for(i=0;i<=my;i++)r+=exp(i*log(m)+i*i*s-m-lgammafn(i+1.));
  return(r);}

void pmp(int q[], int *my, double m[], double s[], int *nn, double res[]){
  int i;
  double ss;
  for(i=0;i<*nn;i++){
    ss=log(s[i]);
    res[i]=mpnc(q[i],m[i],ss)/mpnc(*my,m[i],ss);}}

void dmp(int y[], int *my, double m[], double s[], int *nn,
	 double wt[], double res[]){
  int i;
  double ss;
  for(i=0;i<*nn;i++){
    if(wt[i]>0){
      ss=log(s[i]);
      res[i]=wt[i]*(-m[i]+y[i]*y[i]*ss+y[i]*log(m[i])-lgammafn(y[i]+1)-log(mpnc(*my,m[i],ss)));}
    else res[i]=0;}}

/* power variance function Poisson */

static double pvfc(int y, double m, double s, double f){
  int i,j;
  double r,*c,tmp1,tmp2,tmp3,tmp4;
  c=(double*)R_alloc((size_t)(y*y),sizeof(double));
  tmp1=gammafn(1.-f);
  tmp2=log(m);
  tmp3=log(s+1.);
  tmp4=log(s);
  for(i=0;i<y;i++){
    c[i*(y+1)]=1;
    if(i>0){
      c[i*y]=gammafn(i+1-f)/tmp1;
      if(i>1)
	for(j=1;j<i;j++)c[y*i+j]=c[y*(i-1)+j-1]+c[y*(i-1)+j]*(i-(j+1)*f);}}
  r=0.;
  for(i=1;i<=y;i++){
    r+=c[y*(y-1)+i-1]*exp(i*tmp2+(i*f-y)*tmp3-(i*(f-1))*tmp4);}
  return(r);}

void dpvfp(int y[], double m[], double s[], double f[], int *nn,
	 double wt[], double res[]){
  int i;
  for(i=0;i<*nn;i++){
    if(wt[i]>0){
      if(f[i]==0.)res[i]=dnbinom(y[i],m[i]*s[i],s[i]/(1+s[i]),0);
      else {
	res[i]=wt[i]*exp(-m[i]*((s[i]+1.)*pow((s[i]+1.)/s[i],f[i]-1.)-s[i])/f[i]);
	if(y[i]>0)res[i]*=pvfc(y[i],m[i],s[i],f[i]);
	if(y[i]>1)res[i]/=gammafn(y[i]+1);}}
    else res[i]=0;}}

void ppvfp(int q[], double m[], double s[], double f[], int *nn, double res[]){
  int i,j;
  static int k=1;
  static double wt=1;
  double tmp;
  for(i=0;i<*nn;i++){
    if(f[i]==0.)res[i]=pnbinom(q[i],m[i]*s[i],s[i]/(1+s[i]),1,0);
    else {
      res[i]=0.;
      for(j=0;j<q[i];j++){
	dpvfp(&j,&(m[i]),&(s[i]),&(f[i]),&k,&wt,&tmp);
	res[i]+=tmp;}}}}

/* double binomial */
static double dbnc(int yy, int n, double m, double s){
  int y;
  double r;
  r=0.;
  for(y=0;y<=yy;y++)
    r+=exp(lchoose((double)n,(double)y)+n*(s-1)*log((double)n)
	   +y*s*log(m)+((n-y)*s)*log(1.-m)-(y>0?y*(s-1)*log((double)y):0)
	   -(y<n?(n-y)*(s-1)*log((double)(n-y)):0));
  return(r);}

void pdb(int q[], int n[], double m[], double s[], int *nn, double res[]){
  int i;
  for(i=0;i<*nn;i++)
    res[i]=dbnc(q[i],n[i],m[i],s[i])/dbnc(n[i],n[i],m[i],s[i]);}

void ddb(int y[], int n[], double m[], double s[], int *nn,
	 double wt[], double res[]){
  int i,y2,yy1,yy2;
  for(i=0;i<*nn;i++){
    if(wt[i]>0){
      y2=n[i]-y[i];
      yy1=y[i]>0?y[i]:1;
      yy2=y2>0?y2:1;
      res[i]=wt[i]*(lchoose((double)n[i],(double)y[i])+
		    (s[i]-1)*n[i]*log((double)n[i])+s[i]*y[i]*log(m[i])
		    +s[i]*y2*log(1.-m[i])-(s[i]-1)*y[i]*log((double)yy1)
		    -(s[i]-1)*y2*log((double)yy2)-log(dbnc(n[i],n[i],m[i],s[i])));}
    else res[i]=0;}}

/* multiplicative binomial */
static double mbnc(int yy, int n, double m, double s){
  int y;
  double r;
  r=0.;
  for(y=0;y<=yy;y++)r+=exp(lchoose((double)n,(double)y)+(n-y)*log(1.-m)+y*(log(m)+(n-y)*y*s));
  return(r);}

void pmb(int q[], int n[], double m[], double s[], int *nn, double res[]){
  int i;
  double ss;
  for(i=0;i<*nn;i++){
    ss=log(s[i]);
    res[i]=mbnc(q[i],n[i],m[i],ss)/mbnc(n[i],n[i],m[i],ss);}}

void dmb(int y[], int n[], double m[], double s[], int *nn,
	 double wt[], double res[]){
  int i;
  double ss;
  for(i=0;i<*nn;i++){
    if(wt[i]>0){
      ss=log(s[i]);
      res[i]=wt[i]*(lchoose((double)(n[i]),(double)y[i])+y[i]*log(m[i])
		 +(n[i]-y[i])*(log(1.-m[i])
			       +(n[i]-y[i])*y[i]*ss)-log(mbnc(n[i],n[i],m[i],ss)));}
    else res[i]=0;}}

/* romberg integration routines */
static void interp(double x[], double fx[], int pts, double tab1[],
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

static void evalfn(void fcn(), double a[], double b[], int n,
		   int len, double sum[], double tmpsum[], double zz[],
		   double pnt1[], double pnt2[], double arg1[],
		   double arg2[], double arg3[], double x[])
{
  int i,j,k,nn;

  if(n==1){
    for(k=0;k<len;k++)zz[k]=0.5*(a[k]+b[k]);
    fcn(zz,arg1,arg2,arg3,len,x);
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
      fcn(zz,arg1,arg2,arg3,len,x);
      for(k=0;k<len;k++){
	tmpsum[k]+=x[k];
	zz[k]+=pnt2[k];}
      fcn(zz,arg1,arg2,arg3,len,x);
      for(k=0;k<len;k++){
	tmpsum[k]+=x[k];
	zz[k]+=pnt1[k];}}
    for(k=0;k<len;k++)sum[k]=(sum[k]+(b[k]-a[k])*tmpsum[k]/nn)/3.0;
    return;}}

static void romberg2(void fcn(), double *a, double *b, int len,
		    double *arg1, double *arg2, double *arg3, double eps,
		    int pts, int max, int *err, double sumlen[])
{
  int i,j,j1,finish;
  double errsum,*tab1,*tab2,*x,*fx,*sum,*tmpsum,*zz,*pnt1,*pnt2,*y;

  x=(double*)R_alloc((size_t)(max*len),sizeof(double));
  fx=(double*)R_alloc((size_t)(max*len),sizeof(double));
  sum=(double*)R_alloc((size_t)len,sizeof(double));
  tmpsum=(double*)R_alloc((size_t)len,sizeof(double));
  zz=(double*)R_alloc((size_t)len,sizeof(double));
  pnt1=(double*)R_alloc((size_t)len,sizeof(double));
  pnt2=(double*)R_alloc((size_t)len,sizeof(double));
  tab1=(double*)R_alloc((size_t)pts,sizeof(double));
  tab2=(double*)R_alloc((size_t)pts,sizeof(double));
  y=(double*)R_alloc((size_t)len,sizeof(double));
  if(!x||!fx||!sum||!tmpsum||!zz||!pnt1||!pnt2||!tab1||!tab2||!y){
    *err=1;
    return;}
  *err=0;
  for(i=0;i<len;i++)x[i*max]=1.0;
  for(j=0;j<max;j++){
    j1=j+1;
    evalfn(fcn,a,b,j1,len,sum,tmpsum,zz,pnt1,pnt2,arg1,arg2,arg3,y);
    finish=(j1>=pts?1:0);
    for(i=0;i<len;i++){
      fx[j+i*max]=sum[i];
      if(j1>=pts){
	interp(&x[j1-pts+i*max],&fx[j1-pts+i*max],pts,tab1,tab2,&sumlen[i],&errsum,err);
	if(*err)return;
	if(fabs(errsum)>eps*fabs(sumlen[i]))finish=0;}
      x[j1+i*max]=x[j+i*max]/9.0;
      fx[j1+i*max]=fx[j+i*max];}
    if(finish)return;}
  *err=3;
  return;}

/* simplex distribution */
static void dsimplex(double y[], double m[], double s[], double f[], int len,
	   double res[]){
  int i;
  for(i=0;i<len;i++)res[i]=exp(-pow((y[i]-m[i])/(m[i]*(1-m[i])),2)/(2*y[i]*(1-y[i])*s[i]))/sqrt(2*M_PI*s[i]*pow(y[i]*(1-y[i]),3));}

void psimplex_c(double y[], double m[], double s[], double f[], int *len,
	   double *eps, int *pts, int *max, int *err, double res[]){
  double *x;
  int i;
  x=(double*)R_alloc((size_t)(*len),sizeof(double));
  for(i=0;i<*len;i++)x[i]=0;
  romberg2(dsimplex, x, y, *len, m, s, f, *eps, *pts, *max, err, res);}

/* generalized inverse Gaussian distribution */
static void dginvgauss(double y[], double m[], double s[], double f[], int len,
	   double res[]){
  int i;
  for(i=0;i<len;i++)res[i]=pow(y[i],f[i]-1)*exp(-(1./y[i]+y[i]/pow(m[i],2))/(2.*s[i]))/(pow(m[i],f[i])*(2.*bessel_k(1/(s[i]*m[i]),fabs(f[i]),1.0)));}

void pginvgauss_c(double y[], double m[], double s[], double f[], int *len,
	   double *eps, int *pts, int *max, int *err, double res[]){
  double *x;
  int i;
  x=(double*)R_alloc((size_t)(*len),sizeof(double));
  for(i=0;i<*len;i++)x[i]=0;
  romberg2(dginvgauss, x, y, *len, m, s, f, *eps, *pts, *max, err, res);}

/* power exponential distribution */
static void dpowexp(double y[], double m[], double s[], double f[], int len,
	   double res[]){
  int i;
  double b,ss;
  for(i=0;i<len;i++){
    ss=sqrt(s[i]);
    b=1.+1./(2.*f[i]);
    res[i]=exp(-pow(fabs(y[i]-m[i])/ss,2.*f[i])/2.)/(ss*gammafn(b)*pow(2.,b));}}

void ppowexp_c(double y[], double m[], double s[], double f[], int *len,
	   double *eps, int *pts, int *max, int *err, double res[]){
  double *x;
  int i;
  x=(double*)R_alloc((size_t)(*len),sizeof(double));
  for(i=0;i<*len;i++)x[i]=fabs(y[i]-m[i])+m[i];
  romberg2(dpowexp, m, x, *len, m, s, f, *eps, *pts, *max, err, res);}
