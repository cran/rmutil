#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stddef.h>


/* polynomial interpolation algorithm for function with values fx
 evaluated at x */
void interp_sexp(double x[], double fx[], int pts, double tab1[],
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



static void evalRfn_sexp(SEXP fcn, double a[], double b[], int n, int len,
		double sum[], double tmpsum[], double zz[],
		  double pnt1[], double pnt2[], SEXP envir)
{
  //double *x, nn;
  double nn;
  int i,j,k;


  /* try a state2 def and put in c(1,2,3,9) and see if you can get it*/
  // SEXP abcd2;
  // abcd2 = PROTECT(Rf_allocVector(REALSXP, len));
  // REAL(abcd2)[0] = 1;
  // REAL(abcd2)[1] = 2;
  // REAL(abcd2)[2] = 2.5;
  // REAL(abcd2)[3] = 3.11;
  /*...*/



  SEXP call2, result2, state2, x;
  state2 = Rf_protect(Rf_allocVector(REALSXP, len));




  if(n==1){
    /*one trapezoid*/
    for(k=0;k<len;k++) zz[k]=0.5*(a[k]+b[k]);

    for(k=0;k<len;k++) REAL(state2)[k]=zz[k];      /*.C -> .Call requires this*/
    Rf_protect(call2=Rf_lang2(fcn,state2));        /*.C -> .Call requires this*/
    Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
    Rf_protect(x=Rf_coerceVector(result2,REALSXP));/*.C -> .Call requires this*/

    for(k=0;k<len;k++) sum[k]=(b[k]-a[k])*REAL(x)[k]; /*x[k] -> REAL(x)[k]*/
    Rf_unprotect(4);
    return;
  }else{
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

      for(k=0;k<len;k++) REAL(state2)[k]=zz[k];    /*.C -> .Call requires this*/
    Rf_protect(call2=Rf_lang2(fcn,state2));        /*.C -> .Call requires this*/
    Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
    Rf_protect(x=Rf_coerceVector(result2,REALSXP));/*.C -> .Call requires this*/

      for(k=0;k<len;k++){
        tmpsum[k]+=REAL(x)[k]; /*x[k] -> REAL(x)[k]*/
        zz[k]+=pnt2[k];}
      Rf_unprotect(3);

    for(k=0;k<len;k++) REAL(state2)[k]=zz[k];      /*.C -> .Call requires this*/
    Rf_protect(call2=Rf_lang2(fcn,state2));        /*.C -> .Call requires this*/
    Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
    Rf_protect(x=Rf_coerceVector(result2,REALSXP));/*.C -> .Call requires this*/


      for(k=0;k<len;k++){
        tmpsum[k]+=REAL(x)[k]; /*x[k] -> REAL(x)[k]*/
        zz[k]+=pnt1[k];}
      Rf_unprotect(3);
      }

    Rf_unprotect(1);

    /* calculate total area */
    for(k=0;k<len;k++)sum[k]=(sum[k]+(b[k]-a[k])*tmpsum[k]/nn)/3.0;
    return;}}






SEXP romberg_sexp(SEXP fcn, SEXP a, SEXP b, SEXP len, SEXP eps,
                  SEXP pts, SEXP max, SEXP err, SEXP envir)
{

  int i, j, j1, finish;
  double errsum, *tab1, *tab2, *x, *fx, *sum, *tmpsum, *zz, *pnt1, *pnt2;

  int *PTS = INTEGER(pts);
  int *MAX = INTEGER(max);
  int *LEN = INTEGER(len);
  int *ERR = INTEGER(err);
  double *EPS = REAL(eps);
  double *A = REAL(a);
  double *B = REAL(b);
  double sumlen[*LEN];


  /*allocate vectors according to member of points and length of vector*/
  tab1 = (double*)R_alloc((size_t)(*PTS)         ,sizeof(double));
  tab2 = (double*)R_alloc((size_t)(*PTS)         ,sizeof(double));
  x = (double*)R_alloc((size_t)(*MAX*(*LEN+1)),sizeof(double));
  fx = (double*)R_alloc((size_t)(*MAX*(*LEN+1)),sizeof(double));
  sum = (double*)R_alloc((size_t)(*LEN)         ,sizeof(double));
  tmpsum = (double*)R_alloc((size_t)(*LEN)         ,sizeof(double));
  zz = (double*)R_alloc((size_t)(*LEN)         ,sizeof(double));
  pnt1 = (double*)R_alloc((size_t)(*LEN)         ,sizeof(double));
  pnt2 = (double*)R_alloc((size_t)(*LEN)         ,sizeof(double));


  SEXP ans;
  Rf_protect(ans = Rf_allocVector(REALSXP, *LEN));

  if(!tab1||!tab2||!x||!fx||!sum||!tmpsum||!zz||!pnt1||!pnt2){
    *ERR=1;
    Rf_error("*ERR is now 1 -- Line 162 -- Unable to allocate memory in romberg integration C code");}
  *ERR=0;
  for(i=0;i<*LEN;i++)x[i**MAX]=1.0;

  /* iterate, decreasing step size, until convergence or max number of steps */
  for(j=0;j<*MAX;j++){
    j1=j+1;
    /* evaluate R function and calculate sum of trapezoids */
    evalRfn_sexp(fcn,A,B,j1,*LEN,sum,tmpsum,zz,pnt1,pnt2,envir);
    finish=(j1>=*PTS?1:0);
    /*repeatedly call polynomial interpolation routine*/
    for(i=0;i<*LEN;i++){
      fx[j+i**MAX]=sum[i];
      if(j1>=*PTS){
        interp_sexp(&x[j1-*PTS+i**MAX],&fx[j1-*PTS+i**MAX],*PTS,tab1,tab2,&sumlen[i],&errsum,ERR);
        if(*ERR)Rf_error("*ERR is now 2 -- Line 177 -- Division by zero in romberg integration C code");
        /*  check convergence  */
        if(fabs(errsum)>*EPS*fabs(sumlen[i]))finish=0;}
        /* decrease step size */
        x[j1+i**MAX]=x[j+i**MAX]/9.0;
       fx[j1+i**MAX]=fx[j+i**MAX];}



    if(finish){
      for(i=0;i<*LEN;i++) REAL(ans)[i] = sumlen[i];
      Rf_unprotect(1);
      return ans;}}
  *ERR=3;
  if(*ERR)Rf_error("*ERR is now 3 -- Line 191 -- No convergence in romberg integration C code");
  for(i=0;i<*LEN;i++) REAL(ans)[i] = sumlen[i];
  Rf_unprotect(1);
  return ans;}


  // SEXP call, result;
  // /* try a state2 def and put in c(1,2,3,9) and see if you can get it*/
  // SEXP abcd;
  // abcd = PROTECT(Rf_allocVector(REALSXP, 4));
  // REAL(abcd)[0] = 1;
  // REAL(abcd)[1] = 2;
  // REAL(abcd)[2] = 2.5;
  // REAL(abcd)[3] = 3.11;
  // /*...*/
  //
  // PROTECT(call   = Rf_lang2(fcn,abcd));
  // PROTECT(result = Rf_eval(call,envir));
  // SEXP foo;
  // PROTECT(foo = Rf_coerceVector(result, REALSXP));
  // int len2 = LENGTH(foo);
  // for (int i = 0; i < len2; i++)
  //   if (! R_finite(REAL(foo)[i]))
  //     Rf_error("function returned vector with non-finite element");
  //   UNPROTECT(4);
  //   return foo;
//}
