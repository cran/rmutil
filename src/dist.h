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
void pdp(int q[], int *my, double m[], double s[], int *nn, double res[]);
void ddp(int y[], int *my, double m[], double s[], int *nn,
	 double wt[], double res[]);
void pmp(int q[], int *my, double m[], double s[], int *nn, double res[]);
void dmp(int y[], int *my, double m[], double s[], int *nn,
	 double wt[], double res[]);
void ppvfp(int q[], double m[], double s[], double f[], int *nn,
	 double res[]);
void dpvfp(int y[], double m[], double s[], double f[], int *nn,
	 double wt[], double res[]);
void pdb(int q[], int n[], double m[], double s[], int *nn, double res[]);
void ddb(int y[], int n[], double m[], double s[], int *nn,
	 double wt[], double res[]);
void pmb(int q[], int n[], double m[], double s[], int *nn, double res[]);
void dmb(int y[], int n[], double m[], double s[], int *nn,
	 double wt[], double res[]);

void psimplex_c(double y[], double m[], double s[], double f[], int *len,
	 double *eps, int *pts, int *max, int *err, double res[]);
void pginvgauss_c(double y[], double m[], double s[], double f[], int *len,
	 double *eps, int *pts, int *max, int *err, double res[]);
void ppowexp_c(double y[], double m[], double s[], double f[], int *len,
	 double *eps, int *pts, int *max, int *err, double res[]);
