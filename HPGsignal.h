/*
 * ==========================================================================
 *
 *       Filename:  HPGsignal.h
 *
 *    Description:  header file some common time series data processing routines
 *
 *        Version:  2.0
 *        Created:  11/09/09 21:27:39
 *       Revision:  12/13/10 21:29
 *       Compiler:  gcc
 *
 *         Author:  Henri P. Gavin (hpgavin), h p gavin ~at~ duke ~dot~ e d v
 *        Company:  Duke Univ.
 *
 * ==========================================================================
 */

/*  Copyright (C) 2023 Henri P. Gavin
 
    HPGsignal is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. 
    
    HPGsignal is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with HPGsignal .  If not, see <http://www.gnu.org/licenses/>.
*/

#define MIN(u,v) ((u) < (v) ? (u) : (v))
#define MAX(u,v) ((u) > (v) ? (u) : (v))
#define SGN(u)   ((u) > 0.0 ? 1.0 : -1.0)
#define SQR(u) ((u)*(u))
#define SQR2(u,v) ((u)*(u) + (v)*(v))

#define PI 3.14159265358979323846264338327950288419716939937510

/* BDIFF : backwards differences for a vector
 */
float *bDiff ( float *u, int n, float h );

/* CDIFF : central differences for a vector
 */
float *cDiff ( float *u, int n, float h );

/* TRAPZ - integration of a vector Trapezoidal rule
 */
float trapz ( float *u, int n, float h, int n_int );
	
/* CUMTRAPZ - cumulative integration of a vector Trapezoidal rule
 */
float *cumTrapz ( float *u, int n, float h );

/* CUMTICK - cumulative sum for a vector / Tick's rule	
 */
float *cumTick ( float *u, int n, float h );

/* SMOOTH - convolution smoothing of a vector of floats
 */
float *smooth ( float *u, int n, int k_pts );

/* FIRfilter
 * the causal FIR filter is the discrete time convlution
 *  input:
 *   c  = filter coefficients  c[0], c[1], ... , c[nc-1]
 *   nc = number of filter coefficients 
 *   u  = input data u[0], u[1], ... , u[nu-1]
 *   nu = number of points in u 
 *  output:
 *   y  = filtered response y[0], y[1], ... y[nu-1].
 */
float *FIRfilter ( float *c, int nc, float *u, int nu );

/* * BANDBASS - FIR bandPass filter
 *  u = input data		
 *  n = number of points	
 *  k_pts = number of coefficients,   k_pts < 64	
 *  flo = low   cutoff frequency
 *  fhi = high  cutoff frequency
 */
float *bandPass ( float *u, int n, int k_pts, float flo, float fhi );

/* RMS - root mean square of a vector 
 */
float rms ( float *u, int n );

/* RMS2 - root mean square of a pair of vectors 
 */
float rms2 ( float *u1, float *u2, int n );

/* FIRST_PT  -  subtract the average of the first n_avg points from a vector, x
 */
void firstPoint( float *x, int points, int n_avg );

/* DEBIAS  -  subtract the average value from a vector
 */
void deBias ( float *x, int points );

/* DETREND  -  subtract a straight line from a vector u[i]  1 <= i <= n.
 *	    makes use of diagonalization via orthogonality, n: even.
 */
float *deTrend ( float *u, int n );

float *deTrendSmooth ( float *a, int n, float dt, float T );

/*  DETREND_DISPL
 * deTrend a measured acceleration record, au, without windowing or filtering
 * so that the displacement record (computed with the backwards differences and
 * the trapezoidal rule or Tick's rule) has no drift.   
 * Some padding of short accelerotrams may be required. 
 */
float *deTrendDispl( float *au, int n, float T, float dt );

/* DETREND_LEGENDRE  -  subtract a least squares polynomial line from a vector
 *  u[i]  1 <= i <= n.
 * makes use of diagonalization via orthogonality of Legendre Polynomials,
 * n: even.
 */
float *deTrendLegendre ( float *u, int n, int order );

/* DETREND_LEGENDRE_D   
 * baseLine correction of accel. to remove noise oscillations in displ.
 *
 * subtract a least squares polynomial curve from a vector
 * makes use of diagonalization via orthogonality of Legendre polynomials
 */
void deTrendLegendre_d ( float *accel, float *displ, int n, float sr, int order );

float baseLineSigmoid ( float **a, int n, float dt );

/*----------------------------------------------------------------------------
BASELINE  -  subtracts a straight line passing through first and last points
 u  = vector of data to be processed
 n  = length of the vector u
 p1 = baseLine correction results in u[p1] approx =  0
 pN = baseLine correction results in u[pN] approx = 0
 pA = number of points to average around p1 and around pN
-----------------------------------------------------------------------------*/
float *baseLine ( float *u, int n, int p1, int pN, int pA );

/*----------------------------------------------------------------------------
BASELINE_TAPER - subtracts a straight line passing through first and last points
----------------------------------------------------------------------------*/
void baseLineTaper ( float *x, int N, int points, int Nt );

/*----------------------------------------------------------------------------
END_TAPER - tapers a wave form s.t. x(0)=x'(0)=x''(0)=x'(T)=x''(T)=0, x(T)=1
          T can be 1, 2 or 5 seconds; sr is the sample rate 
-----------------------------------------------------------------------------*/
void endTaper ( float *x, int points, float T, float sr );

void cosTaper ( float *x, int points, int Ni, int Nf );

/*-----------------------------------------------------------------------------
MINMAX  -  make the maximum and minim values equal in magnitude
------------------------------------------------------------------------------*/
void peakPeak ( float *x, int points );

/*-----------------------------------------------------------------------------
MAXABSV  -  find the maximum absolute value of a vector
------------------------------------------------------------------------------*/
float maxAbsV ( float *u, int n, int *idx);

/*-----------------------------------------------------------------------------
MINABSV  -  find the minimum absolute value of a vector
------------------------------------------------------------------------------*/
float minAbsV ( float *u, int n, int *idx );

/*-----------------------------------------------------------------------------
MAXABSV2  -  find the maximum absolute value of a pair of vectors
------------------------------------------------------------------------------*/
float maxAbsV2 ( float *u1, float *u2, int n );


/* -----------------------------------------------------------------------------
 INTERPOLATE - linear 'straight-line' interpolation of a vector
 Interpolate a vector, u, of length Nu sampled at interval Du
 to a vector, y, sampled at interval Dy
 The new sample rate need not be an integer multiple of the original sample rate
----------------------------------------------------------------------------- */
void interpolate ( float *u, float *y, int Nu, float Du, float Dy, int n );

/*----------------------------------------------------------------------------- 
LIMITS - find Max, Min, Avg, RMS, T_min and T_max for a vector of floats x[1..n]
-----------------------------------------------------------------------------*/
void limits ( float *x, int fp, int lp, float dt,
		float *min, float *max, float *avg, float *rms,
		float *t_min, float *t_max );

/*  
 * AVERAGE - compute the average value of a vector from x[i] to x[i+j]
 */
float   average ( float *x, int i, int j );

/*---------------------------------------------------------------------------
RAN2  -  Random Number Generator from Numerical Recipes In C p. 212
----------------------------------------------------------------------------*/
float ran2 ( long *idum );


/*                              FILE fft.c                              */
/* Fast Fourier Transform routines from Numerical Recipes in C, Ch. 12  */
/* by Press et al. Cambridge University Press, 1988                     */
/* The sign of the exponent in the forward FFT is defined to be positve.*/
/* The sign of the exponent in the inverse FFT is defined to be negative*/


/*-----------------------------------------------------------------------------
FOUR1  -   Fast Fourier Transform 

Replaces data [1..2*nn] by its discrete Fourier transform, if isign is 1
or replaces data[1..2*nn by nn times its inverse discrete Fourier transform,
if isign is input as -1.   data is a complex array of length nn,
or equivalently, a real array of length 2*nn. 
nn MUST be an integer power of 2 (this is not checked for!)

Numerical Recipes in C, 2nd ed.  p 507.
------------------------------------------------------------------------------*/
void four1( float *data, int nn, int isign); 


/*------------------------------------------------------------------------------
TWOFFT()  -  The simultaneous FFT of two real valued functions	
		from  Numerical Recipes In C  p 415
------------------------------------------------------------------------------*/
void twofft( float *data1, float *data2, float *fft1, float *fft2, int n);


/*------------------------------------------------------------------------------
REALFT()  -   The FFT of 2n real valued discrete function points

Calculates the Fourier transform of a set of n real-valued data points.
Replaces this data (which is stored in array data[1..n]) by the
positive frequency half of its comples Fourier transform.
The real-valued first and last components of the complex transform are returned
as elements data[1] and data[2], respectively.   n must be a power of 2.
This routine also calculates the inverse transform of a complex data array
if it is the transform of real data. 
The result in this case must be multiplied by 2/n.

Numerical Recipes in C, 2nd ed, p 513.
------------------------------------------------------------------------------*/
void realft(float *data, int n, int isign);


/*------------------------------------------------------------------------------
FOURN()  -  N-dimensional FFT from Numerical Recipes in C
		nn[1..ndim] data[1..nn[1],1..nn[2],...1..nn[ndim]]
------------------------------------------------------------------------------*/
void fourn(float *data, int *nn, int ndim, int isign);


/*------------------------------------------------------------------------------
FFOUR1()  -  1-dimensional FFT from Numerical Recipes in C
           altered for consistency with original FORTRAN, tuned up
------------------------------------------------------------------------------*/
void ffour1 ( float *data, int *nn, int *isign );


/*-----------------------------------------------------------------------------
FTDSP  -  Digital Filtering and Integration Using The  F.F.T.
 y    = vector of data to be processed
 nint = number of integrations ... nint < 0 means differentation
 f_lo = low  frequency bandPass cut-off frequency
 f_hi = high frequency bandPass cut-off frequency
 sr   = sample rate ... 1 / delta_t
 n = number of data points in vector y... a power of 2
------------------------------------------------------------------------------*/
void ftdsp ( float *y, int n, float sr, float f_lo, float f_hi, int nint, float bp_fp[8] );

/* ----------------------------------------------------------------------------
SELECT_VALUE
Returns the kth smallest value in the array arr[1..n]. The input array will
be rearranged to have this value in location arr[k], with all smaller elements
moved to arr[1..k-1] (in arbitrary order) and all larger elements in
arr[k+1..n] (also in arbitrary order).
Adapted from Numerical Recipes in C, 2nd edition, 1992
----------------------------------------------------------------------------*/
float selectValue ( float *arr, unsigned long n, unsigned long k );

#undef SWAP
