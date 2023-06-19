/*
 * ==========================================================================
 *
 *       Filename:  HPGsignal.c
 *
 *    Description:  some common time series data processing routines
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "NRutil.h"
#include "HPGmatrix.h"
#include "HPGsignal.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define DEBUG 0
#define FTDSP_DEBUG   0		/* write intermediate.data files        */

double TEMP;
#define SWAP(a,b) {TEMP=(a);(a)=(b);(b)=TEMP;}


/* BDIFF : backwards differences for a vector
 */
float *bDiff ( float *u, int n, float h )
{
	int     i;
	float   *y;

	y = vector(1,n);	

	y[1] = (u[2] - u[1]) / h;

	for (i=2;i<=n;i++)	y[i] = (u[i] - u[i-1]) / h;

	return y;
}


/* CDIFF : central differences for a vector
 */
float *cDiff ( float *u, int n, float h )
{
	int	i;
	float	*y;

	y = vector(1,n);	

	y[1] = (u[2] - u[1]) / h;

	for (i=2;i<n;i++)	y[i] = 0.5 * (u[i+1] - u[i-1]) / h;

	y[n] = (u[n] - u[n-1]) / h;

	return y;
}


/* TRAPZ - integration of a vector using the Trapezoidal rule
 */
float trapz ( float *u, int n, float h, int n_int )
{
	int	i;
	float	s = 0.0;

	if ( n_int == 1 )	for (i=1;i<=n;i++)	s += u[i];
	if ( n_int == 2 )	for (i=1;i<=n;i++)	s += (n-i+1)*u[i];

	if ( n_int == 1 )	s *= h;
	if ( n_int == 2 )	s *= h*h;

	return s;
}


/* CUMTRAPZ - cumulative integration of a vector Trapezoidal rule
 */
float *cumTrapz ( float *u, int n, float h )
{
	int	i;
	float	*y;

	y = vector(1,n);	

	y[1] = 0.5 * u[1] * h;

	for (i=2;i<=n;i++)	y[i] = y[i-1] + 0.5 * ( u[i] + u[i-1] ) * h;

	return y;
}


/* CUMTICK - cumulative sum for a vector / Tick's rule	
 */
float *cumTick ( float *u, int n, float h  )
{
	int	i;
	float	a1 = 0.3584, a0 = 1.2832,
		*y;

	y = vector(1,n);	

	y[1] =  0.5 * u[1] * h;
	y[2] =  y[1] + 0.5 * (u[1] + u[2]) * h;

	for (i=3;i<=n;i++ )	
		y[i] = y[i-2] + ( a1*u[i] + a0*u[i-1] + a1*u[i-2] ) * h;

	return y;
}


/* SMOOTH - convolution smoothing of a vector of floats
 * http://en.wikipedia.org/wiki/Numerical_smoothing_and_differentiation 
 * smooth a vector of data via convolution with a kernel (1-cos(2*pi*[1:m]/m))/m
 */
float *smooth ( float *u, int n, int m )
{
	float	*h, sum_h=0.0, yi, *v, *y;
	int	i,j;

	if ( m % 2 == 0 ) m += 1;	// make m odd

//	fprintf(stderr," * smooth * n = %d  m = %d\n", n, m ); // debug

	h = vector(1,m);
	y = vector(1,n);
	v = vector(1-m, n+m);

 	for (i=1-m; i<=1;   i++) v[i] = u[-i];		// pad is a reflection
	for (i=1;   i<=n;   i++) v[i] = u[i];
	for (i=n;   i<=n+m; i++) v[i] = u[2*n-i];	// pad is a reflection

	if ( m < 1 )	return u;

	// h = [ -3 12 17 12 -3 ]/35;		% FIR filter coefficients
	// h = [ 0.1 0.25 0.3 0.25 0.1 ];	% FIR filter coefficients

	for (j=1;j<=m;j++)	h[j] = (1.0 + cos(PI*j/m)) / m;
	for (j=1;j<=m;j++)	sum_h += h[j];
	for (j=1;j<=m;j++)      h[j] *= (1.0-1.0/m) / ( 2.0 * sum_h );

	for ( i=1; i<=n; i++ ) {	// smooth the time series
		yi = v[i] / m;
		for (j=1; j<=m; j++)	yi += h[j]*( v[i-j] + v[i+j] );
		y[i] = yi;
	}

//	for (k=0;k<=k_pts;k++) printf("h[%2d] = %8.6f\n", k, h[k] ); // debug

	free_vector( h, 1,m );
	free_vector( v, 1-m, n+m );

	return y;
}

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
float *FIRfilter ( float *c, int nc, float *u, int nu )
{
	float	*y, yi;
	int	i,k,nk;

	y = vector(0,nu-1);	
	for ( i=0; i<nu; i++ )	y[i] = 0.0;

	for ( i=0; i<nu; i++ ) {

		yi = 0;
		nk = MIN(i+1,nc);
		for (k=0; k < nk; k++)    yi += c[k]*u[i-k];
		//yi = c[0]*u[i];
		//nk = MIN( MIN(i+1,nc), MIN(nu-1-i,nc) );
		//for (k=1; k < nk; k++)    yi += c[k]*(u[i-k] + u[i+k]);
		y[i] = yi;
	}
//	for (k=0;k<=k_pts;k++)      printf("h[%2d] = %8.6f\n", k, h[k] );

	return y;
}


/* BANDBASS - symmetric FIR bandPass filter
 *  u = input data    u[1], u[2], ... , u[n-1]
 *  n = number of points  
 *  k_pts = number of coefficients,   k_pts < 64	
 *  flo = low   cutoff frequency
 *  fhi = high  cutoff frequency
 */
#define NH 64
float *bandPass ( float *u, int n, int k_pts, float flo, float fhi )
{
	float	h[NH], *y, yi;
	int	i,k;

	if ( k_pts >= NH ) fprintf(stderr," * bandPass * k_pts > %d\n", NH );

	y = vector(1,n);	
	for (i=1;i<=n;i++)	y[i] = 0.0;

	h[0] = ( fhi - flo ) ;
	for (k=1; k <= k_pts; k++)
		h[k] = ( sin(PI*fhi*k) - sin(PI*flo*k) ) / (PI*k);
        for ( i=1; i<=n; i++ ) {
		if (k_pts<i && i<n-k_pts) {
			yi = h[0]*u[i];
			for (k=1; k < k_pts; k++)
				yi += h[k]*( u[i-k] + u[i+k] );
		} else	yi = u[i];
		y[i] = yi;
	}
//	for (k=0;k<=k_pts;k++)      printf("h[%2d] = %8.6f\n", k, h[k] );

	return y;
}
#undef NH


/* RMS - root mean square of the first n points of a vector 
 */
float rms ( float *u, int n )
{
	int	i;
	float	avg = 0.0,
		rms = 0.0;

	for (i=1;i<=n;i++)	avg += u[i];
	avg /= (float) n;
	for (i=1;i<=n;i++)	rms += (u[i]-avg)*(u[i]-avg);

	rms = sqrt( rms / (float) n );

	return rms;
}


/* RMS2 - root mean square of a pair of vectors 
 */
float rms2 ( float *u1, float *u2, int n )
{
	int	i;
	float	rms=0.0;

	for (i=1;i<=n;i++ )	rms += u1[i]*u1[i] + u2[i]*u2[i];

	rms = sqrt(rms/n);

	return rms;
}


/* FIRST_PT  -  subtract the average of the first n_avg points from a vector
 */
void firstPoint( float *x, int points, int n_avg )
{
        float   x_avg = 0.0;
        int     p;
        
	for (p=1; p <= n_avg; p++)	x_avg += x[p];
	x_avg /= (float) n_avg;
	for (p=1; p <= points; p++)	x[p] -= x_avg;

        return;
}


/* DEBIAS  -  subtract the average value from a vector
 */
void deBias ( float *x, int points )
{
        float   Avg_X = 0.0;            /*  average value */
        int     j;

        for (j=1; j <= points; j++)     Avg_X  += x[j];
        Avg_X /= (float) points;
        for (j=1; j <= points; j++)     x[j] -= Avg_X;

        return;
}


/* DETREND  
 * subtract a straight line from a vector u[i]  1 <= i <= n.
 * makes use of diagonalization via orthogonality, n: even.
 */ 
float *deTrend ( float *u, int n )
{
        float   i, Sum_U = 0.0, Sum_Ui = 0.0, Sum_ii = 0.0, 
                a, b;                           /* slope and intercept  */
        int     p;


	for (p=1; p <= n; p++) {
		i = ((p-n/2.) <= 0) ? p - n/2. - 1.0 : p - n/2.0;
		Sum_U  += u[p];
		Sum_Ui += u[p]*i;
		Sum_ii += i*i;
	}

	a = Sum_Ui / Sum_ii;
	b = Sum_U / (float) n;                     /* average value */

	for (p=1; p <= n; p++) {
		i = ((p-n/2) <= 0) ? p - n/2.0 - 1.0 : p - n/2.0;
		u[p] -= ( a*i + b);
	}

        return u;
}


/* DETREND_SMOOTH
 * detrend a measured acceleration record, a, without windowing or filtering
 * by subtracting a very smoothed (period = 50 sec or more) version of itself
 * Some padding of short accelerotrams may be required. 
 */
float *deTrendSmooth ( float *a, int n, float dt, float T )
{
        float   *as;	/*  smoothed version of a */

        int     i, p;

	as = vector(1,n);

	p = MIN ( floor(T/dt),floor(n/4) );	// number of points in smoothing

	as = smooth ( a, n, p );		// compute trend in accel

	for ( i = 1; i <= n; i++ )	a[i] -= as[i];	// remove trend in displ

	free_vector( as, 1,n );

	return a;
}


/* DETREND_DISPL
 * deTrend a measured acceleration record, a, without windowing or filtering
 * so that the displacement record (computed with the central differences and
 * the trapezoidal rule or Tick's rule) has no drift.   
 * Some extra padding of short accelerotrams may be required. 
 */
float *deTrendDispl ( float *a, int n, float dt, float T )
{
        float   *d,	/*  displacement of uncorr. and corr. accel.	*/
		*ds,	/*  smoothed displacement of uncorrected acceleration */
		*v;	/*  velocity of corrected acceleration		*/

        int     i, p;

	d  = vector(1,n);
	ds = vector(1,n);
	v  = vector(1,n);

	p = MIN ( floor(T/dt),floor(n/4) );	// number of points in smoothing

	a  = deTrendLegendre ( a, n, 2 );	// remove quadratic trend
	a  = baseLine ( a, n, 1, n, 2 );

	d  = cumTrapz ( cumTrapz ( a, n, dt ), n, dt );	// integrate accel to displ
	ds = smooth ( d, n, p );			// compute trend in displ
	for ( i = 1; i<=n; i++ )	d[i] -= ds[i];	// remove trend in displ
	v  = cDiff ( d, n, dt );			// diff. displ to veloc
	for ( i = n-p; i<=n; i++ )	v[i] = 0.0;   // cancel smoothing end-effect
	a  = cDiff ( v, n, dt );			// diff. veloc to accel

//	d  = cumTick ( cumTick ( a, n, dt ), n, dt );	// integrate accel to displ
//	d  = baseLine ( d, n, 1, n, 20 );

	free_vector(  d, 1,n );
	free_vector( ds, 1,n );
	free_vector(  v, 1,n );

	return a;
}


/* DETREND_LEGENDRE  
 * subtract a least squares polynomial line from a vector
 *			u[i]  1 <= i <= n.
 * makes use of diagonalization via orthogonality of Legendre Polynomials,
 * n: even.
 */
float *deTrendLegendre ( float *u, int n, int order )
{
        float   x, *Pm, *Pm_1, *Pm_2, Sum_uPm = 0.0, c[16];
        int     p, m;

	Pm   = vector(1,n);			/*  Legendre polynomials */
	Pm_1 = vector(1,n);			/*  Legendre polynomials */
	Pm_2 = vector(1,n);			/*  Legendre polynomials */


	for ( m=0; m <= order; m++) {
	    Sum_uPm = 0.0;
	    for (p=1; p <= n; p++) {
		x = ((p-n/2.) <= 0) ? p - n/2. - 1.0 : p - n/2.0;
		x *= 2.0/n;
		if ( m == 0 ) {	Pm[p] = 1.0; Pm_1[p] = 0.0; }
		if ( m == 1 ) {	Pm[p] =  x ; Pm_1[p] = 1.0; }
		if ( m > 1 ) 
			Pm[p] = (2.0-1.0/m)*x*Pm_1[p] - (1.0-1.0/m)*Pm_2[p];

		Pm_2[p] = Pm_1[p];
		Pm_1[p] = Pm[p];

		Sum_uPm += u[p] * Pm[p];

	    }
	    c[m] = (2.0*m+1.0) / (n-1.0) * Sum_uPm;

	}

	for ( m=0; m <= order; m++) {
	    for (p=1; p <= n; p++) {
		x = ((p-n/2.) <= 0) ? p - n/2. - 1.0 : p - n/2.0;
		x *= 2.0/n;
		if ( m == 0 ) {	Pm[p] = 1.0; Pm_1[p] = 0.0; };
		if ( m == 1 ) {	Pm[p] =  x ; Pm_1[p] = 1.0; }
		if ( m > 1 ) 
			Pm[p] = (2.0-1.0/m)*x*Pm_1[p] - (1.0-1.0/m)*Pm_2[p];

		Pm_2[p] = Pm_1[p];
		Pm_1[p] = Pm[p];

		u[p] -= c[m]*Pm[p];
	    }
        }

	free_vector(Pm,1,n);
	free_vector(Pm_1,1,n);
	free_vector(Pm_2,1,n);

        return u;
}


/* DETREND_LEGENDRE_D   
 * baseLine correction of accel. to remove noise oscillations in displ.
 * accerlation may be processed to correspond to a residual displacement 
 * 
 * subtract a least squares polynomial line from a vector
 * makes use of diagonalization via orthogonality of Legendre Polynomials,
 * n: even.
 */
void deTrendLegendre_d ( float *accel, float *displ, int n, float sr,int order )
{
	float   *x, **P, **Pp, **Ppp, **PL, **b, Sum_dPk = 0.0, dx_dt;
	int	Nc = 3,		// number of curve-fit constraints
		nPL = order+1+Nc,
		i,j, p, k;

	x    = vector(1,n);	//  time-scaled independent variable, -1<=x<=1 
	P    = matrix(0,order,1,n);	//  Legendre polynomials 
	Pp   = matrix(0,order,1,n);	//  d/dx P(x)  , 1st derivitive of P(x)
	Ppp  = matrix(0,order,1,n);	//  d/dx Pp(x) , 2nd derivitive of P(x)

	PL     = matrix(1,nPL,1,order+1+Nc);	//  system mtx for cstr lst sqr
	b      = matrix(1,nPL,1,1);
	for (i=1; i<=nPL; i++) b[i][1] = 0.0;
	for (i=1; i<=nPL; i++) for (j=1; j<=nPL; j++)	PL[i][j] = 0.0;

	// Legendre polynomial basis is defined for -1 <= x <= +1
	for (p=1; p <= n; p++) {
		x[p] = ((p-n/2.0) <= 0) ? p - n/2.0 - 1.0 : p - n/2.0;
		x[p] *= 2.0 / n;
	}

	// evaluate Legendre polynomials and its first two derivitives
	Legendre( order, x, n, P, Pp, Ppp );


	// fit the displacement with Legendre polynomials
	// compute terms for the system matrix, including constraints
	for ( k=0; k <= order; k++) {

		Sum_dPk = 0.0;
		for (p=1; p <= n; p++)	Sum_dPk += displ[p] * P[k][p];
		b[k+1][1] = 4.0 / (n-1.0) * Sum_dPk;

		// unconstrained solution ...
		// c[k] = (2.0*k+1.0) / (n-1.0) * Sum_dPk;

		PL[k+1][k+1] = 4.0 / (2.0*k+1.0 );

		// at t = 0;  (use Nc = 3)
		PL[k+1][order+1+1]  =  PL[order+1+1][k+1]  =  P[k][1];	//d(0)=0
		PL[k+1][order+1+2]  =  PL[order+1+2][k+1]  =  Pp[k][1];	//v(0)=0
		PL[k+1][order+1+3]  =  PL[order+1+3][k+1]  =  Ppp[k][1];//a(0)=0

		// at t = T;  (use Nc = 6)
//		PL[k+1][order+1+4]  =  PL[order+1+4][k+1]  =  P[k][n];	//d(T)=0
//		PL[k+1][order+1+5]  =  PL[order+1+5][k+1]  =  Pp[k][n]; //v(T)=0
//		PL[k+1][order+1+6]  =  PL[order+1+6][k+1]  =  Ppp[k][n];//a(T)=0

	}

#if DEBUG
	save_matrix( "P.dbg",   P,   0,order, 1,n,   1, "w" );
	save_matrix( "Pp.dbg",  Pp,  0,order, 1,n,   1, "w" );
	save_matrix( "Ppp.dbg", Ppp, 0,order, 1,n,   1, "w" );
	save_matrix( "PL.dbg",  PL,  1,nPL ,  1,nPL, 0, "w" );
	save_matrix( "b.dbg",   b,   1,nPL ,  1,1  , 0, "w" );
	save_vector( "displ.dbg",   displ,   1,n, "w" );
	save_vector( "accel.dbg",   accel,   1,n, "w" );
	save_vector( "x.dbg",       x,       1,n, "w" );
#endif

	// find the curve-fit coefficients and the Lagrange multipliers
	gaussj ( PL, nPL, b, 1 );

	// subtract the acceleration of the fit displacement from accel..dbga
	dx_dt = 4.0 * sr * sr / (n*n);
	for (p=1; p <= n; p++) 
		for ( k=0; k <= order; k++) 
			accel[p] -= b[k+1][1]*Ppp[k][p] * dx_dt;

#if DEBUG
	save_matrix( "c.dbg", b, 1,nPL, 1,1,  0, "w" );
	save_vector( "accel1.dbg", accel, 1,n, "w" );
#endif

	free_vector(x,1,n);

	free_matrix(P,0,order,1,n);
	free_matrix(Pp,0,order,1,n);
	free_matrix(Ppp,0,order,1,n);

	free_matrix(PL,1,nPL,1,nPL);
	free_matrix(b,1,nPL,1,1);

	return;
}


/* BASELINE  
 * subtracts a straight line passing through first and last points
 * u  = vector of data to be processed
 * n  = length of the vector u
 * p1 = baseLine correction results in u[p1] approx =  0
 * pN = baseLine correction results in u[pN] approx = 0
 * pA = number of points to average around p1 and around pN
 */
float *baseLine ( float *u, int n, int p1, int pN, int pA )
{
	float	u1 = 0.0, uN = 0.0, slope ;
	int	p;

	for (p=p1; p<=p1+pA; p++)	u1 += u[p];
	u1 /= (float) (pA+1.0);
	for (p=pN; p>=pN-pA; p--)	uN += u[p];
	uN /= (float) (pA+1.0);

	slope = (float)( uN - u1 ) / (float)( pN - p1 - pA );

	for (p=p1; p<=pN; p++)	u[p] -= ( u1 + slope * (p - p1 - pA/2.0) );
	for (p=pN; p<=n;  p++)	u[p] = 0.0;

	return u;
}


/* BASELINE_TAPER 
 * subtracts a straight line passing through first and last points
 */
void baseLineTaper ( float *u, int points, int N, int Nt )
{
        int     p;
        float   u1, uP;

        u1 = u[1+Nt];      uP = u[points-Nt];
        for (p=1; p<=N; p++)   u[p] -= ( u1 + (uP-u1)*(p-Nt-1)/(points-1) );
        return;
}


/* BASELINE_SIGMOID 
 * subtracts a sigmoid to compute quake
 * ground velocity from bi-axial earthquake ground acceleration 
 */
float baseLineSigmoid ( float **accel, int n, float dt )
{
	float	CAV=0.0, CAV_T=0.0,	// cumulative absolute velocity
		*ac, *vc, *dc,		// acceleration correction vectors
		vT[4], dT[4],
		t0, a, b, z,	// initial time, and envelope parameters
		ts, exp_ts;	// scaled time and exponential of scaled time

	int	iter, k, p,
		p002, p050, p250, p500, p750, p950;


	ac = vector(1,n);
	vc = vector(1,n);
	dc = vector(1,n);

	// first, deTrend from first point to last point in each accel record
	for (k=1; k<=2; k++)	accel[k] = deTrend( accel[k], n );

	// determine duration parameters for the the record
	p002 = p050 = p250 = p500 = p750 = p950 = 1;

	for ( p=1; p <= n; p++ )
		CAV_T += sqrt( SQR(accel[1][p]) + SQR(accel[2][p]) ) * dt;
	for ( p=1; p <= n; p++ ) {
		CAV   += sqrt(SQR(accel[1][p])+SQR(accel[2][p])) * dt / CAV_T; 
		if ( CAV < 0.002 ) p002 = p;
       		if ( CAV < 0.050 ) p050 = p;
		if ( CAV < 0.250 ) p250 = p;
		if ( CAV < 0.500 ) p500 = p;
		if ( CAV < 0.750 ) p750 = p;
		if ( CAV < 0.950 ) p950 = p;
	}       
#if DEBUG
	printf("p002=%d p050=%d p250=%d p500=%d p750=%d p950=%d p1000=%d\n",
		p002, p050, p250, p500, p750, p950, n );
#endif
                
	t0 = p002 * dt;			// initial time of no acceleration
	b  = 0.5*(p750-p500) * dt;	// decay-time constant
	a  = (p500-p002) * dt / b;	// rise-time exponent
	z  = 0.5;
        
	// set of displacement and velocity correction records
	for ( p=1; p<=n; p++ ) {
		ts = (p*dt- a*b - t0) / (z*b);
		exp_ts = exp(-ts);
		ac[p] = 1/(1+exp_ts);		// not used
		vc[p] = exp_ts * pow(1+exp_ts,-2.0) / (z*b);
		dc[p] = exp_ts * pow(1+exp_ts,-3.0) * (exp_ts-1) / (z*b*z*b);
	}       

	// apply the sigmoidal baseLine correction iteratively
	for (iter=1; iter<=10; iter++) {
		for ( k=1; k<=2; k++ ) {
			vT[k] = trapz( accel[k], dt, n, 1 );
			dT[k] = trapz( accel[k], dt ,n, 2 );
			for ( p=1; p<=n; p++ ) 
				accel[k][p] -= (vT[k]*vc[p]+dT[k]*dc[p]);
		}
#if DEBUG
		printf("iter %2d:  vT_1= %8.3f, vT_2= %8.3f,  dT_1= %8.3f,  dT_2= %8.3f\n",
					iter, vT[1], vT[2], dT[1], dT[2] );
#endif
		if ( fabs(vT[1]) < 0.1 && fabs(dT[1]) < 2.0 &&
		     fabs(vT[2]) < 0.1 && fabs(dT[2]) < 2.0 )
			break;
	}

	free_vector(ac,1,n);
	free_vector(vc,1,n);
	free_vector(dc,1,n);

	return CAV_T;

}


/* END_TAPER 
 * taper points (1) to (N) and points (points-N) to (points) of a wave form, x,
 * such that. x(0)=x'(0)=x''(0)=x'(T)=x''(T)=0, x(T)=1
 * The duration of the taper, T, can be 1, 2 or 5 seconds; sr is the sample rate 
 */
void endTaper ( float *x, int points, float sr, float T )
{
        int     p, N;
	float	t, a3,a4,a5, taper;

	if (T==1.0) { a3=10.0;      a4= -15.0;      a5=6.0;      }
	if (T==2.0) { a3= 2.5;      a4=  -1.875;    a5=0.375;    }
	if (T==5.0) { a3= 1.381215; a4=  -0.414365; a5=0.033149; }
        
	N = (int)(T*sr);
	for (p=1; p<=N; p++) {
		t = p/sr;
		taper       = a3*pow(t,3.0) + a4*pow(t,4.0) + a5*pow(t,5.0); 
		x[p]        = x[p]*taper;
		x[points-p] = x[points-p]*taper;
	}
}

/* COS_TAPER 
 * taper points (1) to (Ni) and points (points-Nf) to (points) of a signal, x,
 * by multiplying by a cosine taper function. 
 */
void cosTaper ( float *x, int points, int Ni, int Nf )
{
        int     p; 
    
        for (p = 0; p<=Ni; p++) x[p+1]      *= 0.5*(1-cos(PI*p/(Ni+1)));
        for (p = 0; p<=Nf; p++) x[points-p] *= 0.5*(1+cos(PI*p/(Nf+1)));
}


/* PEAK_PEAK 
 * make the maximum and minim values equal in magnitude
 */
void peakPeak ( float *x, int points )
{
        float   Min_X = 0.0, Max_X = 0.0,       /* min and max values */
                Avg_X = 0.0;                    /* avg of min and max */
        int     j;

        for (j=1; j <= points; j++)     {
		Min_X = MIN ( Min_X, x[j] );
		Max_X = MAX ( Max_X, x[j] );
//                if ( Min_X > x[j] )     Min_X = x[j];
//                if ( Max_X < x[j] )     Max_X = x[j];
        }

        Avg_X  = ( Min_X + Max_X ) / 2.0;

        for (j=1; j <= points; j++)     x[j] -= Avg_X;

        return;
}


/* MAXABSV  -  find the maximum absolute value of the first n points of a vector
 */
float maxAbsV ( float *u, int n, int *idx)
{
	int	i;
	float	y;

	y = fabs(u[1]);
	*idx = 1;
	for (i=2;i<=n;i++) {
		if ( fabs(u[i]) > y ) {
			y = fabs(u[i]);
			*idx = i;
		}
	}
	return y;
}


/* MINABSV  -  find the minimum absolute value of a vector
 */
float minAbsV ( float *u, int n, int *idx )
{
	int	i;
	float	y;

	y = fabs(u[1]);
	*idx = 1;
	for (i=2;i<=n;i++) {
		if ( fabs(u[i]) < y ) {
			y = fabs(u[i]);
			*idx = i;
		}
	}
	return y;
}


/* MAXABSV2  -  find the maximum absolute value of a pair of vectors
 */
float maxAbsV2 ( float *u1, float *u2, int n )
{
	int	i;
	float	max_abs_v, abs_vi;

	max_abs_v = sqrt(u1[1]*u1[1] + u2[1]*u2[1]);
	for (i=2;i<=n;i++) {
		abs_vi = sqrt(u1[i]*u1[i] + u2[i]*u2[i]);
		max_abs_v = (max_abs_v > abs_vi ) ? max_abs_v : abs_vi ;
	}

	return max_abs_v;
}


/* INTERPOLATE - linear 'straight-line' interpolation of a vector
 *	Interpolate a vector, u, of length Nu sampled at interval Du
 *	to a vector, y, sampled at interval Dy
 * The new sample rate need not be an integer multiple of the
 * original sample rate
 */
void interpolate ( float *u, float *y, int Nu, float Du, float Dy, int n )
{

	int	i=1, j=1, Ny;

	float	t_u, t_y;		// original and interpolated time

	Ny = floor(Nu * Du / Dy);	// number of interpolated points 

	t_u = (float) (i-1.0) * Du;

	for (j = 1; j<=Ny; j++) {	// loop along interpolated record 

		t_y = (float) (j-1.0) * Dy; // time value in interpolated record

		while ( t_y > t_u + Du ) {
			++i;	
			// t_u is previous time value in the "u" record 
			t_u = (float) (i-1.0) * Du;
		}

		y[j+n] = u[i] + ( u[i+1]-u[i] ) * (t_y - t_u) / Du;
	}

	return;
}


/* LIMITS - find Max, Min, Avg, RMS, T_min and T_max
 * for a vector of floats x[fp..lp] 
 */
void limits (	float *x,		// vector of data
		int fp, int lp,		// first point and last point used 
		float dt,		// time increment for the vector
		float *min, float *max, float *avg, float *rms,
		float *t_min, float *t_max )
{
	int     t;
	
	*min = *max = x[1];
	*t_min = *t_max = 1;
	*avg = *rms = 0;
	for (t=fp; t<=lp; t++) {
		*avg +=  x[t];
		*rms +=  x[t]*x[t];
//		*min = MIN ( x[t], *min );
//		*max = MAX ( x[t], *max );
		if (x[t] < *min) { *min = x[t]; *t_min = t; }
		if (x[t] > *max) { *max = x[t]; *t_max = t; }
	}
	*avg /= (float)(lp-fp+1);
	*rms = sqrt ( *rms / (float) (lp-fp+1) );
	*t_min *= dt;
	*t_max *= dt;
}

/* 
 * AVERAGE - compute the average value of a vector from x[i] to x[i+j]
 */
float	average ( float *x, int i, int j )
{
	float	avg = 0.0;
	int	k;

	for (k=i; k<=i+j; k++) avg += x[k];

	return	( avg / ((float)(j)) );
}

/* RAN2  -  Random Number Generator from Numerical Recipes In C p. 212
 */
#define M 714025
#define IA 1366
#define IC 150889

float ran2 ( long *idum )
{
        static long iy,ir[98];
        static int iff=0;
        int j;

        if (*idum < 0 || iff == 0) {
                iff=1;
                if ((*idum=(IC-(*idum)) % M) < 0) *idum = -(*idum);
                for (j=1;j<=97;j++) {
                        *idum=(IA*(*idum)+IC) % M;
                        ir[j]=(*idum);
                }
                *idum=(IA*(*idum)+IC) % M;
                iy=(*idum);
        }
        j=1 + 97.0*iy/M;
        if (j > 97 || j < 1) fprintf(stderr,"RAN2: j out of range.\n");
        iy=ir[j];
        *idum=(IA*(*idum)+IC) % M;
        ir[j]=(*idum);
        return (float) iy/M;
}
#undef M
#undef IA
#undef IC


/* Fast Fourier Transform routines from Numerical Recipes in C, Ch. 12	*/
/* by Press et al. Cambridge University Press, 1988			*/
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
void four1( float *data, int nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double	wtemp,wr,wpr,wpi,wi,theta;
	float	tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {				/* Bit Reversal	*/
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep = mmax << 1;
		theta = isign*(6.28318530717959 / mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m=1; m<mmax; m+=2) {
			for (i=m; i<=n; i+=istep) {
				j = i + mmax;		/* Danielson-Lanczos */
				tempr = wr*data[j] - wi*data[j+1];
				tempi = wr*data[j+1] + wi*data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}			/* Trigonometric Recurrence */
			wr = (wtemp=wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}
	return;
}


/*------------------------------------------------------------------------------
TWOFFT()  -  The simultaneous FFT of two real valued functions	
		from  Numerical Recipes In C  p 415
------------------------------------------------------------------------------*/
void twofft( float *data1, float *data2, float *fft1, float *fft2, int n)
{
	int nn3,nn2,jj,j;
	float rep,rem,aip,aim;
	void four1();

	nn3=1+(nn2=2+n+n);
	for (j=1,jj=2;j<=n;j++,jj+=2) {
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	}
	four1(fft1,n,1);
	fft2[1]=fft1[2];
	fft1[2]=fft2[2]=0.0;
	for (j=3;j<=n+1;j+=2) {
		rep=0.5*(fft1[j]+fft1[nn2-j]);
		rem=0.5*(fft1[j]-fft1[nn2-j]);
		aip=0.5*(fft1[j+1]+fft1[nn3-j]);
		aim=0.5*(fft1[j+1]-fft1[nn3-j]);
		fft1[j]=rep;
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j] = -aim;
		fft2[j]=aip;
		fft2[j+1] = -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	}
	return;
}


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
void realft(float *data, int n, int isign)
{
	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;
	void four1();

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r =  c1*(data[i1]+data[i3]);
		h1i =  c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i =  c2*(data[i1]-data[i3]);
		data[i1] =  h1r + wr*h2r - wi*h2i;
		data[i2] =  h1i + wr*h2i + wi*h2r;
		data[i3] =  h1r - wr*h2r + wi*h2i;
		data[i4] = -h1i + wr*h2i + wi*h2r;
		wr = (wtemp=wr)*wpr - wi*wpi + wr;
		wi = wi*wpr + wtemp*wpi + wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1] = c1*((h1r=data[1])+data[2]);
		data[2] = c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}


/*------------------------------------------------------------------------------
FOURN()  -  N-dimensional FFT from Numerical Recipes in C
		nn[1..ndim] data[1..nn[1],1..nn[2],...1..nn[ndim]]
------------------------------------------------------------------------------*/
void fourn(float *data, int *nn, int ndim, int isign)
{
	int i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	int ibit,idim,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;	/* for trig. recurrences */

	ntot=1;
	for (idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];	/* total number of complex values */
	nprev=1;
	for (idim=ndim;idim>=1;idim--) { /* main loop over dimensions */
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {	/* bit reversal */
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;			/* Danielson-Lanczos	*/
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;	/* D.-L. */
						k2=k1+ifp1;
						tempr=wr*data[k2]-wi*data[k2+1];
						tempi=wr*data[k2+1]+wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}


/* altered for consistency with original FORTRAN, tuned up	*/

void ffour1 ( float *data, int *nn, int *isign )
{
    int n, mmax, m, j, i;
    double wtemp, wr, wpr, wpi, wi, theta, wpin;
    float tempr, tempi, datar, datai, data1r, data1i;
    n = *nn * 2;
    j = 0;
    for (i = 0; i < n; i += 2) {
        if (j > i) {            /* could use j>i+1 to help compiler analysis */
            SWAP(data[j], data[i]);
            SWAP(data[j + 1], data[i + 1]);
        }
        m = *nn;
        while (m >= 2 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    theta = 3.141592653589795 * .5;
    if (*isign < 0)
        theta = -theta;
    wpin = 0;                   /* sin(+-PI) */
    for (mmax = 2; n > mmax; mmax *= 2) {
        wpi = wpin;
        wpin = sin(theta);
        wpr = 1 - wpin * wpin - wpin * wpin;    /* cos(theta*2) */
        theta *= .5;
        wr = 1;
        wi = 0;
        for (m = 0; m < mmax; m += 2) {
            j = m + mmax;
            tempr = (float) wr *(data1r = data[j]);
            tempi = (float) wi *(data1i = data[j + 1]);
            for (i = m; i < n - mmax * 2; i += mmax * 2) {
                /* mixed precision not significantly more
                 * accurate here; if removing float casts,
                 * tempr and tempi should be double */
                tempr -= tempi;
                tempi = (float) wr *data1i + (float) wi *data1r;
                /* don't expect compiler to analyze j > i+1 */
                data1r = data[j + mmax * 2];
                data1i = data[j + mmax * 2 + 1];
                data[i] = (datar = data[i]) + tempr;
                data[i + 1] = (datai = data[i + 1]) + tempi;
                data[j] = datar - tempr;
                data[j + 1] = datai - tempi;
                tempr = (float) wr *data1r;
                tempi = (float) wi *data1i;
                j += mmax * 2;
            }
            tempr -= tempi;
            tempi = (float) wr *data1i + (float) wi *data1r;
            data[i] = (datar = data[i]) + tempr;
            data[i + 1] = (datai = data[i + 1]) + tempi;
            data[j] = datar - tempr;
            data[j + 1] = datai - tempi;
            wr = (wtemp = wr) * wpr - wi * wpi;
            wi = wtemp * wpi + wi * wpr;
        }
    }
}


/* FTDSP  -  Digital Filtering and Integration Using The  F.F.T.
 * y    = vector of data to be processed
 * nint = number of integrations ... nint < 0 means differentation
 * f_lo = low  frequency bandPass cut-off frequency
 * f_hi = high frequency bandPass cut-off frequency
 * sr   = sample rate ... 1 / delta_t
 * n = number of data points in vector y... n need not be a power of 2
 * bp_fp[8] = returned band pass filter parameters
 */
void ftdsp ( float y[], int n, float dt, float f_lo, float f_hi, int nint, float bp_fp[8] )
{
	float	w,		// angular frequency (rad/sec) 
		temp,
		*Y,		// vector for FFT calculations, length = 2^q
		taper = 1.0;	// taper the time- and frequency-domain

	int	N = 256,	// the lowest power of 2 greater than n
		i1,		// point of first data in vector Y
		binLo, binHi,	// low and high freq. bin numbers
		Pw,		// nummber of points in window transition
		nBandLo, nBandHi, // points in frequency-domain window taper
		i, k;

#if FTDSP_DEBUG
	time_t	now;
	FILE	*fp_debug;
#endif

	// y = deTrendSmooth ( y, n, dt, 1.0/f_lo ); // smooth baseline :)
	// y = deTrend(y,n);		// detrend or baseLine correction ??

 	// number of points in the time-domain taper
	Pw = (int) floor( MIN( 0.10/(f_lo*dt), (n/50) ));
	if ( f_lo*n*dt > 10 )  Pw = (int) floor( 0.1/(f_lo*dt) );
	for (i=1;i<=Pw;i++) {		// taper the window in the time-domain 
		taper  = 0.5 * (1.0 - cos(0.5*PI*(i-1)/Pw));
		y[i]   *= taper;
		y[n-i] *= taper;
	}

	// allocate memory for data of length, N, equal to a power of 2, N > n
	N  = (int) pow( 2.0 , ceil(log(n)/log(2.0)));
	i1 = (int) floor((N-n)/2.0);

	Y = vector( 1, N );		// working vector for FFT computations

	// place the data, y, in the center of the power-of-two vector, Y
	for (i = 1;      i <= i1; i++) Y[i] = 0.0;
	for (i = 1;      i <= n;  i++) Y[i+i1-1] = y[i];
	for (i = i1+n-1; i <= N;  i++) Y[i] = 0.0;

#if FTDSP_DEBUG
	fp_debug = fopen("ftdsp_debug_file_1.dbg","w");
	(void) time(&now);
	fprintf(fp_debug,"%% ftdsp_debug_file_1 ... TAPER(DETREND(y)) ... %s",ctime(&now));
	for (k=1; k<=N; k++) fprintf(fp_debug,"%11.4e\n", Y[k] );
	fclose(fp_debug);
	fprintf(stderr,"n = %d  N = %d  i1 = %d \n", n, N, i1 );
	fprintf(stderr," ... 1\n");
#endif

	binLo = 1 + (int) (2*N*f_lo*dt);
	binHi = 2 + (int) (2*N*f_hi*dt);
	if (!(binLo%2)) --binLo;		// binLo is even, make it odd
	if (  binHi%2)  ++binHi;		// binHi is odd,  make it even

	// width of the frequency-domain transition band
	nBandLo = nBandHi = (int) abs(binHi - binLo)/10.0;
	nBandLo = (nBandLo < binLo) ? nBandLo : binLo; // for pass, f_lo < f_hi
	nBandHi = (nBandHi < binHi) ? nBandHi : binHi; // for cut,  f_hi < f_lo

	bp_fp[0] = (float)(Pw);	// number of points in time-domain taper
	bp_fp[1] = (float)(binLo-nBandLo)*0.5/N/dt; // band-pass transition band
	bp_fp[2] = (float)(binLo+nBandLo)*0.5/N/dt; // band-pass transition band
	bp_fp[3] = (float)(binHi-nBandHi)*0.5/N/dt; // band-pass transition band
	bp_fp[4] = (float)(binHi+nBandHi)*0.5/N/dt; // band-pass transition band
	bp_fp[5] = (float)(N);	// number of points in the FFT

	realft ( Y, N, 1 );				/* forward real FFT  */

#if FTDSP_DEBUG
	fp_debug = fopen("ftdsp_debug_file_2.dbg","w");
	(void) time(&now);
	fprintf(fp_debug,"%% ftdsp_debug_file_2 ... REALFT(TAPER(DETREND(y))) ... %s\n",ctime(&now));
	fprintf(fp_debug,"%%  k   frequency   real(y[k])  imag(y[k]) \n");
	for (k=1; k<=N; k+=2) {
		w = PI*(k-1)/n/dt;
		fprintf(fp_debug,"%4d %11.3e %11.3e %11.3e\n", k, w/(2*PI), Y[k], Y[k+1] ) ;
	}
	fclose(fp_debug);

	fprintf(stderr,"binLo = %d  binHi = %d  nBandLo = %d  nBandHi = %d nint = %d\n",
				binLo, binHi, nBandLo, nBandHi, nint );
	fprintf(stderr," ... 2\n");
#endif


	for (k=3; k<=N; k+=2) {	/*  multiply by the integration window */

		w = PI*(k-1)/n/dt;

		for ( i= 1; i <= nint; i++ ) {		/* integrate */
			temp	=  Y[k];
			Y[k]	= -Y[k+1]/w;
			Y[k+1]	=  temp/w;
		}
		for ( i=-1; i >= nint; i-- ) {		/* differentiate */
			temp	=  Y[k];
			Y[k]	=  Y[k+1]*w;
			Y[k+1]	= -temp*w;
		}
	}
	for(i=1;i<=nint;i++)	Y[1] /= 1.0;		// PI/m/dt; // ??
	for(i=1;i<=nint;i++)	Y[2] /= (PI/dt);

#if FTDSP_DEBUG
	fp_debug = fopen("ftdsp_debug_file_3.dbg","w");
	(void) time(&now);
	fprintf(fp_debug,"%% ftdsp_debug_file_3 ... FILTERED(REALFT(TAPER(DETREND(y)))) ... %s\n",ctime(&now));
	fprintf(fp_debug,"%% k   freq.    real(y[k])   imag(y[k]) \n");
	for (k=1; k<=N; k+=2) {
		w = PI*(k-1)/n/dt;
		fprintf(fp_debug,"%d %e %e %e\n", k, w/(2*PI), Y[k], Y[k+1] ) ;
	}
	fclose(fp_debug);
	fprintf(stderr," ... 3\n");
#endif

	if (binLo > 1)		Y[1] = 0.0;
	if (binHi < n+2)	Y[2] = 0.0;

 	// taper the window in the frequency domain 
	for ( k = binLo - nBandLo; k <= binLo + nBandLo; k += 2 ) {
		taper = 0.5 * ( 1.0 + sin((k-binLo)*0.5*PI/nBandLo) ); // 1 + sin
                Y[k]   *= taper;
                Y[k+1] *= taper;
	}
	for ( k = binHi - nBandHi; k <= binHi + nBandHi; k += 2 ) {
		taper = 0.5 * ( 1.0 - sin((k-binHi)*0.5*PI/nBandHi) ); // 1 - sin
                Y[k]   *= taper;
                Y[k+1] *= taper;
	}

	if ( binLo < binHi ) {			// band pass filter 
		for (k = 3; k <= binLo - nBandLo; k += 2 ) Y[k] = Y[k+1] = 0.0;
		for (k = N; k >= binHi + nBandHi; k -= 2 ) Y[k] = Y[k-1] = 0.0;
	} else {					// band cut filter
		for (k=binHi+nBandHi; k<=binLo-nBandLo; k+=2) Y[k]=Y[k+1] = 0.0;
	}

#if FTDSP_DEBUG
	fp_debug = fopen("ftdsp_debug_file_4.dbg","w");
	(void) time(&now);
	fprintf(fp_debug,"%% ftdsp_debug_file_4 ... WINDOWED(FILTERED(REALFT(TAPER(DETREND(y))))) ... %s\n",ctime(&now));
	fprintf(fp_debug,"%% k   freq.    real(Y[k])   imag(Y[k]) \n");
	for (k=1; k<=N; k+=2) {
		w = PI*(k-1)/n/dt;
		fprintf(fp_debug,"%d %e %e %e\n", k, w/(2*PI), Y[k], Y[k+1] ) ;
	}
	fclose(fp_debug);
	fprintf(stderr," ... 4\n");
	fflush(stderr);
#endif

	realft ( Y, N, -1 );				/* inverse real FFT  */
	for (k=1; k<=N; k++) Y[k] *= (2.0/N);

#if FTDSP_DEBUG       
	fp_debug = fopen("ftdsp_debug_file_5.dbg","w");
	(void) time(&now);
	fprintf(fp_debug,"%% ftdsp_debug_file_5 ... INVREALFT(WINDOWED(FILTERED(REALFT(TAPER(DETREND(y)))))) ... %s\n",ctime(&now));
	fprintf(fp_debug,"%% y[k] \n");
	for (k=1; k<=N; k++) fprintf(fp_debug,"%11.4e\n", Y[k] );
	fclose(fp_debug);
#endif

	for ( i = 1; i <= n; i++ ) y[i] = Y[i+i1-1];

	free_vector ( Y, 1, N );

	return;
}

#undef FTDSP_DEBUG

/* SELECT_VALUE
 * Returns the kth smallest value in the array arr[1..n]. The input array will
 * be rearranged to have this value in location arr[k], with all smaller
 * elements moved to arr[1..k-1] (in arbitrary order) and all larger elements
 * in arr[k+1..n] (also in arbitrary order).
 * Adapted from Numerical Recipes in C, 2nd edition, 1992
 */
float selectValue ( float *arr, unsigned long n, unsigned long k )
{
	unsigned long i,r,j,l,mid;
	float a;

	l = 1;
	r = n;
	for (;;) {
	    if (r <= l+1) {	// Active partition contains 1 or 2 elements.
		if (r == l+1 && arr[r] < arr[l]) { SWAP(arr[l],arr[r]) }
 				// Case of 2 elements.
		return arr[k];
	    } else {
// Choose median of left, center, and right elements as partitioning element a. 
// Also rearrange so that arr[l] <= arr[l+1], arr[r] >= arr[l+1].
		mid = (l+r) >> 1;
		SWAP(arr[mid],arr[l+1])
		if (arr[l] > arr[r])		{ SWAP(arr[l],arr[r])   }
		if (arr[l+1] > arr[r])		{ SWAP(arr[l+1],arr[r]) }
		if (arr[l] > arr[l+1]) 		{ SWAP(arr[l],arr[l+1]) }
		i=l+1;			// Initialize pointers for partitioning.
		j=r;
		a=arr[l+1];			// Partitioning element.
		for (;;) {			// Beginning of innermost loop.
		    do i++; while (arr[i] < a); // Scan up   find element > a.
		    do j--; while (arr[j] > a); // Scan down find element < a.
		    if (j < i) break; // Pointers crossed. Partitioning done.
		    SWAP(arr[i],arr[j])
		}				// End of innermost loop.
		arr[l+1]=arr[j];		// Insert partitioning element.
		arr[j]=a;
		if (j >= k) r=j-1;	// Keep active the partition that 
		if (j <= k) l=i;	// contains the kth element.
	    }
	}
}
#undef SWAP
#undef DEBUG

