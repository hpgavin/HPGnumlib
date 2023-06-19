/* HPGrandom.c - library for random number generator functions */

/* Copyright (C) 2012 Henri P. Gavin
 
    HPGrandom is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. 
    
    HPGrandom is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with HPGrandom .  If not, see <http://www.gnu.org/licenses/>.
*/


#include "HPGrandom.h"
#include "HPGcholesky.h"
#include "NRutil.h"

#define DEBUG 0


/* ---------------------------------------------------------------------------
 * INITIALIZE_RANDOM
 * initialize the sequence of pseudo-random numbers with usec time of day
 * ---------------------------------------------------------------------------
 */
void initialize_random ( unsigned int seed )
{
	struct  timeval tv;
	struct  timezone tz;

	if ( seed <= 0 ) {
		gettimeofday(&tv, &tz);
		seed = (unsigned int)tv.tv_usec*100 + (unsigned int)time(NULL);
	}
	if (DEBUG) printf(" random number generator seed = %d \n", seed );
	srandom(seed);
}


/* ----------------------------------------------------------------------------
 * INVNORM
 * inverse of the standard normal cumulative distribution function
 *  z = sqrt(2) * erfinv ( 2*u - 1 );
 * ----------------------------------------------------------------------------
 */
double invnorm ( double u )
{
	double	z;

	z = Sqrt2 * erfinv ( 2*u - 1 );

	return z;
}


/* ---------------------------------------------------------------------------
 * ERFINV
 * Newton iterations for the inverse error function 
 * ---------------------------------------------------------------------------
 */
#define MAXIT 100
#define TOL 1e-12

double erfinv ( double y )
{
	int	iterations = 0;	
	double	x_old = 1.0, x_new = 0.0;
	
	if ( y <= -1.0 ) 	return -HUGE_VAL ;	/* minus "infinity" */

	if ( y >=  1.0 )	return  HUGE_VAL ;	/* plus  "infinity" */

	x_old = 1.0;
	if ( y >= 0.0 ) x_new =  sqrt ( -log( 1 - fabs(y))); 
	if ( y <= 0.0 ) x_new = -sqrt ( -log( 1 - fabs(y))); 

	while ( fabs(erf(x_new) - y )  >  TOL*fabs(y) ) {
		x_old = x_new;
		x_new = x_old - (erf(x_old) - y) * exp(x_old*x_old) * SqrtPi2;
		if (++iterations > MAXIT) {
			fprintf(stderr,"erfinv: iteration limit exceeded\n");
			break;
		}
	}
	if (DEBUG) fprintf(stderr," iteration = %4d \n", iterations );

	return x_new;
}
#undef MAXIT
#undef TOL


/* --------------------------------------------------------------------------
 * LOGNORMAL
 * generate n lognormal variables X with mean M, coefficient of variation COV,
 * and correlation Cz such that X is within +/- S standard deviations of M
 * --------------------------------------------------------------------------
 */
void lognormal (
 float  *X,	// the computed log normal random variables
 float  *M,	// the means of the log normal random variables 
 float  *COV,	// the coefficient of variation of the log normal variables
 float	**Cz,	// the lower Cholesky factor of the correl. mtrx 
 float  *Czd,	// the diagonal of the lower Cholesky factor of the correl mtrx
 float	S,	// the number of allowabe standard deviations from mean
 int    n  )
{
	int	i;
	float	*Zu, *Zc, *Vlnx; 
	int	ok = 0,
		reduce = 0;	// Cz is already reduced to Cholesky factor

	Vlnx = vector(1,n);
	Zu   = vector(1,n);
	Zc   = vector(1,n);

	for (i=1;i<=n;i++)      Vlnx[i] = log ( COV[i]*COV[i] + 1.0 );

	while (ok==0) {

	    for (i=1;i<=n;i++) 
		Zu[i] = Sqrt2 * erfinv( 2.0*((float)random()/RAND_MAX) - 1.0 );


	    cholesky_correlate ( Cz, Czd, n, Zu, Zc, reduce );


	    for (i=1; i<=n; i++)
	      X[i] = M[i] * exp ( -0.5*Vlnx[i] + Zc[i] * sqrt(Vlnx[i]) );

            if (DEBUG) 
	      for (i=1;i<=n;i++) {
		printf("   M[%d]=%7.3f  COV[%d]=%7.3f", i, M[i], i, COV[i] );
		printf("   Zu[%d]=%7.3f  Zc[%d]=%7.3f  X[%d]=%7.3f \n", 
						i, Zu[i], i, Zc[i], i, X[i] );
	      }

	    ok = 1;
	    for (i=1; i<=n; i++) if ( Zc[i] >  S )	ok = 0;
	    for (i=1; i<=n; i++) if ( Zc[i] < -S )	ok = 0;

	}

	free_vector(Vlnx,1,n);
	free_vector(Zu,1,n);
	free_vector(Zc,1,n);
}


/* --------------------------------------------------------------------------
 * INVERSE_NORMAL -
 * returns the standard normal random deviate of prob p 
 * --------------------------------------------------------------------------
 */

#define C0 2.515517
#define C1 0.802853
#define C2 0.010328
#define D1 1.432788
#define D2 0.189269
#define D3 0.001308

/*
#define C0   2.7062028586760798
#define C1   1.1372722114755480
#define C2  -0.1824423322365641
#define D1   1.8094831320366702
#define D2   0.1184943586698157
#define D3  -0.0452019252916463
*/


float inverse_normal( float p )
{
	float	t,z;		/* p = Phi(z)   Phi: Normal CDF */
	int	sign = 1;

	if (p > 0.5) {
		p=1-p;
		sign = -1;
	}
	t = sqrt(-2*log(p));
	z = -t + (C0 + C1*t + C2*t*t) / (1 + D1*t + D2*t*t + D3*t*t*t);
	z *= sign;
	return(z);
}

#undef C0
#undef C1
#undef C2
#undef D1
#undef D2
#undef D3




/* ---------------------------------------------------------------------------
 * LTQNORM
 * Lower tail quantile for standard normal distribution function.
 *
 * 	http://home.online.no/~pjacklam/notes/invnorm/
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http://www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 * ---------------------------------------------------------------------------
 */


/* Coefficients in rational approximations. */
static const double a[] =
{
	-3.969683028665376e+01,
	 2.209460984245205e+02,
	-2.759285104469687e+02,
	 1.383577518672690e+02,
	-3.066479806614716e+01,
	 2.506628277459239e+00
};

static const double b[] =
{
	-5.447609879822406e+01,
	 1.615858368580409e+02,
	-1.556989798598866e+02,
	 6.680131188771972e+01,
	-1.328068155288572e+01
};

static const double c[] =
{
	-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
	 4.374664141464968e+00,
	 2.938163982698783e+00
};

static const double d[] =
{
	7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00
};

#define LOW 0.02425
#define HIGH 0.97575

double ltqnorm ( double p )
{
	double q, r;

	if ( p < 0 || p > 1 )	return 0.0;

	if (p == 0)	return -HUGE_VAL /* minus "infinity" */;

	if (p == 1)	return HUGE_VAL /* "infinity" */;

	if (p < LOW)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2.0*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else if (p > HIGH)
	{
		/* Rational approximation for upper region */
		q  = sqrt(-2.0*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else
	{
		/* Rational approximation for central region */
    		q = p - 0.5;
    		r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}


#undef DEBUG
