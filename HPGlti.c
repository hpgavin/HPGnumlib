/*******************************************************************************
HPGlti.c
routines for linear time invariant system simulation
*******************************************************************************/

/*  Copyright (C) 2023 Henri P. Gavin
 
    HPGlti is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. 
    
    HPGlti is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with HPGlti .  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>               // standard input/output library routines
#include <string.h>              // standard string handling library
#include <math.h>                // standard mathematics library

// local libraries ...
#include "HPGlti.h"              // header for linear time invariant system 
#include "NRutil.h"              // memory allocation routines


/* DT_SIM - evaluate x at time (k+1) 
  x(k+1) = A * x(k) + B * y(k)   ... for the Discrete Time Dynamics
-----------------------------------------------------------------------------*/
void DT_sim ( float **A, float **B, int n,  int m, float *x1, float *xc, float *y )
{
  int  i, j;

  for (i=1; i<=n; i++) {
    x1[i] = 0.0;
    for (j=1; j<=n; j++)  x1[i] += A[i][j] * xc[j];
    for (j=1; j<=m; j++)  x1[i] += B[i][j] * y[j];
  }
  return;
}



/* STATE_DERIV  -  calculate the state derivitive of a LTI system    17jul01
    dx/dt = A x + B y
-----------------------------------------------------------------------------*/
void state_deriv ( float **A, float **B, float *x, float *u, float *dxdt, int n, int m )
{
  int i,j;

  for (i=1; i<=n; i++) dxdt[i] = 0.0;
  for (i=1; i<=n; i++) for (j=1; j<=n; j++)  dxdt[i] += A[i][j]*x[j];
  for (i=1; i<=n; i++) for (j=1; j<=m; j++)  dxdt[i] += B[i][j]*u[j];

  return;
}


/* RK4  -  4th Order Runge-Kutta, ``Numerical Recipes In C,'' 2001-07-17 
---------------------------------------------------------------------------*/
void rk4 ( float *y, float *dydx, int n, float h, float vi, float *yout, 
  void (*derivs)() )
{
  int  i;
  float  hh, h6, *dym, *dyt, *yt;

  dym = vector(1,n);
  dyt = vector(1,n);
  yt  = vector(1,n);

  hh = h/2.0;
  h6 = h/6.0;

  for(i=1;i<=n;i++)
    yt[i] = y[i]+hh*dydx[i];

  (*derivs)(vi,yt,dyt);
  for(i=1;i<=n;i++)
    yt[i] = y[i]+hh*dyt[i];

  (*derivs)(vi,yt,dym);
  for(i=1;i<=n;i++) {
    yt[i] = y[i]+h*dym[i];
    dym[i] += dyt[i];
  }

  (*derivs)(vi,yt,dyt);
  for(i=1;i<=n;i++)
    yout[i] = y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);

  free_vector(yt,1,n);
  free_vector(dyt,1,n);
  free_vector(dym,1,n);
  return;
}


/* LSYM  -  simulate discrete-time linear time invariant system response  
---------------------------------------------------------------------------*/
void lsym( float **A, float **B, float **C, float **D, float *x, float **u, float **y, int n, int m, int r, int K )
{
  int i,j,k;
  float *x1;

  x1 = vector(1,n);

  for (k=1; k<=K; k++) for (i=1; i<=r; i++)  y[i][k] = 0.0;

  for (k=1; k<=K; k++) {
    for (i=1; i<=n; i++) x1[i] = 0.0;
    for (i=1; i<=n; i++) for (j=1; j<=n; j++)  x1[i]    += A[i][j] * x[j];
    for (i=1; i<=n; i++) for (j=1; j<=r; j++)  x1[i]    += B[i][j] * u[j][k];
    for (i=1; i<=m; i++) for (j=1; j<=n; j++)   y[i][k] += C[i][j] * x[j];
    for (i=1; i<=m; i++) for (j=1; j<=r; j++)   y[i][k] += D[i][j] * u[j][k];
    for (i=1; i<=n; i++)   x[i] = x1[i];
  }

  free_vector(x1, 1, n);

}


/* C2D - convert Continuous time dynamics to Discrete time dynamics

 Ak = A;
 k_factorial = 1;

 max_iter = 100;
 tolerance = 1e-6;

 for k = 1:max_iter

    expmA = expmA + Ak/k_factorial;

    if all( abs(Ak(find(expmA~=0))) ./ abs(expmA(find(expmA~=0)))/k_factorial < tolerance )
	break;
    end
    
    Ak = Ak*A;
    k_factorial = k_factorial*(k+1);
 end

								7 June 2007
---------------------------------------------------------------------------*/
void c2d( float	**A, float **B, float sr, int n, int m )
{
	float	k_factorial = 1.0,
		tlrnc = 0.000001,		/* convergence tolerance */
		**Mdt,
		**expmM,
		**Mk,
		**MkM,
		**matrix();

	int	i,j,q,r,k,
		converged,
		max_iter = 100;

	/* allocate memory for working matrices	*/
	Mdt   = matrix(1,n+m,1,n+m);
	expmM = matrix(1,n+m,1,n+m);
	Mk    = matrix(1,n+m,1,n+m);
	MkM   = matrix(1,n+m,1,n+m);

	for (i=1;i<=n;i++) {	 /* Mdt = [ Ac Bc ; 0 0 ];	*/
		for (j=1;j<=n;j++)	Mdt[i][j]    = A[i][j] / sr;
		for (j=1;j<=m;j++)	Mdt[i][n+j] = B[i][j] / sr;
	}
	for (i=n+1;i<=n+m;i++) for (j=1;j<=n+m;j++) Mdt[i][j]   = 0.0;

	/* compute the matrix exponential of a real square matrix M */
	for (i=1;i<=n+m;i++) {
		for (j=1;j<=n+m;j++)
			expmM[i][j]   = 0.0;
		expmM[i][i] = 1.0;
	}
	for (i=1;i<=n+m;i++)
		for (j=1;j<=n+m;j++)
			Mk[i][j] = Mdt[i][j];

	for ( k=1; k <= max_iter; k++ ) {

		for (i=1;i<=n+m;i++) 	 
			for (j=1;j<=n+m;j++)
				expmM[i][j] +=  Mk[i][j] / k_factorial;

		converged = 1;
		/* check for convergence on every term of expmM */
		for (i=1;i<=n+m;i++) 	 
		    for (j=1;j<=n+m;j++)
			if ( expmM[i][j] != 0.0 && converged ) 
			    if ( fabs(Mk[i][j]/expmM[i][j]/k_factorial) > tlrnc)
					converged = 0;
					/* term [i][j] has not converged */
		if (converged) break;


		for (q=1;q<=n+m;q++) {		/* Mk = Mk * Mdt; */
			for (r=1;r<=n+m;r++) {
				MkM[q][r] = 0.0;
				for (i=1;i<=n+m;i++)
					MkM[q][r] += Mk[q][i] * Mdt[i][r];
			}
		}
		for (i=1;i<=n+m;i++) for (j=1;j<=n+m;j++) Mk[i][j] = MkM[i][j];
	
		k_factorial = k_factorial*(k+1);
	}

	for (i=1;i<=n;i++) {	 /* expmM = [ Ad Bd ; 0 0 ];	*/
		for (j=1;j<=n;j++)	A[i][j]  = expmM[i][j];
		for (j=1;j<=m;j++)	B[i][j]  = expmM[i][n+j];
	}

	free_matrix(Mdt,1,n+m,1,n+m);
	free_matrix(expmM,1,n+m,1,n+m);
	free_matrix(Mk,1,n+m,1,n+m);
	free_matrix(MkM,1,n+m,1,n+m);

	return;
}


