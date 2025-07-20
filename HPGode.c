/*   HPGode.c
 *
 *    ordinary differential equation solvers
 *
 *   Author:  Henri P. Gavin (hpgavin), h p gavin ~at~ duke ~dot~ e d v
 *  Company:  Duke Univ.
 *
 * ==========================================================================
 */

/*  Copyright (C) 2001 Henri P. Gavin
 
    HPGode is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. 
    
    HPGode is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with HPGode.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h> 

#include "HPGode.h"
#include "NRutil.h"

/*
  ode4u - solve a system of nonhomogeneous ordinary differential equations 
  using the 4th order Runge-Kutta method.  
 
   Input Variable      Description                         
   --------------      ----------- 
   void (*ode_fctn)( float *dxdt, float t, float *x, float *u, float *y, float *c )
                     : a function which provides the state derivative given 
                       the state x, the time t, the input u, and other constants
       time          : a vector of instants in time at which the solution 
                       is computed 
        x0           : the initial value of the states 
        c            : vector of other constants 
        u            : m-by-p matrix of forcing for the ode, where m
                       m = number of inputs  and  n = length(time);
        N            : number of time values
        n            : number of states
        m            : number of inputs 
        l            : number of outputs 
 
  Output Variable        Description
  --------------         -----------
       time          : is returned, un-changed
       x_sol         : the solution to the differential equation at 
                       times given in the vector time.
 
  Henri Gavin, Civil and Environmental Engineering, Duke Univsersity, Jan. 2001
  WH Press, et al, Numerical Recipes in C, 1992, Section 16.1
*/

void ode4u( 
  void (*ode_fctn)(float *, float,  float *, float *, float *, float *), // the ode function 
  float *time,        // a vector of time values
  float *x0,          // the initial state vector
  float *u[],         // the vector of inputs at each instant in time
  float *x_drv[],     // the vector of the state derivitive at each instant in time
  float *x_sol[],     // the vector of the solution at each instant in time
  float *y[],         // system outputs 
  float *c,           // a vector of other constant values
  int N,              // the number of time values
  int n,              // the dimension of the state vector
  int m,              // the number of inputs 
  int l               // the number of outputs
) 
{
  float t, dt, dt2;
  float *xi, *ui, *yi, *dxdt1, *dxdt2, *dxdt3, *dxdt4;
  int     i, p;
//float *vector();
//void    free_vector();

  xi = vector(1,n);
  yi = vector(1,l);
  ui = vector(1,m);

  dxdt1 = vector(1,n);
  dxdt2 = vector(1,n);
  dxdt3 = vector(1,n);
  dxdt4 = vector(1,n);

  // initialization
  for(i=1;i<=m;i++) ui[i] = u[i][1];             // initial input
  (*ode_fctn)( dxdt1, time[1], x0, ui, yi, c );  // initial ode 
  for(i=1;i<=n;i++) x_sol[i][1] = x0[i];         // initial state
  for(i=1;i<=n;i++) x_drv[i][1] = dxdt1[i];      // initial dxdt
  for(i=1;i<=l;i++) y[i][1] = yi[i];             // initial output

  for ( p=1; p<N; p++ ) {

    t   = time[p];
    dt  = time[p+1]-t;
    dt2 = dt/2.0;

    for(i=1;i<=m;i++) ui[i] = u[i][p];
    for(i=1;i<=n;i++) xi[i] = x0[i] + dxdt1[i]*dt2;
    for(i=1;i<=m;i++) ui[i] = (u[i][p] + u[i][p+1])/2.0;

    (*ode_fctn)( dxdt2, t+dt2, xi, ui, yi, c );

    for(i=1;i<=n;i++) xi[i] = x0[i] + dxdt2[i]*dt2;

    (*ode_fctn)( dxdt3, t+dt2, xi, ui, yi, c );

    for(i=1;i<=n;i++) xi[i] = x0[i] + dxdt3[i]*dt;
    for(i=1;i<=m;i++) ui[i] = u[i][p+1];

    (*ode_fctn)( dxdt4, time[p+1] , xi, ui, yi, c );

    for(i=1;i<=n;i++) // update x0
      x0[i] += (dxdt1[i] + 2.0*(dxdt2[i] + dxdt3[i]) + dxdt4[i])*dt/6.0;  // RK4

    (*ode_fctn)( dxdt1, time[p+1], x0, ui, yi, c );

    for(i=1;i<=n;i++) x_sol[i][p+1] = x0[i];    // store state
    for(i=1;i<=n;i++) x_drv[i][p+1] = dxdt1[i]; // store state derivitive
    for(i=1;i<=l;i++) y[i][p+1]     = yi[i];    // store output

    for(i=1;i<=n;i++) if(x0[i] > 1e30 || x0[i] < -1e30) break;

  }

  free_vector(dxdt1,1,n);
  free_vector(dxdt2,1,n);
  free_vector(dxdt3,1,n);
  free_vector(dxdt4,1,n);

  free_vector(xi,1,n);
  free_vector(ui,1,m);
  free_vector(yi,1,l);

  return;

}
