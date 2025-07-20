/*  HPGode.h -  header files for functions in HPGode.c
 *
 *  Henri P. Gavin (hpgavin), h p gavin ~at~ duke ~dot~ e d v
 *
 * ==========================================================================
 */

/*
 Copyright (C) 2001 Henri P. Gavin
 
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

/*
  ode4u - solve a system of nonhomogeneous ordinary differential equations 
  using the 4th order Runge-Kutta method.  
  
   Input Variable      Description                         
   --------------      ----------- 
   void (*ode_fctn)( float  *dxdt, float  t, float  *x, float  *u, float  *params)
      *ode_fctn      : a function which provides the state derivative given 
                       the state x, the time t, the input u, and other  params
       time          : a vector of points in time at which the solution 
                       is computed 
        x0           : the initial value of the states 
      params         : vector of other parameters for dxdt
        u            : m-by-p matrix of forcing for the ode, where m
                       m = number of inputs  and  n = length(time);
      points         : number of data points
        n            : size of state vector
        m            : size of the inpute vector
 
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
  float  *time,      // a vector of time values
  float  *x0,        // the initial state vector
  float  *u[],       // the vector of inputs at each instant in time
  float  *x_sol[],   // the vector of the solution at each instant in time
  float  *x_drv[],   // the vector of the state derivitive at each instant in time
  float  *y[],       // outputs
  float  *c,         // a vector of constants
  int N,             // the number if time values
  int n,             // the dimension of the state vector
  int m,             // the number of inputs 
  int l 
);

