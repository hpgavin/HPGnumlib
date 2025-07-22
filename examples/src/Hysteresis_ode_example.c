/* =====================================================================================
 *   Hysteresis_ode_example.c 
 *  
 *   example of solving a system of ode's using ode4u.c
 *
 *   2025-07-20
 *
 * gcc -O -o Hysteresis_ode_example  Mathiew_ode_example.c  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <math.h>

#include "../../NRutil.h"
#include "../../HPGode.h"
#include "../../HPGutil.h"

void Hysteresis_sys(float *dxdt, float t, float *x, float *ui, float *y, float *c)
// Hysteresis_sys
// compute the state derivitive of a Hysteresis nonlinear dynamical system 
// a mass sliding on a moving horizontal surface
{
  float   cG, cU, cD,                // constants
          aS,                        // forcing - surface acceleration 
          r, dr, z,                  // states 
          fi;                        // sliding force

  // constants
  cG = c[1];                         // natural frequency         rad/s
  cU = c[2];                         // damping ratio                .
  cD = c[3];                         // friction coefficient         .

  // states
  r  = x[1];                         // position                     m
  dr = x[2];                         // velocity                     m/s
  z  = x[3];                         // friction force / cu          

  // forcing
  aS = ui[1];                        // ground acceleration 

  // analysis
  fi = cU*cG*z;                      // friction force

  // o.d.e's
  dxdt[1] =  dr; 
  dxdt[2] = -fi - aS;
  dxdt[3] =  ( 1.0 - pow(z,3.0)*sgn(dr) ) * dr / cD;  // sliding force derivitive

  // other responses of interest
  y[1] = fi; 

  return;
}


int main (argc,argv)
int     argc;
char    *argv[];
{
  int    n=3,  m=1, l=1;             // number of states, inputs, outputs
  int    i,j;

  // constants
  float   cG  = 9.81,                // gravitational acceleration  m/s^2
          cU = 0.10,                 // coefficient of friction 
          cD = 0.001;                // pre-slip displacement        m

  // time 
  float   T  = 2.0,                  // total time duration         s
          dt =   0.01;               // time step                   s
  int     N  = floor(T/dt);          // number of time points

  // forcing
  float  Tp = 0.5,                   // forcing pulse period         s
         Ap = 1.0;                   // forcing pulse amplitude     m/s^2
  float *aS, **ui; 

  float   *c;                        // constants
  float   *time;                     // time 
  float   *x0;                       // initial states
  float  **x_drv;                    // state solution derivitives
  float  **x_sol;                    // state solution
  float  **y;                        // output

  FILE    *fp;                       // output file pointer
 
  // memory allocation 
  c     = vector(1,3);               // constants
  time  = vector(1,N);               // time
   aS   = vector(1,N);               // forcing pin acceleration 
   ui   = matrix(1,m,1,N);           // forcing pin acceleration 
  x0    = vector(1,n);               // initial state vector
  x_drv = matrix(1,n,1,N);           // state derivitive solution
  x_sol = matrix(1,n,1,N);           // state solution 
  y     = matrix(1,l,1,N);           // output

  // constants 
  c[1] = cG; 
  c[2] = cU; 
  c[3] = cD;

  // time, s 
  for(j=1;j<=N;j++) time[j] = (j-1)*dt; 

  // forcing 
  for(j=1;j<=N;j++)              aS[j] = 0.0;
  for(j=1;j<=floor(Tp/dt);j++)   aS[j] = -Ap*sin(2.0*PI*time[j]/Tp); 

  for(j=1;j<=N;j++) ui[1][j] = aS[j];  // input to the ode solver

  // initial state
  x0[1] = 0.0;         //  r
  x0[2] = 0.0;         // dr
  x0[3] = 0.0;         //  z


  // solve the system of o.d.e's 
  ode4u( Hysteresis_sys, time, x0, ui, x_drv, x_sol, y, c, N, n,m,l );

  // write results to a file
  j = 1;
  fp = fopen("Hysteresis_sys_data.txt","w"); 
  fprintf(fp,"%% (1) time\t(2) aS\t\t");
  for(i=1;i<=n;i++) fprintf(fp,"(%d) x_drv_%d\t",    2+i , i ); 
  for(i=1;i<=n;i++) fprintf(fp,"(%d) x_sol_%d\t",  2+n+i , i ); 
  for(i=1;i<=l;i++) fprintf(fp,"(%d) y_%d\t",    2+2*n+i , i ); 
  fprintf(fp,"\n");
  for(j=1;j<=N;j++) {
    fprintf(fp,"%12.5e\t", time[j] ); 
    fprintf(fp,"%12.5e\t", aS[j] ); 
    for(i=1;i<=n;i++) fprintf(fp,"%12.5e\t", x_drv[i][j] ); 
    for(i=1;i<=n;i++) fprintf(fp,"%12.5e\t", x_sol[i][j] ); 
    for(i=1;i<=l;i++) fprintf(fp,"%12.5e\t", y[i][j] ); 
    fprintf(fp,"\n");
  }

  fclose(fp);  // close the output data file

  // de-allocate memory
  free_vector(x0,1,n); 
  free_matrix(x_sol,1,n,1,N);
  free_matrix(x_drv,1,n,1,N);
  free_matrix(y,1,l,1,N); 

}
