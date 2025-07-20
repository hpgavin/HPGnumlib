/* =====================================================================================
 *   Mathieu_ode_example.c 
 *  
 *   example of solving a system of ode's using ode4u.c
 *
 *   2025-07-20
 *
 * gcc -O -o Mathieu_ode_example  Mathiew_ode_example.c  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "NRutil.h"
#include "HPGode.h"
#include "HPGsignal.h" 

#define PI 3.14159265358979323846264338327950288419716939937510


void Mathieu_sys( float *dzdt, float t, float *x, float *ddu, float *y, float *c)
// Mathiew_sys
// compute the state derivitive of a Mathieu nonlinear dynamical system 
// a pendulum with vertical forcing at the pivot is an example of a Mathieu equation 
{
  double  cg, cm, cl, cc,            // constants
          q, dq, du, u,              // states 
          x1, x2, dx1, dx1,          // coordinates
          T, V;                      // energy

  // constants
  cg = c[1];                         // gravitational acceleration  m/s^2
  cm = c[2];                         // pendulum mass              kg
  cl = c[3];                         // pendulum length             m
  cc = c[4];                         // pendulum damping           N.m/(rad/s)

  // states
  q  = x[1];                         // swing angle, rad
  dq = x[2];                         // swing angle velocity, rad/s
  du = x[3];                         // pin vertical velocity
  u  = x[4];                         // pin vertical position

  // forcing
  ddu = u[1];                        // pivot acceleration forcing

  // analysis
  x1  =  cl * sin(q);                // horizontal position
  x2  =  cl * cos(q) - c.l + u;      //   vertical position

  dx1 =  cl * cos(q)*dq;             // horizontal velocity
  dx2 = -cl * sin(q)*dq + du;        //   vertical velocity

  T =  0.5*cm*(dx1*dx1 + dx2*dx2);   //  kinetic  energy
  V = -cm*cg*x2;                     // potential  energy 

  // o.d.e's
  dzdt[1] =  dq; 
  dzdt[2] = -cm * cl * (cg + ddu) * sin(q) - cc*dq; 
  dzdt[3] =  ddu; 
  dzdt[4] =   du;

  //  other responses of interest
  y[1] = T; 
  y[2] = V;

  return;
}


int main (argc,argv)
int     argc;
char    *argv[];
{
  FILE    *fp;               // output file pointer
 
  int    n=4,  m=1, l=2;     // number of states, inputs, outputs
  int    i,j,k;

  // constants
  double  cg = 9.81,         // gravitational acceleration  m/s^2
          cm = 1.0,          // pendulum mass               kg 
          cl = 1.0,          // pendulum length             m
          cc = 0.001;        // pendulum damping            N.m/(rad/s)

  // time 
  double  T  = 220;          // total time duration
          dt =   0.05;       // time step
  int     N  = floor(T/dt);  // number of time points

  // forcing
  double  w = 2*sqrt(cg/cl) * 1.00;    // forcing frequency  rad/s
  double **u, **du, **ddu;       // forcing

  double *t;                     // time 
  double *c;                     // constants
  double *z0;                    // initial states
  double **y;                    // output

  // memory allocation 
  c     = dvector(1,4);          // constants
    u   = dmatrix(1,m,1,N);      // forcing pin displacement 
   du   = dmatrix(1,m,1,N);      // forcing pin velocity
  ddu   = dmatrix(1,m,1,N);      // forcing pin acceleration 
  z0    = dvector(1,n);          // initial state vector
  z_drv = dmatrix(1,n,1,N);      // state derivitive solution
  z_sol = dmatrix(1,n,1,N);      // state solution 
  y     = dmatrix(1,l,1,N);      // output

  // constants 
  c[1] = cg; 
  c[2] = cm; 
  c[3] = cl;
  c[4] = cc;

  // time, s 
  for(j=1;j<=N;j++) t[j] = (j-1)*dt; 

  // forcing 
  for(j=1;j<=N;j++) 
    u[1][j] = 0.200 * cl * sin(w*t[j]) * pow( sin(PI*t[j]/T) , 2 );
   du = cDiff(u,N,dt);
  ddu = cDiff(du,N,dt);

  // initial state
  z0[1] = 0.01*2*PI;  
  z0[2] = 0.0;
  z0[3] = 0.0;
  z0[4] = 0.0;


  // solve the system of o.d.e's 
  ode4u( Mathieu_sys, time, z0, u, z_drv, z_sol, y, c, N, n,m,l );


  // write results to a file
  fp = fopen("Mathieu_sys_data.txt","w"); 
  for(j=1;j<=N;j++) {
    fprintf(fp,"%12.5e\t", t[j] ); 
    fprintf(fp,"%12.5e\t%12.5e\t%12.5e\t", u[1][j], du[1][j], ddu[1][j] ); 
    for(i=1;i<=n;i++) fprintf(fp,"%12.5e\t", x_drv[i][j] ); 
    for(i=1;i<=n;i++) fprintf(fp,"%12.5e\t", x_sol[i][j] ); 
    for(i=1;i<=l;i++) fprintf(fp,"%12.5e\t", y[i][j] ); 
    fprintf(fp,"\n");
  }

  fclose(fp);  // close the output data file

  // de-allocate memory
  free_dvector(z0,1,n); 
  free_dmatrix(x_sol,1,n,1,N);
  free_dmatrix(x_drv,1,n,1,N);
  free_dmatrix(u,1,m,1,N);
  free_dmatrix(du,1,m,1,N);
  free_dmatrix(ddu,1,m,1,N);
  free_dmatrix(y,1,l,1,N); 

}
