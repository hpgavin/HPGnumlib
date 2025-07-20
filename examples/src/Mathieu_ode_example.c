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
#include <math.h>

#include "../../NRutil.h"
#include "../../HPGode.h"
#include "../../HPGsignal.h"
#include "../../HPGmatrix.h"

void Mathieu_sys(float *dxdt, float t, float *x, float *ui, float *y, float *c)
// Mathiew_sys
// compute the state derivitive of a Mathieu nonlinear dynamical system 
// a pendulum with vertical forcing at the pivot is an example of a Mathieu system 
{
  float   cg, cm, cl, ck, cc,        // constants
          ddu,                       // forcing
          q, dq, du, u,              // states 
          x1, x2, dx1, dx2,          // coordinates
          T, V;                      // energy

  // constants
  cg = c[1];                         // gravitational acceleration  m/s^2
  cm = c[2];                         // pendulum mass              kg
  cl = c[3];                         // pendulum length             m
  ck = c[4];                         // pendulum spring            N.m/rad
  cc = c[5];                         // pendulum damping           N.m/(rad/s)

  // states
  q  = x[1];                         // swing angle, rad
  dq = x[2];                         // swing angle velocity, rad/s
  du = x[3];                         // pin vertical velocity
  u  = x[4];                         // pin vertical position

  // forcing
  ddu = ui[1];                       // vertical pivot acceleration 

  // analysis
  x1  =  cl * sin(q);                // horizontal position
  x2  =  cl * cos(q) - cl + u;       //   vertical position

  dx1 =  cl * cos(q)*dq;             // horizontal velocity
  dx2 = -cl * sin(q)*dq + du;        //   vertical velocity

  T =  0.5*cm*(dx1*dx1 + dx2*dx2);   //    kinetic energy
  V =  0.5*ck*q*q - cm*cg*x2;        //  potential energy 

  // o.d.e's
  dxdt[1] =  dq; 
  dxdt[2] = -cm * cl * (cg + ddu) * sin(q) - ck*q - cc*dq; 
  dxdt[3] =  ddu; 
  dxdt[4] =   du;

  // other responses of interest
  y[1] = T; 
  y[2] = V;

  return;
}


int main (argc,argv)
int     argc;
char    *argv[];
{
  int    n=4,  m=1, l=2;             // number of states, inputs, outputs
  int    i,j;

  // constants
  float   cg = 9.81,                 // gravitational acceleration  m/s^2
          cm = 1.0,                  // pendulum mass               kg 
          cl = 1.0,                  // pendulum length             m
          ck = 0.200,                // pendulum stiffness        N.m/rad
          cc = 0.001;                // pendulum damping          N.m/(rad/s)

  // time 
  float   T  = 200,                  // total time duration         s
          dt =   0.01;               // time step                   s
  int     N  = floor(T/dt);          // number of time points

  // forcing
  float   w = 2*sqrt(cg/cl) * 1.00;  // forcing frequency  rad/s
  float  *u, *du, *ddu, **ui;        // forcing: pin postn, veloc, accel

  float   *c;                        // constants
  float   *time;                     // time 
  float   *x0;                       // initial states
  float  **x_drv;                    // state solution derivitives
  float  **x_sol;                    // state solution
  float  **y;                        // output

  FILE    *fp;                       // output file pointer
 
  // memory allocation 
  c     = vector(1,5);               // constants
  time  = vector(1,N);               // time
     u  = vector(1,N);               // forcing pin displacement 
    du  = vector(1,N);               // forcing pin velocity
   ddu  = vector(1,N);               // forcing pin acceleration 
     ui = matrix(1,m,1,N);           // forcing pin acceleration 
  x0    = vector(1,n);               // initial state vector
  x_drv = matrix(1,n,1,N);           // state derivitive solution
  x_sol = matrix(1,n,1,N);           // state solution 
  y     = matrix(1,l,1,N);           // output

  // constants 
  c[1] = cg; 
  c[2] = cm; 
  c[3] = cl;
  c[4] = ck;
  c[5] = cc;

  // time, s 
  for(j=1;j<=N;j++) time[j] = (j-1)*dt; 

  // forcing 
  for(j=1;j<=N;j++) 
    u[j] = 0.2 * cl * sin(w*time[j]) * pow( sin(PI*time[j]/T) , 2.0 );
   du = cDiff(u,N,dt);
  ddu = cDiff(du,N,dt);
  for(j=1;j<=N;j++) ui[1][j] = ddu[j];  // input to the ode solver

  // initial state
  x0[1] = 0.01*2*PI;   //  q
  x0[2] = 0.0;         // dq
  x0[3] = 0.0;         // du
  x0[4] = 0.0;         //  u


  // solve the system of o.d.e's 
  ode4u( Mathieu_sys, time, x0, ui, x_drv, x_sol, y, c, N, n,m,l );


  // write results to a file
  j = 1;
  fp = fopen("Mathieu_sys_data.txt","w"); 
  fprintf(fp,"%% (1) time\t(2) u\t\t(3) du\t\t(4) ddu\t\t");
  for(i=1;i<=n;i++) fprintf(fp,"(%d) x_drv_%d\t",    4+i , i ); 
  for(i=1;i<=n;i++) fprintf(fp,"(%d) x_sol_%d\t",  4+n+i , i ); 
  for(i=1;i<=l;i++) fprintf(fp,"(%d) y_%d\t",    4+2*n+i , i ); 
  fprintf(fp,"\n");
  for(j=1;j<=N;j++) {
    fprintf(fp,"%12.5e\t", time[j] ); 
    fprintf(fp,"%12.5e\t%12.5e\t%12.5e\t", u[j], du[j], ddu[j] ); 
    for(i=1;i<=n;i++) fprintf(fp,"%12.5e\t", x_drv[i][j] ); 
    for(i=1;i<=n;i++) fprintf(fp,"%12.5e\t", x_sol[i][j] ); 
    for(i=1;i<=l;i++) fprintf(fp,"%12.5e\t", y[i][j] ); 
    fprintf(fp,"\n");
  }

  fclose(fp);  // close the output data file

  // de-allocate memory
  free_vector(  u,1,N);
  free_vector( du,1,N);
  free_vector(ddu,1,N);
  free_vector(x0,1,n); 
  free_matrix(x_sol,1,n,1,N);
  free_matrix(x_drv,1,n,1,N);
  free_matrix(y,1,l,1,N); 

}
