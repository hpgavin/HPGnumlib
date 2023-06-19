/* HPGlti.h
   routines for linear time invariant system simulation
*******************************************************************************/


/*---------------------------------------------------------------------------
C2D - convert Continuous time dynamics to Discrete time dynamics


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
void c2d( float **A, float **B, float sr, int n, int m );


/* LSYM  -  simulate discrete-time linear time invariant system response  
---------------------------------------------------------------------------*/
void lsym( float **A, float **B, float **C, float **D, float *x, float **u, float **y, int n, int m, int r, int K );


/* DT_SIM - evaluate x at time (k+1) 
  x(k+1) = A * x(k) + B * y(k)   ... for the Discrete Time Dynamics
-----------------------------------------------------------------------------*/
void DT_sim ( float **A, float **B, int _N,  int _L, float *xc1, float *xc, float *y );


/* STATE_DERIV  -  calculate the state derivitive of a LTI system    17jul01
    dx/dt = A x + B y
-----------------------------------------------------------------------------*/
void state_deriv ( float **A, float **B, float *x, float *u, float *dxdt, int n, int m );


/* RK4  -  4th Order Runge-Kutta, ``Numerical Recipes In C,'' 2001-07-17 
---------------------------------------------------------------------------*/
void rk4 ( float *y, float *dydx, int n, float h, float vi, float *yout, void (*derivs)() );


