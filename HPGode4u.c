/*
  Solve a system of nonhomogeneous ordinary differential equations using the 
  4th order Runge-Kutta method.  
 
   Input Variable      Description                         
   --------------      ----------- 
       dxdt          : a function of the form dxdt(t,x,u,params,dxdt)
                       which provides the state derivative given 
                       the state, x, and the time, t  
       time          : a vector of points in time at which the solution 
                       is computed 
        x0           : the initial value of the states 
      params         : vector of other parameters for dxdt
        u            : m-by-p matrix of forcing for the ode, where m
                       m = number of inputs  and  n = length(time);
      points         : number of data points
        n            : size of state vector
	m	     : size of the inpute vector
 
  Output Variable        Description
  --------------         -----------
       time          : is returned, un-changed
       x_sol         : the solution to the differential equation at 
                       times given in the vector time.
 
  Henri Gavin, Civil and Environmental Engineering, Duke Univsersity, Jan. 2001
  WH Press, et al, Numerical Recipes in C, 1992, Section 16.1
*/

void HPGode4u( dxdt, time, x0, u, params, points, n, m, x_sol )
void	(*dxdt)();
float	*time, *x0, *params, **u, **x_sol;
int	points, n, m;
{
	float	t, dt, dt2,
		*k1, *k2, *k3, *k4, 
		*xx, *uu,
		*vector();
	int	i,p;	
	void	free_vector();


	k1 = vector(1,n);
	k2 = vector(1,n);
	k3 = vector(1,n);
	k4 = vector(1,n);
	xx = vector(1,n);
	uu = vector(1,m);

	for (i=1; i<=n; i++) {
		x_sol[1][i] = x0[i];	/* initial conditions */
		k1[i] = k2[i] = k3[i] = k4[i] = 0.0;
	}

	for ( p=2; p<=points; p++ ) {

		t      = time[p]; 
		dt     = t - time[p-1];	/* the time step for this interval */
		dt2    = dt/2;		/* half of the time step	*/

		for (i=1;i<=n;i++)	xx[i] = x0[i];
		for (i=1;i<=m;i++)	uu[i] = u[p-1][i];
		(*dxdt)( t - dt,  xx, uu, params, k1, n,m );

		for (i=1;i<=n;i++)	xx[i] = x0[i]+k1[i]*dt2;
		for (i=1;i<=m;i++)	uu[i] = (u[p-1][i]+u[p][i])/2.0;
		(*dxdt)( t - dt2, xx, uu, params, k2, n,m );

		for (i=1;i<=n;i++)	xx[i] = x0[i]+k2[i]*dt2;
		(*dxdt)( t - dt2, xx, uu, params, k3, n,m );

		for (i=1;i<=n;i++)	xx[i] = x0[i]+k3[i]*dt;
		for (i=1;i<=m;i++)	uu[i] = u[p][i];
		(*dxdt)( t,       xx, uu, params, k4, n,m );

		for (i=1; i<=n; i++) {
			x0[i] = x0[i] + ( k1[i] + 2*(k2[i] + k3[i]) + k4[i] )*dt/6; 
			x_sol[p][i] = x0[i];
		}
	}

	free_vector(k1,1,n);
	free_vector(k2,1,n);
	free_vector(k3,1,n);
	free_vector(k4,1,n);
	free_vector(xx,1,n);
	free_vector(uu,1,m);

	return;
}
