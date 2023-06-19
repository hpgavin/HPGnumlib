/* HPGcholesky.c  
 *
 * Cholesky decomposition and related routines adapted from ...
 * Press, Numerical Recipes in C, 2nd ed., Cambridge Univ. Pr., 1988, Ch.2   
 */

#include "HPGcholesky.h"
#include "NRutil.h"

/* CHOLESKCY_DCMP
 *
 * Given a positive-definite symmetric matrix A[1..n][1..n], this
 * routine constructs its Cholesky decomposition, [A] = [L][L]'.
 *
 * If (reduce > 0) ...
 * On input, only the upper triangle of  [A]  need be given; it is not modified.
 * The Cholesky factor, [L], is returned in the lower triangle of [A],
 * except for its diagonal elements, which are returned in p[1..n].
 * 
 * If (solve > 0)
 * Solves the set of n linear equations [A]{x} = {b}, where [A] is a 
 * positive-definite symmetric matrix.  A[1..n][1..n] and p[1..n] are input
 * as the output of the routine choldc(). 
 * Only the lower triangle of [A] is accessed.  
 * b[1..n] is input as the right-hand side vector.  
 * The solution vector is returned in x[1..n].   
 * [A], n, and {p} are not modified and can be left in place for successive
 * calls with different right-hand sideds {b}. 
 * {b} is not modified unless you identify b and x in the calling sequence, 
 * which is allowed.
 * 
 *   [A] is a diagonally dominant matrix of dimension [1..n][1..n]. 
 *   {b} is a r.h.s. vector of dimension [1..n].
 *   {b} is updated using [LU] and then back-substitution is done to obtain {x}.
 *   
 *
 * if i>j: L[j][i]=L[i][j]=A[i][j]; L[i][i] = p[i]; if i<=j: A[i][j] = A[i][j]
 * ----------------------------------------------------------------------------
 */
void cholesky_dcmp (
 double **A,	// the system matrix, and it's Cholesky decomposition
 int n,		// the dimension of the matrix
 double *p,	// the vector of the diagonal of the Cholesky decomp'n
 double *b,	// the right hand side vector
 double *x,	// the solution vector
 int reduce,	// 1: do a forward reduction; 0: don't do the reduction
 int solve,	// 1: do a back substitution; 0: do no back substitution
 int *pd	// 1: positive definite; 0: not positive definite
){
	float	sum;
	int	i,j,k;

	if ( reduce ) {		// forward reduction of [A]
	  *pd = 1;
	  for (i=1; i<=n; i++) {
	    for (j=i; j<=n; j++) {
		for ( sum = A[i][j],k=i-1; k>=1; k-- )	sum -= A[i][k]*A[j][k];
		if ( i == j ) {
		    if ( sum <= 0.0 ) {
		      fprintf(stderr," choldc: matrix not positive-definite ");
		      fprintf(stderr," ( p[%d]*p[%d] = %11.4e )\n", i, i, sum );
		      *pd = 0;
		      return;
		    }
		    p[i] = sqrt(sum);
		} else	A[j][i] = sum / p[i];
	    }
	  }
	}			// the forward recution of [A] is now complete

	if ( solve ) { 		// back substitution to solve for {x}
	    for (i=1; i<=n; i++) {	// solve L y = b, storing y in x 
		for (sum=b[i],k=i-1; k>=1; k--)	sum -= A[i][k]*x[k];
		x[i]=sum/p[i];
	    }
	    for (i=n;i>=1;i--) {	// solve L' x = y
		for (sum=x[i],k=i+1; k<=n; k++)	sum -= A[k][i]*x[k];
		x[i]=sum/p[i];
	    }
	}

	return;
}


/* CHOLESKY_CORRELATE
 * Computes random variables X with a given correlation (or covariance)
 *  matrix C from a set of uncorrelated Guassean random variables Z
 *  C is the n-by-n correlation (or covariance) matrix
 *  Z is the n-by-1 matrix of uncorrleated Gaussean random variables
 *  X is the n-by-1 matrix of correlated Gaussean random variables
 *  The upper triangle of C, including the diagonal is returned unchanged
 *  The lower triangle of C is returned as the Cholesky factor of C
 * ----------------------------------------------------------------------------
 */
void cholesky_correlate (
 float **C,	// lower triangle of Cz is Cholesky factor of correl. matx
 float *Czd,	// diagonal of Cholesky factor of correl. matx
 int n,		// dimension of the correlation matrix
 float *Zu,	// uncorrelated standard normal random variables
 float *Zc,	// correlated standard normal random variables
 int reduce )	// 1: Cz is not yet factored, 0: Cz is already factored
{
	int	i,j, pd;
	double	*v;

	v = dvector(1,n);

	if ( reduce ) 
		cholesky_dcmp ( (double **) C, n, (double *) Czd, v, v, reduce, 0, &pd );

	for (i=1; i<=n; i++) {
		Zc[i] = Czd[i]*Zu[i];
		for (j=i-1; j>=1; j--)
			Zc[i] += C[i][j]*Zu[j];
	}

	free_dvector(v,1,n);
}


/* CHOLESKY_INV  - 
 * Inverts the Cholesky factor, [L], given  A  as the output of cholesky_dcmp.
 * The lower triangle of  [A]  is replaced by [L]^(-1).
 * ----------------------------------------------------------------------------
 */
void cholesky_inv ( double **A, int n, double *p )
{
	float	sum;
	int	i,j,k;

	for (i=1; i <= n; i++) {
		A[i][i] = 1.0/p[i]; 
		for (j = i+1; j <= n; j++) {
			sum = 0.0;
			for (k=1; k < j; k++)	sum -= A[j][k]*A[k][i];
			A[j][i] = sum/p[j];
		}
	}
	return;
}


/* CHOLESKY_SYM  -  
 * Standardizes the generalized eigen-problem [A]{x} = {d} [B] {x} into
 * standard form, C (L'x) = d (L'x).  The matrix [B] must be positive-definite.
 * Matrix [C] is output, which has the same eigen-values as the generalized
 * eigen-problem.  
 * The eigen-vectors of the generalized eigen-problem can be recovered
 * from [B] = [L][L]' and the eigen-vectors of [C].   [A], [B], and [C] are 
 * of dimension [1..n][1..n].  [A] must be symmetric, in this code.   
 * Only the upper triangle of [A] is accessed.
 * L[j][i]=L[i][j]=B[i][j] i>j	L[i][i] = p[i]		B[i][j] = B[i][j] i<=j
 * ----------------------------------------------------------------------------
 */ 
void cholesky_sym ( double **A, double **B, int n, double **C, double *p )
{
	double	**Y, *v, sum;
	int	i,j,k, pd;

	Y = dmatrix(1,n,1,n);
	v = dvector(1,n);

	cholesky_dcmp ( B, n, p, v, v, 1, 0, &pd  );	// B = L L'

	for (j=1; j<=n; j++) {				// Y L' = A	
		for(i=j; i<=n; i++) {
			sum = 0.0;
			if (j > 1) for (k=1; k<j; k++) sum += Y[i][k]*B[j][k];
			Y[i][j] = ( A[j][i] - sum ) / p[j];
		}
	}
	for (j=1; j<=n; j++) {				// L C = Y
		for (i=j; i<=n; i++) {
			sum = 0.0;
			if (i > 1) for (k=1; k<i; k++) sum += C[k][j]*B[i][k];
			C[i][j] = ( Y[i][j] - sum ) / p[i];
			C[j][i] = C[i][j];
		}
	}
	free_dmatrix(Y,1,n,1,n);
	free_dvector(v,1,n);
	return;
}


/* CHOLESKY_VEC  -  
 * Extracts eigen-vectors of a generalized eigen-problem, A x =  d B x,
 * from the eigen-vectors of the associated symmetric eigen-problem,
 * C(L'x)=d(L'x).
 * L[1..n][1..n] is the Cholesky factor of B[1..n][1..n].	
 * B is positive-definite.
 * On input, B has already been factored by cholesky_dcmp() in cholesky_sym() 
 * and L'x is x.  
 * On output L'x is replaced by x. Dimensions:  A, B, and C are [1..n][1..n]
 * L[j][i]=L[i][j]=B[i][j] i>j	L[i][i] = p[i]		B[i][j] = B[i][j] i<=j
 * ---------------------------------------------------------------------------
 */
void cholesky_vec ( double **B, int n, double **x, double *p )
{
	float	sum;
	int	i, j, k;

	for (j=1; j<=n; j++) {
	    for (i=0; i<n; i++) {
		sum = 0.0;
		if (i > 1) for (k=i-1; k >=0; k--) sum += B[n-k][n-i]*x[n-k][j];
		x[n-i][j] = ( x[n-i][j] -  sum ) / p[n-i];
	    }
	}
	return;
}


/* CHOLESKY_MPROVE  -  
 * Improves a solution vector x[1..n] of the linear set of equations 
 * [A]{x} = {b}.  The matrix A[1..n][1..n], and the vectors b[1..n] and x[1..n] 
 * are input, as is the dimension n.   The matrix [A] is the Cholesky - 
 * decomposition of the original system matrix, as returned by cholesky_dcmp(). 
 * Also input is the diagonal vector, {p}, of the triangular Cholesky factor,
 * [L].
 * On output, only {x} is modified to an improved set of values.  
 * Numerical Recipes In C, 2nd ed. Cambridge Univ Press, 1992, Ch. 2,  p55ff 
 * ---------------------------------------------------------------------------
 */
void cholesky_mprove (
 double **A, int n,
 double *p, double *b, double *x, double *err )
{
	float	sdp;		/* accumulate the r.h.s. in float precision */	
	double	*r;		/* the residual error			*/
	int	j,i, pd;

	r = dvector(1,n);

	for (i=1;i<=n;i++) {		/* calculate the r.h.s. of 	*/
		sdp = -b[i];		/* [A]{r} = [A]{x+r} - {b}	*/
		for (j=1;j<=n;j++) {	/* A in upper triangle only	*/
			if ( i <= j )   sdp += A[i][j] * x[j];
			else		sdp += A[j][i] * x[j];
		}
		r[i] = sdp;
	}
	*err = 0.0;			// the RMS error of the solution
	cholesky_dcmp ( A, n, p, r, r, 0, 1, &pd ); // solve for the error term
	for (i=1;i<=n;i++) {
		x[i] -= r[i];		// subtract it from the old solution 
		*err += r[i]*r[i];
	}
	*err = sqrt ( *err / (float) n );
	free_dvector(r,1,n);
	return;
}

