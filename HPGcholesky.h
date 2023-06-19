/* HPGcholeskly.h
 *
 * header file for functions in HPGcholesky.c
 *
 * Cholesky decomposition and related routines adapted from ...
 * Press, Numerical Recipes in C, 2nd ed., Cambridge Univ. Pr., 1988, Ch.2    
 */

#include <stdio.h>
#include <math.h>


void cholesky_dcmp (
 double **A,	// the system matrix, and it's Cholesky decomposition
 int n,		// the dimension of the matrix
 double *p,	// the vector of the diagonal of the Cholesky decomp'n
 double *b,	// the right hand side vector
 double *x,	// the solution vector
 int reduce,	// 1: do a forward reduction; 0: don't do the reduction
 int solve,	// 1: do a back substitution; 0: do no back substitution
 int *pd	// 1: positive definite; 0: not positive definite
);
 

void cholesky_correlate (
 float **C,	// lower triangle of Cz is Cholesky factor of correl. matx
 float *Czd,	// diagonal of Cholesky factor of correl. matx
 int n,		// dimension of the correlation matrix
 float *Zu,	// uncorrelated standard normal random variables
 float *Zc,	// correlated standard normal random variables
 int reduce	// 1: Cz is not yet factored, 0: Cz is already factored
);

void cholesky_inv ( double **A, int n, double *p );

void cholesky_sym ( double **A, double **B, int n, double **C, double *p );

void cholesky_vec ( double **B, int n, double **x, double *p );

void cholesky_mprove (
	double **A, int n, double *p, double *b, double *x, double *err );
