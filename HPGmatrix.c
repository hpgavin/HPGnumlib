/*   HPGmatrix.c
 *
 *    Matrix math functions
 *
 *	 Author:  Henri P. Gavin (hpgavin), h p gavin ~at~ duke ~dot~ e d v
 *	Company:  Duke Univ.
 *
 * ==========================================================================
 */

/*  Copyright (C) 2012 Henri P. Gavin
 
    HPGmatrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. 
    
    HPGmatrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with HPGmatrix.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdio.h>

#include "HPGmatrix.h"
#include "NRutil.h"

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define TOL 1e-9	/* machine tolerance for finding singular values */

#define sind(x) (sin(fmod((x),360) * M_PI / 180))
#define cosd(x) (cos(fmod((x),360) * M_PI / 180))


/* 
 * GAUSSJ
 * Linear equation solution by Gauss-Jordan elimination, [A][X]=[B] above. A[1..n][1..n]
 * is the input matrix. B[1..n][1..m] is input containing the m right-hand side vectors. On
 * output, a is replaced by its matrix inverse, and B is replaced by the corresponding set of solution
 * vectors.
 */
void gaussj(float **A, int n, float **B, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol=1,irow=1,j,k,l,ll;
	float big,dum,pivinv,temp;

/* The integer arrays ipiv, indxr, and indxc are used for bookkeeping on the pivoting. */

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;

/*  This is the main loop over the columns to be reduced. */

	for (i=1;i<=n;i++) {
		big=0.0;

/*  This is the outer loop for the search for a pivot element. */

		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(A[j][k]) >= big) {
							big=fabs(A[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) NRerror("gaussj: Singular Matrix-1");
				}
		++(ipiv[icol]);

/* We now have the pivot element, so we interchange rows, if needed, to put the pivot
 * element on the diagonal. The columns are not physically interchanged, only relabeled:
 * indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
 * indxr[i] is the row in which that pivot element was originally located. If indxr[i] =
 * indxc[i] there is an implied column interchange. With this form of bookkeeping, the
 * solution b's will end up in the correct order, and the inverse matrix will be scrambled
 * by columns.
 */

		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(A[irow][l],A[icol][l])
			for (l=1;l<=m;l++) SWAP(B[irow][l],B[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;

/* We are now ready to divide the pivot row by the by the pivot element, located at irow,icol */

		if (A[icol][icol] == 0.0) NRerror("gaussj: Singular Matrix-2");
		pivinv=1.0/A[icol][icol];
		A[icol][icol]=1.0;
		for (l=1;l<=n;l++) A[icol][l] *= pivinv;
		for (l=1;l<=m;l++) B[icol][l] *= pivinv;

/*  Next, we reduce the rows ... except for the pivot one, of course. */

		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=A[ll][icol];
				A[ll][icol]=0.0;
				for (l=1;l<=n;l++) A[ll][l] -= A[icol][l]*dum;
				for (l=1;l<=m;l++) B[ll][l] -= B[icol][l]*dum;
			}
	}

/* This is the end of the main loop over columns of the reduction. It only remains to unscram-
 * ble the solution in view of the column interchanges. We do this by interchanging pairs of
 * columns in the reverse order that the permutation was built up.
 */

	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(A[k][indxr[l]],A[k][indxc[l]]);
	}
/*  And we are done. */
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
#undef SWAP


/*
 * LU_DCMP  
 * Solves [A]{x} = {b} simply and efficiently by performing an 
 * LU - decomposition of [A].  No pivoting is performed. 
 * [A] is a diagonally dominant matrix of dimension [1..n][1..n]. 
 * {b} is a r.h.s. vector of dimension [1..n].
 * {b} is updated using [LU] and then back-substitution is done to obtain {x}.  
 * {b} is replaced by {x} and [A] is replaced by the LU - reduction of itself.
 *
 *  usage:  double **A, *b;
 *	  int   n, reduce, solve, pd; 
 *	  lu_dcmp ( A, n, b, reduce, solve, &pd );		     5may98
 */
void lu_dcmp (
	double **A,	// the system matrix, and it's LU- reduction
	int n,		// the dimension of the matrix
	double *b, 	// the right hand side vector, and the solution vector 
	int reduce,	// 1: do a forward reduction; 0: don't do the reduction 
	int solve,	// 1: do a back substitution for {x};  0: do no bk-sub'n
	int *pd		// 1: positive diagonal  and  successful LU decomp'n	
){
	double	pivot;		/* a diagonal element of [A]		*/
	int	i, j, k;

	*pd = 1;
	if ( reduce ) {			/* forward reduction of [A]	*/

	    for (k=1; k <= n; k++) {
		if ( 0.0 == (pivot = A[k][k]) ) {
		    fprintf(stderr," lu_dcmp: zero found on the diagonal\n");
		    fprintf(stderr," A[%d][%d] = %11.4e\n", k, k, A[k][k] );
		    *pd = 0;
		    return;
		}
		for (i = k+1; i <= n; i++) {
		    A[i][k] /= pivot;
		    for (j=k+1; j <= n; j++)	A[i][j] -= A[i][k]*A[k][j];
		}
	    }
	}		/* the forward reduction of [A] is now complete	*/

	if ( solve ) {		/* back substitution to solve for {x}	*/

	    /* {b} is run through the same forward reduction as was [A]	*/

	    for (k=1; k <= n; k++)
		for (i=k+1; i <= n; i++)	b[i] -= A[i][k]*b[k];

	    /* now back substitution is conducted on {b};  [A] is preserved */

	    for (j=n; j >= 2; j--)
		for (i=1; i <= j-1; i++)	b[i] -= b[j]*A[i][j]/A[j][j];

	    /* finally we solve for the {x} vector			*/

	    for (i=1; i<=n; i++)		b[i] /= A[i][i];
	} 

	/* TAA DAAAAAAAA! {b} is now {x} and is ready to be returned	*/

	return;
}


/*
 * LDL_DCMP  -  Solves [A]{x} = {b} simply and efficiently by performing an 
 * L D L' - decomposition of [A].  No pivoting is performed.  
 * [A] is a symmetric diagonally-dominant matrix of dimension [1..n][1..n]. 
 * {b} is a r.h.s. vector of dimension [1..n].
 * {b} is updated using L D L' and then back-substitution is done to obtain {x}
 * {b} is returned unchanged.  ldl_dcmp(A,n,d,x,x,1,1,&pd) is valid.  
 *     The lower triangle of [A] is replaced by the lower triangle L of the 
 *     L D L' reduction.  The diagonal of D is returned in the vector {d}
 *
 * usage: double **A, *d, *b, *x;
 *	int   n, reduce, solve, pd;
 *	ldl_dcmp ( A, n, d, b, x, reduce, solve, &pd );
 *
 * H.P. Gavin, Civil Engineering, Duke University, hpgavin@duke.edu  9 Oct 2001
 * Bathe, Finite Element Procecures in Engineering Analysis, Prentice Hall, 1982
 */
void ldl_dcmp (
	double **A,	// the system matrix, and L of the L D L' decomp.
	int n,		// the dimension of the matrix	
	double *d,	// diagonal of D in the  L D L' - decomp'n
	double *b,	// the right hand side vector	
	double *x,	// the solution vector	
	int reduce,	// 1: do a forward reduction of A; 0: don't
	int solve,	// 1: do a back substitution for {x}; 0: don't
	int *pd		// 1: definite matrix and successful L D L' decomp'n
){
	int	i, j, k, m;
	*pd = 0;	/* number of negative elements on the diagonal of D */

	if ( reduce ) {		/* forward column-wise reduction of [A]	*/

	    for (j=1; j<=n; j++) {

	    	for (m=1, i=1; i < j; i++) 	/* scan the sky-line	*/
			if ( A[i][j] == 0.0 ) ++m; else	break;
		
		for (i=m; i < j; i++) {
			A[j][i] = A[i][j];
			for (k=m; k < i; k++) A[j][i] -= A[j][k]*A[i][k];
		}

		d[j] = A[j][j];
	    	for (i=m; i < j; i++)	d[j] -= A[j][i]*A[j][i]/d[i];

	    	for (i=m; i < j; i++)	A[j][i] /= d[i];
		if ( d[j] == 0.0 ) {
		    fprintf(stderr," ldl_dcmp(): zero found on diagonal ...\n");
		    fprintf(stderr," d[%d] = %11.4e\n", j, d[j] );
		    return;
		}
		if ( d[j] < 0.0 ) (*pd)--;
	    }
	}		/* the forward reduction of [A] is now complete	*/

	if ( solve ) {		/* back substitution to solve for {x}   */

		/* {x} is run through the same forward reduction as was [A] */

	    for (i=1; i <= n; i++) {
		x[i] = b[i];
		for (j=1; j < i; j++)		x[i] -= A[i][j]*x[j];
	    }

	    for (i=1; i <= n; i++)		x[i] /= d[i];

	    /* now back substitution is conducted on {x};  [A] is preserved */

	    for (i=n; i > 1; i--) 
		for (j=1; j < i; j++)		x[j] -= A[i][j]*x[i];

	} 
	return;
}


/*
 * LDL_MPROVE  Improves a solution vector x[1..n] of the linear set of equations
 * [A]{x} = {b}.  The matrix A[1..n][1..n], and the vectors b[1..n] and x[1..n]
 * are input, as is the dimension n.   The matrix [A] is the L D L'
 * decomposition of the original system matrix, as returned by ldl_dcmp().
 * Also input is the diagonal vector, {d} of [D] of the L D L' decompositon.
 * On output, only {x} is modified to an improved set of values.
 *
 * usage: double **A, *d, *b, *x, rms_resid;
 * 	int   n, ok;
 * 	ldl_mprove ( A, n, d, b, x, &rms_resid, &ok );
 * 
 * H.P. Gavin, Civil Engineering, Duke University, hpgavin@duke.edu  4 May 2001
 */
void ldl_mprove(
	double **A, int n, double *d, double *b, double *x, 
	double *rms_resid, int *ok
){
	double  sdp;		/* accumulate the r.h.s. in double precision */
	double  *resid,		/* the residual error		  	*/
		rms_resid_new=0, /* the RMS error of the mprvd solution	*/
		*dvector();	/* allocate memory for a vector	of doubles */
	int	j,i, pd;
	void	ldl_dcmp(),
		free_dvector();

	resid = dvector(1,n);

	for (i=1;i<=n;i++) {		/* calculate the r.h.s. of      */
		sdp = b[i];		/* [A]{r} = {b} - [A]{x+r}      */
		for (j=1;j<=n;j++) {	/* A in upper triangle only     */
			if ( i <= j )   sdp -= A[i][j] * x[j];
			else		sdp -= A[j][i] * x[j];
		}
		resid[i] = sdp;
	}

	/* solve for the error term */
	ldl_dcmp ( A, n, d, resid, resid, 0, 1, &pd );


	for (i=1;i<=n;i++)	rms_resid_new += resid[i]*resid[i];

	rms_resid_new = sqrt ( rms_resid_new / (double) n );

	*ok = 0;
	if ( rms_resid_new / *rms_resid < 0.90 ) { 	/* good improvement */
				/* subtract the error from the old solution */
		for (i=1;i<=n;i++)	x[i] += resid[i];
		*rms_resid = rms_resid_new;
		*ok = 1;	/* the solution has improved		*/
	}

	free_dvector(resid,1,n);
	return;
}


/*
 * LDL_DCMP_PM  -  Solves partitioned matrix equations
 *
 *           [A_qq]{x_q} + [A_qr]{x_r} = {b_q}
 *           [A_rq]{x_q} + [A_rr]{x_r} = {b_r}+{c_r}
 *           where {b_q}, {b_r}, and {x_r} are known and 
 *           where {x_q} and {c_r} are unknown
 * 
 * via L D L' - decomposition of [A_qq].  No pivoting is performed.  
 * [A] is a symmetric diagonally-dominant matrix of dimension [1..n][1..n]. 
 * {b} is a r.h.s. vector of dimension [1..n].
 * {b} is updated using L D L' and then back-substitution is done to obtain {x}
 * {b_q} and {b_r}  are returned unchanged. 
 * {c_r} is returned as a vector of [1..n] with {c_q}=0.
 * {q} is a vector of the indexes of known values {b_q}
 * {r} is a vector of the indexes of known values {x_r}
 *     The lower triangle of [A_qq] is replaced by the lower triangle L of its 
 *     L D L' reduction.  The diagonal of D is returned in the vector {d}
 *
 * usage: double **A, *d, *b, *x;
 *	int   n, reduce, solve, pd;
 *	ldl_dcmp_pm ( A, n, d, b, x, c, q, r, reduce, solve, &pd );
 *
 * H.P. Gavin, Civil Engineering, Duke University, hpgavin@duke.edu
 * Bathe, Finite Element Procecures in Engineering Analysis, Prentice Hall, 1982
 * 2014-05-14 
 */
void ldl_dcmp_pm (
	double **A,	/**< the system matrix, and L of the L D L' decomp.*/
	int n,		/**< the dimension of the matrix		*/
	double *d,	/**< diagonal of D in the  L D L' - decomp'n    */
	double *b,	/**< the right hand side vector			*/
	double *x,	/**< part of the solution vector		*/
	double *c,	/**< the part of the solution vector in the rhs */
	int *q,		/**< q[j]=1 if  b[j] is known; q[j]=0 otherwise	*/
	int *r,		/**< r[j]=1 if  x[j] is known; r[j]=0 otherwise	*/
	int reduce,	/**< 1: do a forward reduction of A; 0: don't   */
	int solve,	/**< 1: do a back substitution for {x}; 0: don't */
	int *pd		/**< 1: definite matrix and successful L D L' decomp'n*/
){
	int	i, j, k, m;
	*pd = 0;	/* number of negative elements on the diagonal of D */

	if ( reduce ) {		/* forward column-wise reduction of [A]	*/

	    for (j=1; j<=n; j++) {

	      d[j] = 0.0;

	      if ( q[j] ) { /* reduce column j, except where q[i]==0	*/	

		for (m=1, i=1; i < j; i++) 	/* scan the sky-line	*/
			if ( A[i][j] == 0.0 ) ++m; else	break;
		
		for (i=m; i < j; i++) {
		    if ( q[i] ) {
			A[j][i] = A[i][j];
			for (k=m; k < i; k++)
				if ( q[k] )
					A[j][i] -= A[j][k]*A[i][k];
		    }
		}

		d[j] = A[j][j];
	    	for (i=m; i < j; i++) if ( q[i] ) d[j] -= A[j][i]*A[j][i]/d[i];
	    	for (i=m; i < j; i++) if ( q[i] ) A[j][i] /= d[i];

		if ( d[j] == 0.0 ) {
		 fprintf(stderr," ldl_dcmp_pm(): zero found on diagonal ...\n");
		 fprintf(stderr," d[%d] = %11.4e\n", j, d[j] );
		 return;
		}
		if ( d[j] < 0.0 ) (*pd)--;
	      }
	    }
		
	}		/* the forward reduction of [A] is now complete	*/

	if ( solve ) {		/* back substitution to solve for {x}   */

	    for (i=1; i <= n; i++) {
		if ( q[i] ) {
			x[i] = b[i];
			for (j=1; j<= n; j++) if ( r[j] ) x[i] -= A[i][j]*x[j];
		}
	    }

		/* {x} is run through the same forward reduction as was [A] */
	    for (i=1; i <= n; i++) 
		if ( q[i] ) for (j=1;j<i;j++) if ( q[j] ) x[i] -= A[i][j]*x[j];

	    for (i=1; i <= n; i++)	if ( q[i] )	x[i] /= d[i];

	    /* now back substitution is conducted on {x};  [A] is preserved */

	    for (i=n; i > 1; i--) 
		if ( q[i] )
			for (j=1; j < i; j++)
				if ( q[j] )
					x[j] -= A[i][j]*x[i];
	
	    /* finally, evaluate c_r	*/

	    for (i=1; i<=n; i++) {
		c[i] = 0.0;
		if ( r[i] ) {
			c[i] = -b[i]; // changed from 0.0 to -b[i]; 2014-05-14
			for (j=1; j<=n; j++)	c[i] += A[i][j]*x[j];
		}
	    }
				
	}
	return;
}


/*
 * LDL_MPROVE_PM 
 * Improves a solution vector x[1..n] of the partitioned set of linear equations
 *           [A_qq]{x_q} + [A_qr]{x_r} = {b_q}
 *           [A_rq]{x_q} + [A_rr]{x_r} = {b_r}+{c_r}
 *           where {b_q}, {b_r}, and {x_r} are known and 
 *           where {x_q} and {c_r} are unknown
 * by reducing the residual r_q
 *           A_qq r_q = {b_q} - [A_qq]{x_q+r_q} + [A_qr]{x_r} 
 * The matrix A[1..n][1..n], and the vectors b[1..n] and x[1..n]
 * are input, as is the dimension n.   The matrix [A_qq] is the L D L'
 * decomposition of the original system matrix, as returned by ldl_dcmp_pm().
 * Also input is the diagonal vector, {d} of [D] of the L D L' decompositon.
 * On output, only {x} and {c} are modified to an improved set of values.
 * The partial right-hand-side vectors, {b_q} and {b_r}, are returned unchanged.
 * Further, the calculations in ldl_mprove_pm do not involve b_r.  
 * 
 * usage: double **A, *d, *b, *x, rms_resid;
 * 	int   n, ok, *q, *r;
 *	ldl_mprove_pm ( A, n, d, b, x, q, r, &rms_resid, &ok );
 *
 * H.P. Gavin, Civil Engineering, Duke University, hpgavin@duke.edu 
 * 2001-05-01, 2014-05-14 
 */
void ldl_mprove_pm (
	double **A,	/**< the system matrix, and L of the L D L' decomp.*/
	int n,		/**< the dimension of the matrix		*/
	double *d,	/**< diagonal of D in the  L D L' - decomp'n    */
	double *b,	/**< the right hand side vector			*/
	double *x,	/**< part of the solution vector		*/
	double *c,	/**< the part of the solution vector in the rhs */
	int *q,		/**< q[j]=1 if  b[j] is known; q[j]=0 otherwise	*/
	int *r,		/**< r[j]=1 if  x[j] is known; r[j]=0 otherwise	*/
	double *rms_resid, /**< root-mean-square of residual error	*/
	int *ok		/**< 1: >10% reduction in rms_resid; 0: not	*/
){
	double  sdp;		// accumulate the r.h.s. in double precision
	double  *dx,		// the residual error
		*dc,		// update to partial r.h.s. vector, c
		rms_resid_new=0.0, // the RMS error of the mprvd solution
		*dvector();	// allocate memory for a vector	of doubles
	int	j,i, pd;
	void	ldl_dcmp(),
		free_dvector();

	dx  = dvector(1,n);
	dc  = dvector(1,n);

	for (i=1;i<=n;i++)	dx[i] = 0.0;

	// calculate the r.h.s. of ...
	//  [A_qq]{dx_q} = {b_q} - [A_qq]*{x_q} - [A_qr]*{x_r}      
	//  {dx_r} is left unchanged at 0.0;
	for (i=1;i<=n;i++) {	
	    if ( q[i] ) {
		sdp = b[i];	
		for (j=1;j<=n;j++) {
			if ( q[j] ) {	// A_qq in upper triangle only
				if ( i <= j )   sdp -= A[i][j] * x[j];
				else		sdp -= A[j][i] * x[j];
			}
		}
		for (j=1;j<=n;j++) if ( r[j] )	sdp -= A[i][j] * x[j];
		dx[i] = sdp;
	    } // else dx[i] = 0.0; // x[i];
	}

	// solve for the residual error term, A is already factored
	ldl_dcmp_pm ( A, n, d, dx, dx, dc, q,r, 0, 1, &pd );

	for (i=1;i<=n;i++) if ( q[i] )	rms_resid_new += dx[i]*dx[i];

	rms_resid_new = sqrt ( rms_resid_new / (double) n );

	*ok = 0;
	if ( rms_resid_new / *rms_resid < 0.90 ) { /*  enough improvement    */
		for (i=1;i<=n;i++) {	/*  update the solution 2014-05-14   */
		    	if ( q[i] )	x[i] += dx[i];
			if ( r[i] )	c[i] += dc[i];
		}
		*rms_resid = rms_resid_new;	/* return the new residual   */
		*ok = 1;			/* the solution has improved */
	}

	free_dvector(dx,1,n);
	free_dvector(dc,1,n);
	return;
}


/*
 * PSB_UPDATE
 * Update secant stiffness matrix via the Powell-Symmetric-Broyden update eqn.
 *
 *       B = B - (f*d' + d*f') / (d' * d) + f'*d * d*d' / (d' * d)^2 ;
 *
 * H.P. Gavin, Civil Engineering, Duke University, hpgavin@duke.edu  24 Oct 2012
 */
void PSB_update ( 
	double **B,	/**< secant stiffness matrix            */
	double *f,	/**< out-of-balance force vector        */
	double *d,	/**< incremental displacement vector    */
	int n )		/**< matrix dimension is n-by-n         */
{
	int	i, j;
	double	dtd = 0.0, ftd = 0.0, dtd2 = 0.0;

	for (i=1; i<=n; i++)	dtd += d[i]*d[i];
	dtd2 = dtd*dtd;

	for (i=1; i<=n; i++)	ftd += f[i]*d[i];

	for (i=1; i<=n; i++)	/*  update upper triangle of B[i][j] */
	    for (j=i; j<=n; j++) 
		B[i][j] -= ( (f[i]*d[j] + f[j]*d[i])/dtd - ftd*d[i]*d[j]/dtd2 );
	
}


/*
 * PSEUDO_INV - calculate the pseudo-inverse of A ,
 * 	     Ai = inv ( A'*A + beta * trace(A'*A) * I ) * A' 
 *	     beta is a regularization factor, which should be small (1e-10)
 *	     A is m by n      Ai is m by n			      8oct01
 */
void pseudo_inv(
	double **A, double **Ai, int n, int m, double beta, int verbose
){
	double	*diag, *b, *x, **AtA, **AtAi, tmp, tr_AtA=0.0,
		*dvector(), **dmatrix(), error;
	int     i,j,k, ok; 
	void	ldl_dcmp(), ldl_mprove(), free_dvector(), free_dmatrix();

	diag = dvector(1,n);
	b    = dvector(1,n);
	x    = dvector(1,n);
	AtA  = dmatrix(1,n,1,n);
	AtAi = dmatrix(1,n,1,n);

	if (beta>1) fprintf(stderr," pseudo_inv: warning beta = %lf\n", beta);

	for (i=1; i<=n; i++) {
		diag[i] = x[i] = b[i] = 0.0;
		for (j=i; j<=n; j++)    AtA[i][j] = AtA[j][i] = 0.0;
	}

	for (i=1; i<=n; i++) {			/* compute A' * A */
		for (j=1; j<=n; j++) {
			tmp = 0.0;
			for (k=1; k<=m; k++) tmp += A[k][i] * A[k][j];
			AtA[i][j] = tmp;
		}
	}
	for (i=1; i<=n; i++)			    /* make symmetric */
		for (j=i; j<=n; j++)
			AtA[i][j]=AtA[j][i] = 0.5*(AtA[i][j] + AtA[j][i]);

	tr_AtA = 0.0;
	for (i=1; i<=n; i++) tr_AtA += AtA[i][i];       /* trace of AtA */
	for (i=1; i<=n; i++) AtA[i][i] += beta*tr_AtA;	/* add beta I */

	ldl_dcmp ( AtA, n, diag, b, x, 1, 0, &ok );	/*  L D L'  decomp */

	for (j=1; j<=n; j++) {				/* compute inv(AtA) */

		for (k=1; k<=n; k++)  b[k] = 0.0;
		b[j] = 1.0;
		ldl_dcmp( AtA, n, diag, b, x, 0, 1, &ok ); /* L D L' bksbtn */

		if ( verbose )
		 fprintf(stdout,"  RMS matrix error:"); /*improve the solution*/
		error = 1.0; ok = 1;
		do {
			ldl_mprove ( AtA, n, diag, b, x, &error, &ok );
			if ( verbose ) fprintf(stdout,"%9.2e", error );
		} while ( ok );
		if ( verbose ) fprintf(stdout,"\n");

		for (k=1; k<=n; k++)  AtAi[k][j] = x[k];  /* save inv(AtA) */
	}

	for (i=1; i<=n; i++)			    /* make symmetric */
		for (j=i; j<=n; j++)
			AtAi[i][j]=AtAi[j][i] = 0.5*(AtAi[i][j] + AtAi[j][i]);

	for (i=1; i<=n; i++) {			/* compute inv(A'*A)*A'	*/
		for (j=1; j<=m; j++) {
			tmp = 0.0;
			for (k=1; k<=n; k++)    tmp += AtAi[i][k]*A[j][k];
			Ai[i][j] = tmp;
		}
	}

	free_dmatrix (AtAi, 1,n,1,n);
	free_dmatrix (AtA,  1,n,1,n);
	free_dvector (x, 1,n);
	free_dvector (b, 1,n);
	free_dvector (diag, 1,n);

	return;
}


/* 
 * PRODABj -  matrix-matrix multiplication for symmetric A	      27apr01
 *		 u = A * B(:,j)
 */
void prodABj ( double **A, double **B, double *u, int n, int j )
{
	int     i, k;

	for (i=1; i<=n; i++)    u[i] = 0.0;

	for (i=1; i<=n; i++) {
		for (k=1; k<=n; k++) {
			if ( i <= k )   u[i] += A[i][k]*B[k][j];
			else	    u[i] += A[k][i]*B[k][j];
		}
	}
	return;
}


/* 
 * prodAB - matrix-matrix multiplication      C = A * B			27apr01
 */
void prodAB ( double **A, double **B, double **C, int I, int J, int K )
{
	int     i, j, k;

	for (i=1; i<=I; i++)
		for (k=1; k<=K; k++)
			C[i][k] = 0.0;

	for (i=1; i<=I; i++)
		for (k=1; k<=K; k++)
			for (j=1; j<=J; j++)
				C[i][k] += A[i][j]*B[j][k];
	return;
}


/*
 * INVAB  -  calculate product inv(A) * B  
 *	 A is n by n      B is n by m				    6jun07
 */
void invAB( double **A, double **B, int n, int m, double **AiB, int *ok, int verbose)
{
	double  *diag, *b, *x, error;
	int     i,j,k;

	diag = dvector(1,n);
	x    = dvector(1,n);
	b    = dvector(1,n);

	for (i=1; i<=n; i++) diag[i] = x[i] = 0.0;

	ldl_dcmp( A, n, diag, b, x, 1, 0, ok );	 /*   L D L'  decomp */
	if ( *ok < 0 ) {
		fprintf(stderr," Make sure that all six");
		fprintf(stderr," rigid body translations are restrained!\n");
	}

	for (j=1; j<=m; j++) {

		for (k=1; k<=n; k++)  b[k] = B[k][j];
		ldl_dcmp( A, n, diag, b, x, 0, 1, ok ); /*   L D L'  bksbtn */

		if ( verbose ) fprintf(stdout,"    LDL' RMS matrix precision:");
		error = *ok = 1;
		do {				    /* improve the solution*/
			ldl_mprove ( A, n, diag, b, x, &error, ok );
			if ( verbose ) fprintf(stdout,"%9.2e", error );
		} while ( *ok );
		if ( verbose ) fprintf(stdout,"\n");

		for (i=1; i<=n; i++)    AiB[i][j] = x[i];
	}

	free_dvector(diag,1,n);
	free_dvector(x,1,n);
	free_dvector(b,1,n);
}


/*
 * XTinvAY  -  calculate quadratic form with inverse matrix   X' * inv(A) * Y   
 *	   A is n by n    X is n by m     Y is n by m		    15sep01
 */
void xtinvAy(
	double **X, double **A, double **Y, int n, int m, double **Ac, int verbose
){
	double  *diag, *x, *y, error;
	int     i,j,k, ok;

	diag = dvector(1,n);
	x    = dvector(1,n);
	y    = dvector(1,n);

	for (i=1; i<=n; i++) diag[i] = x[i] = 0.0;

	ldl_dcmp( A, n, diag, y, x, 1, 0, &ok );	/*   L D L'  decomp */

	for (j=1; j<=m; j++) {

		for (k=1; k<=n; k++)  y[k] = Y[k][j];
		ldl_dcmp( A, n, diag, y, x, 0, 1, &ok ); /*   L D L'  bksbtn */

		if ( verbose ) fprintf(stdout,"    LDL' RMS matrix precision:");
		error = ok = 1;
		do {				    /* improve the solution*/
			ldl_mprove ( A, n, diag, y, x, &error, &ok );
			if ( verbose ) fprintf(stdout,"%9.2e", error );
		} while ( ok );
		if ( verbose ) fprintf(stdout,"\n");

		for (i=1; i<=m; i++) {
			Ac[i][j] = 0.0;
			for (k=1; k<=n; k++) Ac[i][j] += X[k][i] * x[k];
		}
	}

	free_dvector(diag,1,n);
	free_dvector(x,1,n);
	free_dvector(y,1,n);
}


/*  COORD_XFRM - coordinate transform of a matrix of column 2-vectors
 * 
 * Rr  = [ cosd(theta) -sind(theta) ; sind(theta) cosd(theta) ]*[ Rx ; Ry ];
 */
void coord_xfrm ( float **Rr, float **R, float theta, int n )
{
	float	R1, R2;
	int	i;

	for ( i = 1; i<=n; i++ ) {
		R1 =  cosd(theta)*R[1][i] - sind(theta)*R[2][i];
		R2 =  sind(theta)*R[1][i] + cosd(theta)*R[2][i];
		Rr[1][i] = R1;
		Rr[2][i] = R2;
	}
	return;
}


/*
 * xtAx - carry out matrix-matrix-matrix multiplication for symmetric A  7nov02
 *	 C = X' A X     C is J by J      X is N by J     A is N by N      
 */ 
void xtAx(double **A, double **X, double **C, int N, int J)
{

	double  **AX;
	int     i,j,k;

	AX = dmatrix(1,N,1,J);

	for (i=1; i<=J; i++)    for (j=1; j<=J; j++)     C[i][j] = 0.0;
	for (i=1; i<=N; i++)    for (j=1; j<=J; j++)    AX[i][j] = 0.0;
		
	for (i=1; i<=N; i++) {	  /*  use upper triangle of A */
		for (j=1; j<=J; j++) {
			for (k=1; k<=N; k++) {
				if ( i <= k )   AX[i][j] += A[i][k] * X[k][j];
				else	    AX[i][j] += A[k][i] * X[k][j];
			}
		} 
	}       
		
	for (i=1; i<=J; i++)    
		for (j=1; j<=J; j++)
			for (k=1; k<=N; k++) 
				C[i][j] += X[k][i] * AX[k][j];
			
	for (i=1; i<=J; i++)	    /*  make  C  symmetric */
		for (j=i; j<=J; j++)
			C[i][j] = C[j][i] = 0.5 * ( C[i][j] + C[j][i] );

	free_dmatrix(AX,1,N,1,J);
	return;
}


/* 
 * xtAy - carry out vector-matrix-vector multiplication for symmetric A  7apr94
 */
double xtAy(double *x, double **A, double *y, int n, double *d)
{

	double  xtAy = 0.0;
	int     i,j;

	for (i=1; i<=n; i++) {				/*  d = A y  */
		d[i]  = 0.0;    
		for (j=1; j<=n; j++) {	  //  A in upper triangle only 
			if ( i <= j )   d[i] += A[i][j] * y[j];
			else	    d[i] += A[j][i] * y[j];
		}
	}
	for (i=1; i<=n; i++)    xtAy += x[i] * d[i];	/*  xAy = x' A y  */
	return ( xtAy );
}


/*
 * invAXinvA -  calculate quadratic form with inverse matrix 
 *	   replace X with inv(A) * X * inv(A) 
 *	   A is n by n and symmetric   X is n by n and symmetric    15sep01
 */ 
void invAXinvA ( double **A, double **X, int n, int verbose )
{
	double  *diag, *b, *x, **Ai, **XAi, Aij, error;
	int     i,j,k, ok;

	diag = dvector(1,n);
	x    = dvector(1,n);
	b    = dvector(1,n);
	Ai   = dmatrix(1,n,1,n);
	XAi  = dmatrix(1,n,1,n);

	for (i=1; i<=n; i++) {
		diag[i] = x[i] = b[i] = 0.0;
		for (j=1; j<=n; j++)    XAi[i][j] = Ai[i][j] = 0.0;
	}

	ldl_dcmp ( A, n, diag, b, x, 1, 0, &ok );       /*   L D L'  decomp */

	for (j=1; j<=n; j++) {			  /*  compute inv(A)  */

		for (k=1; k<=n; k++)  b[k] = 0.0;
		b[j] = 1.0;
		ldl_dcmp( A, n, diag, b, x, 0, 1, &ok ); /*   L D L'  bksbtn */

		if ( verbose ) fprintf(stdout,"    LDL' RMS matrix precision:");
		error = ok = 1;
		do {				    /* improve the solution*/
			ldl_mprove ( A, n, diag, b, x, &error, &ok );
			if ( verbose ) fprintf(stdout,"%9.2e", error );
		} while ( ok );
		if ( verbose ) fprintf(stdout,"\n");

		for (k=1; k<=n; k++)  Ai[j][k] = x[k];  /*  save inv(A) */
	}

	for (i=1; i<=n; i++)			    /*  make symmetric */
		for (j=i; j<=n; j++)
			Ai[i][j] = Ai[j][i] = 0.5 * ( Ai[i][j] + Ai[j][i] );

	for (i=1; i<=n; i++) {		  /*  compute X * inv(A)   */
		for (j=1; j<=n; j++) {
			Aij = 0.0;
			for (k=1; k<=n; k++)    Aij += X[i][k]*Ai[k][j];
			XAi[i][j] = Aij;
		}
	}

	for (i=1; i<=n; i++) {	  /*  compute inv(A) * X * inv(A)  */
		for (j=1; j<=n; j++) {
			Aij = 0.0;
			for (k=1; k<=n; k++)    Aij += Ai[i][k] * XAi[k][j];
			X[i][j] = Aij;
		}
	}
	for (i=1; i<=n; i++)			    /*  make symmetric */
		for (j=i; j<=n; j++)
			X[i][j] = X[j][i] = 0.5 * ( X[i][j] + X[j][i] );

	free_dvector ( diag, 1,n );
	free_dvector ( x, 1,n );
	free_dvector ( b, 1,n );
	free_dmatrix ( Ai, 1,n,1,n );
	free_dmatrix ( XAi, 1,n,1,n );

}


/* 
 *  RELATIVE_NORM -  compute the relative 2-norm between two vectors       26dec01 
 *       compute the relative 2-norm between two vectors N and D
 *	       return  ( sqrt(sum(N[i]*N[i]) / sqrt(D[i]*D[i]) )
 *
 */
double relative_norm( double *N, double *D, int n )
{
	double  nN = 0.0, nD = 0.0;
	int     i;

	for (i=1; i<=n; i++)    nN += ( N[i]*N[i] );
	for (i=1; i<=n; i++)    nD += ( D[i]*D[i] );

	return ( sqrt(nN) / sqrt(nD) );
}


/*
 *  Legendre 
 *  compute matrix of the Legendre polynomials and its first two derivitives
 */
void Legendre( int order, float *t, int n, float **P, float **Pp, float **Ppp )
{
	int k,p;

// save_vector( n, t, "t.dat");

	for (p=1; p <= n; p++) {

	 P[0][p]   =  1.00;
	 P[1][p]   =  t[p];
	 P[2][p]   =  1.50 * t[p]*t[p] - 0.50;
	 P[3][p]   =  2.50 * t[p]*t[p]*t[p] - 1.50 * t[p];

	 Pp[0][p]  =  0.00;
	 Pp[1][p]  =  1.00;
	 Pp[2][p]  =  3.00 * t[p];
	 Pp[3][p]  =  7.50 * t[p]*t[p] - 1.50;

	 Ppp[0][p] =  0.00;
	 Ppp[1][p] =  0.00;
	 Ppp[2][p] =  3.00; 
	 Ppp[3][p] = 15.00 * t[p]; 

	 for ( k=4; k <= order; k++) {

	  P[k][p]   = (2.0-1.0/k)*t[p]*P[k-1][p] - (1.0-1.0/k)*P[k-2][p];

	  Pp[k][p]  = (2.0-1.0/k)*(P[k-1][p] + t[p]*Pp[k-1][p]) - (1.0-1.0/k)*Pp[k-2][p];

	  Ppp[k][p] = (2.0-1.0/k)*(2*Pp[k-1][p] + t[p]*Ppp[k-1][p]) - (1.0-1.0/k)*Ppp[k-2][p];

	 }
	}

	return;
}


/*------------------------------------------------------------------------------
SVDCMP -  singular value decomposition from Numerical Recipes in C 
------------------------------------------------------------------------------*/
static float at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static float maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? \
	(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void svdcmp( float **a, int m, int n, float *w, float **v )
{
	int flag,i,its,j,jj,k,l,nm;
	float c,f,h,s,x,y,z;
	float anorm=0.0,g=0.0,scale=0.0;
	float *rv1;

	if (m < n) fprintf(stderr,"SVDCMP: You must augment A with extra zero rows");
	rv1=vector(1,n);
		/*	Householder  Reduction  To  Bidiagonal  Form */
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++)
							s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++)
							a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++)
							s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++)
							a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
		/*	Accumulation  Of  Right-Hand  Transformations */
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++)
						s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++)
						v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
		/*	Accumulation  Of  Left-Hand  Transformations */
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++)	
						s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++)
						a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
		/*	Diagonalization  Of  The  Bidiagonal  Form  */
	for (k=n;k>=1;k--) {
		for (its=1;its<=300;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					if (fabs(f)+anorm != anorm) {
						g=w[i];
						h=PYTHAG(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						for (j=1;j<=m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 300) fprintf(stderr,"No cnvrgnce in 300 SVDCMP itrs");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=PYTHAG(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}	
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);
	return;
}

/*----------------------------------------------------------------------------
SVBKSB  -  singular value backsubstitution from Numerical Recipes in C
----------------------------------------------------------------------------*/
void svbksb(u,w,v,m,n,b,x)
float **u,w[],**v,b[],x[];
int m,n;
{
	int jj,j,i;
	float s,*tmp,*vector();
	void free_vector();

	tmp=vector(1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	free_vector(tmp,1,n);
	return;
}


	/*		function svdcmp_coefs()			*/
	/* Calculates the least sqares solution to the problem	*/
	/* [A] {x} = {b}.  [A] is m by n with	m >= n.		*/
	/* Uses a singular value decomposition in which the	*/
	/* near zero diagonal elements are zeroed artificially.	*/

void svdcmp_coefs( float **A, int m, int n, float *x, float *b)
{
	float	*w, **v,	/* [A] = [u] [diag w] [v]	*/
		wmin, wmax=0.0,	/* min & max singular values	*/
		**matrix(),*vector();
	int	i,j;
	void	svdcmp(), svbksb(),
		free_matrix(),free_vector();

	w = vector(1,n);
	v = matrix(1,n,1,n);

	svdcmp (A,m,n,w,v);

	for (j=1;j<=n;j++) if (w[j] > wmax) wmax = w[j];
	wmin = wmax*TOL;
	i=0;
	for (j=1;j<=n;j++) if (w[j] < wmin) {w[j]=0.0; ++i;}
	if (i > 0) printf ("  warning: %d singular values\n", i);

	svbksb (A,w,v,m,n,b,x);

	free_vector(w,1,n);
	free_matrix(v,1,n,1,n);
	return;
}

/*-----------------------------------------------------------------------------
LST_SQR_FIT  -  fits a power polyomial to data using SVD
------------------------------------------------------------------------------*
void    lst_sqr_fit (data,bias,points,N)
float   *data, *bias;
int     points, N;
{
        float   t, **time,
                power(),
                **matrix();
        int     i,j;

        void    svdcmp_coefs(),
                free_matrix();

        time = matrix(1,points,1,N);

        for (i=1;i<=points;i++) {
                t = (float) i;
                time[i][1] = 1.0;
                for (j=2;j<=N;j++) time[i][j] = time[i][j-1]*t;
        }

        svdcmp_coefs ( time, points, N, bias, data );

        printf ("  y =%14.6e\n",bias[1]);
        for (j=2;j<=N;j++) printf("    +%14.6e*x^%d\n",bias[j],j-1);

        for (i=1;i<=points;i++) {
                t = (float) i;
                data[i] -= bias[1];
                for (j=2;j<=N;j++) data[i] -= bias[j]*power(t,j-1);
        }

        free_matrix(time,1,points,1,N);
}

                //      Put  The  Bias  Trend  Back

        for (i=1;i<=2*m;i++) {
                data = (float) i;
                y[i] += bias[1];
                for (j=2;j<=N;j++) y[i] += bias[j]*power(data,j-1);
        }

		*/


/*-----------------------------------------------------------------------------
MPROVE  -  Improves a solution vector x[1..n] of the linear set of equations
[A]{x} = {b}.  The matrix A[1..n][1..n], and the vectors b[1..n] and x[1..n]
are input, as is the dimension n.  Also input is the singular value
decomposition of a, as returned by svdcmp.
------------------------------------------------------------------------------*/
void mprove(a,u,s,v,n,b,x,err)
float	**a, **u, *s, **v, *b,*x,*err;
int	n;
{
	int	i,j;
	double	sdp;
	float	*r, *vector();
	void	svbksb(), free_vector();

	*err = 0.0;

	r = vector(1,n);
	for (i=1;i<=n;i++) {
		sdp = -b[i];
		for (j=1;j<=n;j++) sdp += a[i][j]*x[j];
		r[i] = sdp;
	}
	svbksb(u,s,v,n,n,r,r);
	for (i=1;i<=n;i++) {
		x[i] -= r[i];
                *err += r[i]*r[i];
        }
        *err = sqrt ( *err / (float) n );
	free_vector(r,1,n);
}

