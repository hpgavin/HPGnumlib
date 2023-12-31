/** @file
	Memory allocation functions from Numerical Recipes in C, by Press,
	Cambridge University Press, 1988
	http://www.nr.com/public-domain.html
*/

#ifndef _NRUTIL_H_
#define _NRUTIL_H_



#include <stdint.h>    // various integer types with specified widths

typedef struct F_COMPLEX {float r,i;} f_complex; // also in complex.h

/*
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
*/

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void NRerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
int32_t *i32vector(long nl, long nh);
uint16_t *u16vector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **subMatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_i32vector(int32_t *v, long nl, long nh);
void free_u16vector(uint16_t *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_subMatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

f_complex *Cvector(int nl, int nh);
f_complex **Cmatrix(int nrl, int nrh, int ncl, int nch);
void free_Cvector(f_complex *v, int nl, int nh);
void free_Cmatrix(f_complex **m, int nrl, int nrh, int ncl, int nch);

void NRerror(char error_text[]);

float  ***D3matrix(int nrl,int nrh, int ncl, int nch, int nzl, int nzh);
double ***D3dmatrix(int nrl,int nrh, int ncl, int nch, int nzl, int nzh);
void   free_D3matrix(float ***m, int nrl, int nrh, int ncl, int nch, int nzl, int nzh);
void free_D3dmatrix(double ***m, int nrl, int nrh, int ncl, int nch, int nzl, int nzh);

void show_vector(float *A, int n);
void show_dvector(double *A, int n );
void show_matrix(float **A, int m, int n );
void show_dmatrix(double **A, int m, int n );

void save_vector(char filename[], float *V, int nl, int nh, const char *mode);
void save_dvector(char filename[], double *V, int nl, int nh, const char *mode);
void save_ivector(char filename[], int *V, int nl, int nh, const char *mode);

void save_matrix ( char filename[], float **A, int ml, int mh, int nl, int nh, int transpose, const char *mode );
void save_dmatrix ( char filename[], double **A, int ml, int mh, int nl, int nh, int transpose, const char *mode );

void save_ut_matrix ( char filename[], float **A, int n, const char *mode );
void save_ut_dmatrix ( char filename[], double **A, int n, const char *mode );



#else /* ANSI */
/* traditional - K&R */

void NRerror();
float *vector();
float **matrix();
float **subMatrix();
float **convert_matrix();
float ***f3tensor();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
unsigned char *cvector();
unsigned long *lvector();
void free_vector();
void free_dvector();
void free_ivector();
void free_cvector();
void free_lvector();
void free_matrix();
void free_subMatrix();
void free_convert_matrix();
void free_dmatrix();
void free_imatrix();
void free_f3tensor();

#endif /* ANSI */

#endif /* _NRUTIL_H */

