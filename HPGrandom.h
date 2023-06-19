/*
 * HPGrandom.h
 *
 * header file for Random.c and related programs
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#define PI       3.14159265358979323846264338327950288419716939937510
#define Sqrt2    1.41421356237309504880168872420969808
#define SqrtPi2  0.88622692545275801364908374167057259


void initialize_random ( unsigned int seed );

double erfinv ( double y ); 

double invnorm ( double u );


void lognormal (
 float  *X,     // the computed log normal random variables
 float  *M,     // the means of the log normal random variables 
 float  *COV,   // the coefficient of variation of the log normal variables
 float  **Cz,   // the lower Cholesky factor of the correl. mtrx 
 float  *Czd,   // the diagonal of the lower Cholesky factor of the correl mtrx
 float  S,      // the number of allowabe standard deviations from mean
 int    n  
);


float inverse_normal ( float p );

double ltqnorm ( double p );

