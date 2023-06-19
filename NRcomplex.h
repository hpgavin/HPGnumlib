/*			FILE NRcomplex.h				*/
/* complex math arithmetic from Numerical Recipes in C  pp 710-712	*/

#ifndef _NRUTIL_H_
typedef struct F_COMPLEX {float r,i;} f_complex; 
#endif

#ifndef _NRCOMPLEX_H_
#define _NRCOMPLEX_H_

f_complex Cadd( f_complex a, f_complex b);

f_complex Csub( f_complex a, f_complex b);

f_complex Cmul( f_complex a, f_complex b);

f_complex Cdiv( f_complex a, f_complex b);

f_complex Complex( float re, float im);

float Cabs( f_complex z );

float Carg( f_complex z );

f_complex Conjg( f_complex z );

f_complex Csqrt( f_complex z );

f_complex RCmul( float x, f_complex a);

f_complex Cpow( f_complex z, int n);

float Re( f_complex z );

float Im( f_complex z);

#endif
