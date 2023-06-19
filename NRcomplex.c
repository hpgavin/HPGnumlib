/*			FILE NRcomplex.c				*/
/* complex math arithmetic from Numerical Recipes in C  pp 710-712	*/

#include <math.h>
#include <stdio.h>
#include "NRcomplex.h"

// typedef struct F_COMPLEX {float r,i;} f_complex; /* */

f_complex Cadd( f_complex a, f_complex b)
{
	f_complex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

f_complex Csub( f_complex a, f_complex b)
{
	f_complex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}

f_complex Cmul( f_complex a, f_complex b)
{
	f_complex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

f_complex Cdiv( f_complex a, f_complex b)
{
	f_complex c;
	float r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

f_complex Complex( float re, float im)
{
	f_complex c;
	c.r=re;
	c.i=im;
	return c;
}

float Cabs( f_complex z )
{
	float x,y,ans,temp;

	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp = y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

float Carg( f_complex z )
{
	float	x,y,ans;

	ans = atan2 ( z.i, z.r );

	return ans;
}

f_complex Conjg( f_complex z )
{
	f_complex c;
	c.r = z.r;
	c.i = -z.i;
	return c;
}

f_complex Csqrt( f_complex z )
{
	f_complex c;
	float x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=c.r=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y ) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

f_complex RCmul( float x, f_complex a)
{
	f_complex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}

f_complex Cpow( f_complex z, int n)
{
	f_complex c;
	if (n < 0)
		fprintf(stderr,"power must be a positive integer");
	for ( c = Complex(1.0,0.0); n > 0; --n)
		c = Cmul(c,z);
	return c;
}

float Re( f_complex z )
{
	return z.r;
}

float Im( f_complex z)
{
	return z.i;
}
