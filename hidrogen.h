/* 02-03-2018 */
/* alex */
/* hidrogen.h */
#ifndef _HIDROGEN_H
#define _HIDROGEN_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#undef PI
#define PI 3.14159265358979324

typedef double real;
typedef double complex comp;
real pol_legendre_r(unsigned int n, real x);
real pol_legendre(unsigned int n, real x);
real pol_legendre_ass(unsigned int l, int m, real x);
real pol_laguerre(int n, real x);
real pol_laguerre_ass(unsigned int n, unsigned int m, real x);
comp harmonic_esferic(unsigned int l, int m, real theta, real phi);
real harmonic_esferic_re(unsigned int l, int m, real theta, real phi);
real harmonic_esferic_im(unsigned int l, int m, real theta, real phi);
real harmonic_esferic_factor(unsigned int l, int m, real theta);
real harmonic_esferic_mod(unsigned int l, int m, real theta);
real harmonic_esferic_modQ(unsigned int l, int m, real theta);
real part_radial(unsigned int n, unsigned int l, real r);
real part_radialQ(unsigned int n, unsigned int l, real r);
comp psi_hidrogen(unsigned int n, unsigned int l, int m, real r, real theta, real phi);
#endif
