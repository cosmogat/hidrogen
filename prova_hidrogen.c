/* 28-02-2018 */
/* alex */
/* hidrogen.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "hidrogen.h"

int main () {
  /* real x, y, z, r, t, p, h = 0.5, radial, angular, prob; */
  /* int n, l, m, i, j, k, N = 1, Ni = floor(PI / h) , Nj = floor((2.0 * PI) / h), Nk = floor(20.0 / h); */
  /* char nom[20]; */
  real h_pol = 0.01, x;
  int i;
  /* real h = 0.05, t, p, x, y, z, harm_modQ, l, m; */
  /* int i, j, Ni = floor(PI / h) , Nj = floor((2.0 * PI) / h); */

  FILE * fit;

  /* for (n = 1; n <= N; n++) { */
  /*   for (l = 0; l < n; l++) { */
  /*     for (m = -l; m <= l; m++) { */
  /* 	sprintf(nom, "dad%d_%d_%d.dat", n, l, m); */
  /* 	fit = fopen(nom, "w"); */
  /* 	for (k = 0; k < Nk; k++) { */
  /* 	  r = ((real) k) * h; */
  /* 	  radial = part_radial(n, l, r); */
  /* 	  for (i = 0; i <= Ni; i++) { */
  /* 	    t = ((real) i) * h; */
  /* 	    angular = harmonic_esferic_mod(l, m, t); */
  /* 	    for (j = 0; j <= Nj; j++) { */
  /* 	      p = ((real) j) * h; */
  /* 	      prob = r * r * radial * radial * angular * angular; */
  /* 	      x = sin(t) * cos(p); */
  /* 	      y = sin(t) * sin(p); */
  /* 	      z = cos(t);	       */
  /* 	      fprintf(fit, "%f %f %f %f %f %f\n", x, y, z, x * prob, y * prob, z * prob); */
  /* 	    } */
  /* 	  } */
  /* 	} */
  /* 	fclose(fit); */
  /*     } */
  /*   } */
  /* } */
  
  fit = fopen("dades.dat", "w");
  /* for (i = 0; i <= 200; i++) { */
  /*   x = -1.0 + (i * h_pol); */
  /*   // fprintf(fit, "%f %f %f %f %f %f\n", x, pol_legendre_r(0, x), pol_legendre_r(1, x), pol_legendre_r(2, x), pol_legendre_r(3, x), pol_legendre_r(4, x)); */
  /*   // fprintf(fit, "%f %f %f %f %f %f %f\n", x, pol_legendre(0, x), pol_legendre(1, x), pol_legendre(2, x), pol_legendre(3, x), pol_legendre(4, x), pol_legendre(5, x)); */
  /*   // fprintf(fit, "%f %f %f %f %f %f\n", x, pol_legendre_ass(2, 2, x), pol_legendre_ass(3, 2, x), pol_legendre_ass(4, 2, x), pol_legendre_ass(5, 2, x), pol_legendre_ass(6, 2, x)); */
  /*   // fprintf(fit, "%f %f %f %f %f %f %f\n", x, pol_legendre_ass(3, -3, x), pol_legendre_ass(3, -2, x), pol_legendre_ass(3, -1, x), (15.0 / 720.0) * pow(1.0 - (x*x), 1.5), (15.0 / 120.0) * x * (1.0 - (x * x)), (3.0 / 24.0) * ((5.0 * x * x) - 1.0) * sqrt(1.0 - (x * x))); */
  /*   // fprintf(fit, "%f %f %f %f %f %f %f\n", x, pol_legendre_ass(3, 0, x), pol_legendre_ass(3, 1, x), pol_legendre_ass(3, 2, x), 0.5 * ((5.0 * x * x * x) - (3.0 * x)), -1.5 * ((5.0 * x * x) - 1.0) * sqrt(1.0 - (x * x)), 15.0 * x * (1.0 - (x * x))); */
  /*   // fprintf(fit, "%f %f %f %f %f %f %f\n", x, pol_laguerre(0, x), pol_laguerre(1, x), pol_laguerre(2, x), pol_laguerre(3, x), 0.5 * ((x * x) - (4.0 * x) + 2.0), (-(x * x * x) + (9.0 * x * x) - (18.0 * x) + 6.0) / 6.0); */
  /* } */
  /* for (i = 0; i <= 1200; i++) { */
  /*   x = -2.0 + (i * h_pol); */
  /*   fprintf(fit, "%f %f %f %f %f %f %f %f %f %f %f\n", x, pol_laguerre_ass(0, 0, x), pol_laguerre_ass(1, 0, x), pol_laguerre_ass(1, 1, x), pol_laguerre_ass(2, 0, x), pol_laguerre_ass(2, 1, x), pol_laguerre_ass(2, 2, x), pol_laguerre_ass(3, 0, x), pol_laguerre_ass(3, 1, x), pol_laguerre_ass(3, 2, x), pol_laguerre_ass(3, 3, x)); */
  /* } */
  for (i = 0; i <= 1200; i++) {
    x = 0.0 + (i * h_pol);
    fprintf(fit, "%f ", x);
    fprintf(fit, "%f ", x * x * part_radial(1, 0, x) * part_radial(1, 0, x));
    fprintf(fit, "%f ", x * x * part_radial(2, 0, x) * part_radial(2, 0, x));
    fprintf(fit, "%f ", x * x * part_radial(2, 1, x) * part_radial(2, 1, x));
    fprintf(fit, "%f ", x * x * part_radial(3, 0, x) * part_radial(3, 0, x));
    fprintf(fit, "%f ", x * x * part_radial(3, 1, x) * part_radial(3, 1, x));
    fprintf(fit, "%f ", x * x * part_radial(3, 2, x) * part_radial(3, 2, x));
    fprintf(fit, "\n");
  }
  /* l = 3; */
  /* m = -2; */
  /* for (i = 0; i <= Ni; i++) { */
  /*   t = ((real) i) * h; */
  /*   harm_modQ = harmonic_esferic_modQ(l, m, t); */
  /*   for (j = 0; j <= Nj; j++) { */
  /*     p = ((real) j) * h; */
  /*     x = sin(t) * cos(p); */
  /*     y = sin(t) * sin(p); */
  /*     z = cos(t); */
  /*     fprintf(fit, "%f %f %f %f %f %f\n", x, y, z, x * harm_modQ, y * harm_modQ, z * harm_modQ); */
  /*   } */
  /* } */
  fclose(fit);
  return 0;
}
