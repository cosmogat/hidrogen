/* 02-03-2018 */
/* alex */
/* hidrogen.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "hidrogen.h"


real pol_legendre_r(unsigned int n, real x) {
  real i_r = (real) n;
  if (n == 0)
    return 1.0;
  else if (n == 1)
    return x;
  else
    return ((((2 * i_r) - 1.0) / i_r) * x * pol_legendre_r(n - 1, x)) - (((i_r - 1) / i_r) * pol_legendre_r(n - 2, x));
}

real pol_legendre(unsigned int n, real x) {
  unsigned int i;
  real pol[n + 1], i_r;
  if (n == 0)
    return 1.0;
  pol[0] = 1.0;
  pol[1] = x;
  for (i = 2; i <= n; i++) {
    i_r = (real) i;
    pol[i] = ((((2.0 * i_r) - 1.0) * x * pol[i - 1]) - ((i_r - 1.0) * pol[i - 2])) / i_r;
  }
  return pol[n];
}

real pol_legendre_ass(unsigned int l, int m, real x) {
  unsigned int  m_abs = abs(m), i, it_l, it_m;
  int fac = 1;
  real pol_l[l - m_abs + 1], pol_m, l_r, m_r, fac_r, arrel;
  if ((x < -1.0) || (x > 1.0)) {
    fprintf(stderr, "Els polinomis de Legendre associats estan definits per a x entre -1.0 i 1.0\n");
    exit(-1);
  }
  if (m_abs > l)
    return 0.0;
  pol_m = 1.0;
  arrel = sqrt(1.0 - (x * x));
  /* {P_a}^a fins a = m */
  for (it_m = 1; it_m <= m_abs; it_m++) {
    m_r = (real) it_m;
    pol_m *= (1.0 - (2.0 * m_r)) * arrel;
  }
  /* {P_b}^m, fins b = l*/
  pol_l[0] = pol_m;
  if (m_abs != l) {
    m_r = (real) m_abs;
    pol_l[1] = ((2.0 * m_abs) + 1.0) * x * pol_l[0];
    for (it_l = (m_abs + 2); it_l <= l; it_l++) {
      l_r = (real) it_l;
      m_r = (real) m_abs;
      i = it_l - m_abs;
      pol_l[i] = ((((2.0 * l_r) - 1.0) / (l_r - m_r)) * x * pol_l[i - 1]) - (((l_r + m_r - 1.0) / (l_r - m_r)) * pol_l[i - 2]);
    }
  }
  /* Si m < 0, {P_l}^{-m} */
  if (m < 0) {
    for (i = (l - m_abs + 1); i <= (l + m_abs); i++)
      fac *= i;
    if ((m_abs % 2) == 0)
      fac_r = (real) fac;
    else
      fac_r = (real) (-fac);
    return pol_l[l - m_abs] / fac_r;
  }
  return pol_l[l - m_abs];
}

real pol_laguerre(int n, real x) {
  int i;
  real pol[n + 1], i_r;
  if (n == 0)
    return 1.0;
  if (n < 0)
    return exp(x) * pol_laguerre((-n) - 1, -x);
  pol[0] = 1.0;
  pol[1] = 1.0 - x;
  for (i = 2; i <= n; i++) {
    i_r = (real) i;
    pol[i] = ((((2.0 * i_r) - 1.0 - x) * pol[i - 1]) - ((i_r - 1.0) * pol[i - 2])) / i_r;
  }
  return pol[n];
}

real pol_laguerre_ass(unsigned int n, unsigned int m, real x) {
  unsigned int i;
  real pol[n + 1], i_r, m_r = (real) m;
  if (n == 0)
    return 1.0;
  if (m > n)
    return 0.0;
  pol[0] = 1.0;
  pol[1] = 1.0 + m_r - x;
  for (i = 2; i <= n; i++) {
    i_r = (real) i;
    pol[i] = ((((2.0 * i_r) - 1.0 + m_r - x) * pol[i - 1]) - ((i_r - 1.0 + m_r) * pol[i - 2])) / i_r;
  }
  return pol[n];
}

comp harmonic_esferic(unsigned int l, int m, real theta, real phi) {
  unsigned int i, fac = 1;
  real fac_r, l_r = (real) l, part_re, part_im;
  comp ret;
  for (i = (l - m + 1); i <= (l + m); i++)
    fac *= i;
  fac_r = (real) fac;
  fac_r = sqrt(((2.0 * l_r) + 1.0) / (fac_r * 4.0 * PI)) * pol_legendre_ass(l, m, cos(theta));
  part_re = fac_r * cos(phi);
  part_im = fac_r * sin(phi);
  ret = part_re + (part_im * I);
  return ret;
}

real harmonic_esferic_re(unsigned int l, int m, real theta, real phi) {
  unsigned int i, fac = 1;
  real fac_r, l_r = (real) l;
  for (i = (l - m + 1); i <= (l + m); i++)
    fac *= i;
  fac_r = (real) fac;
  fac_r = sqrt(((2.0 * l_r) + 1.0) / (fac_r * 4.0 * PI)) * pol_legendre_ass(l, m, cos(theta));
  return (fac_r * cos(phi));
}

real harmonic_esferic_im(unsigned int l, int m, real theta, real phi) {
  unsigned int i, fac = 1;
  real fac_r, l_r = (real) l;
  for (i = (l - m + 1); i <= (l + m); i++)
    fac *= i;
  fac_r = (real) fac;
  fac_r = sqrt(((2.0 * l_r) + 1.0) / (fac_r * 4.0 * PI)) * pol_legendre_ass(l, m, cos(theta));
  return (fac_r * sin(phi));
}

real harmonic_esferic_factor(unsigned int l, int m, real theta) {
  unsigned int i, fac = 1;
  real fac_r, l_r = (real) l;
  for (i = (l - m + 1); i <= (l + m); i++)
    fac *= i;
  fac_r = (real) fac;
  return sqrt(((2.0 * l_r) + 1.0) / (fac_r * 4.0 * PI)) * pol_legendre_ass(l, m, cos(theta));
}

real harmonic_esferic_mod(unsigned int l, int m, real theta) {
  real he_fac = harmonic_esferic_factor(l, m, theta);
  if (he_fac < 0)
    return -he_fac;
  else
    return he_fac;
}

real harmonic_esferic_modQ(unsigned int l, int m, real theta) {
  unsigned int i, fac = 1;
  real fac_r, l_r = (real) l, legen = pol_legendre_ass(l, m, cos(theta));
  for (i = (l - m + 1); i <= (l + m); i++)
    fac *= i;
  fac_r = (real) fac;
  return ((2.0 * l_r) + 1.0) / (fac_r * 4.0 * PI) * legen * legen;
}

real part_radial(unsigned int n, unsigned int l, real r) {
  unsigned int i, fac1 = 1, fac2 = 1;
  real aa0, rho, fac, n_r = (real) n;
  if (n == 0) {
    fprintf(stderr, "El nombre quàntic n ha de ser major que 1.\n");
    exit(-1);
  }
  else if (l >= n) {
    fprintf(stderr, "El nombre quàntic l ha de ser menor o igual que n - 1.\n");
    exit(-1);
  }
  aa0 = 0.5293211126080956;
  rho = (2.0 * r) / (n_r * aa0);
  for (i = 1; i <= (n - l - 1); i++)
    fac1 *= i;
  fac2 = fac1;
  for (i = n - l; i <= (n + l); i++)
    fac2 *= i;
  fac = ((real) fac1) / (2.0 * n_r * ((real) fac2));
  fac = sqrt(pow((2.0 / (n_r * aa0)), 3) * fac);
  return fac * exp(-0.5 * rho) * pow(rho, l) * pol_laguerre_ass(n - l - 1, (2 * l) + 1, rho);
}

real part_radialQ(unsigned int n, unsigned int l, real r) {
  unsigned int i, fac1 = 1, fac2 = 1;
  real aa0, rho, fac, n_r = (real) n, lague;
  if (n == 0) {
    fprintf(stderr, "El nombre quàntic n ha de ser major que 1.\n");
    exit(-1);
  }
  else if (l >= n) {
    fprintf(stderr, "El nombre quàntic l ha de ser menor o igual que n - 1.\n");
    exit(-1);
  }
  aa0 = 0.5293211126080956;
  rho = (2.0 * r) / (n_r * aa0);
  for (i = 1; i <= (n - l - 1); i++)
    fac1 *= i;
  fac2 = fac1;
  for (i = n - l; i <= (n + l); i++)
    fac2 *= i;
  fac = ((real) fac1) / (2.0 * n_r * ((real) fac2));
  fac *= pow((2.0 / (n_r * aa0)), 3);
  lague =  pol_laguerre_ass(n - l - 1, (2 * l) + 1, rho);
  return fac * exp(-rho) * pow(rho, 2 * l) * lague * lague;
}

comp psi_hidrogen(unsigned int n, unsigned int l, int m, real r, real theta, real phi) {
  real radial = part_radial(n, l, r);
  real angular = harmonic_esferic_factor(l, m, theta);
  real fac = radial * angular;
  real part_re, part_im;
  comp ret;
  part_re = fac * cos(phi);
  part_im = fac * sin(phi);
  ret = part_re + (part_im * I);
  return ret;
}
