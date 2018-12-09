/* misc.c - miscellaneous helper functions.
 *
 * Copyright (C) 2013 INRIA
 *
 * This file is part of CMH.
 *
 * CMH is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the License, or (at your
 * option) any later version.
 *
 * CMH is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/ .
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpc.h>

#include "macros.h"
#include "misc.h"


/* computes the integer ceil(log_2(int)) */
 
unsigned long int ceilLog2OfInteger(unsigned long int i)
{
  int bits = 0;
  while (i > 255) {
    bits = bits + 8;
    i >>= 8;
  }
  while (i > 0) {
    bits++;
    i >>= 1;
  }
  return bits;
}







/*
 * The following function returns an integer n such that:
 *
 *  |a-b|/|a| <= 1/2^n   and    |a-b|/|b| <= 1/2^n
 *
 */
int
creldist (mpc_t a, mpc_t b)
{
   int A, B, D;
   mpc_t d;

   cinit(d, cprec(a));
   csub(d, a, b);
   if ( ! fsgn(MPC_RE(a)) )
      if ( ! fsgn(MPC_IM(a)) )
         A = INT_MAX;
      else
         A = fexp(MPC_IM(a));
   else
      if ( ! fsgn(MPC_IM(a)) )
         A = fexp(MPC_RE(a));
      else
         A = max(fexp(MPC_RE(a)), fexp(MPC_IM(a)));
    
   if ( ! fsgn(MPC_RE(b)) )
      if ( ! fsgn(MPC_IM(b)) )
         B = INT_MAX;
      else
         B = fexp(MPC_IM(b));
   else
      if ( ! fsgn(MPC_IM(b)) )
         B = fexp(MPC_RE(b));
      else
         B = max(fexp(MPC_RE(b)), fexp(MPC_IM(b)));
   A = min(A, B);


   if ( ! fsgn(MPC_RE(d)) )
      if ( ! fsgn(MPC_IM(d)) ) {
         cclear(d);
         return INT_MAX;
      }
      else {
         D = fexp(MPC_IM(d));
         cclear(d);
         return A-D+2;
      }
   else
      if ( ! fsgn(MPC_IM(d)) ) {
         D = fexp(MPC_RE(d));
         cclear(d);
         return A-D+2;
      }
      else {
         D = max(fexp(MPC_RE(d)), fexp(MPC_IM(d)));
         cclear(d);
         return A-D+2;
      }
}


/* This takes two complex numbers a and b. It returns a precision p such
 * that |(a-b)/X| <= 2^-p, with X the complex number of largest module
 * between a and b.
 */
mpfr_prec_t agreeing_bits (mpc_srcptr a, mpc_srcptr b)
{
    if (!mpfr_number_p(mpc_realref(a))) return 0;
    if (!mpfr_number_p(mpc_imagref(a))) return 0;
    if (!mpfr_number_p(mpc_realref(b))) return 0;
    if (!mpfr_number_p(mpc_imagref(b))) return 0;

    mpc_t d;
    mpc_init2(d, mpc_get_prec(a));
    mpc_sub(d, a, b, MPC_RNDNN);
    mpfr_exp_t e, f;

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif
    e = MAX(mpfr_get_exp(mpc_realref(d)), mpfr_get_exp(mpc_imagref(d)));
    f = MAX(
            MAX(mpfr_get_exp(mpc_realref(a)), mpfr_get_exp(mpc_imagref(a))),
            MAX(mpfr_get_exp(mpc_realref(b)), mpfr_get_exp(mpc_imagref(b))));
    mpc_clear(d);

    /* Write a = aM * 2^aE, with aM a complex mantissa. Ditto for b, and d.
     *
     * Suppose a is normalized (either its real or imaginary part). Thus
     * we have \sqrt(2) >= |aM| >= 1/2. Now we suppose
     * |a-b|<=sqrt(2)*2^e. Within a and b, the one of largest module has
     * at least |X| >= 1/2 * 2^f, where f is the exponent. Therefore we
     * have |(a-b)/X|<=2\sqrt(2)*2^(e-f)<=2^-(f-e-2)
     */
    if (f <= e + 2)
        return 0;
    return f-e-2;
}

/*mpfr_prec_t agreeing_bits_n (mpc_t * f, mpc_t * g, int n)
{
    mpfr_prec_t ok = MPFR_PREC_MAX;
    for(int i = 0 ; ok > 0 && i < n ; i++) {
        mpfr_prec_t x = agreeing_bits(f[i], g[i]);
        if (x < ok) ok = x;
    }
    return ok;
}
*/
/* which of a or -a is closest to b ? We need to be more stable than
 * when just taking signs, because a and b may be imaginary, for
 * instance. */
int does_negation_bring_closer(mpc_t a, mpc_t b)
{
    mp_prec_t prec = cprec(a);
    mpfr_t zzf;
    mp_exp_t e[2][2];
    finit(zzf, prec);
    fsub(zzf, MPC_RE(b), MPC_RE(a)); e[0][0] = mpfr_get_exp(zzf);
    fadd(zzf, MPC_RE(b), MPC_RE(a)); e[1][0] = mpfr_get_exp(zzf);
    fsub(zzf, MPC_IM(b), MPC_IM(a)); e[0][1] = mpfr_get_exp(zzf);
    fadd(zzf, MPC_IM(b), MPC_IM(a)); e[1][1] = mpfr_get_exp(zzf);
    fclear(zzf);
    if (e[0][1] > e[0][0]) e[0][0] = e[0][1];
    if (e[1][1] > e[1][0]) e[1][0] = e[1][1];
    return e[0][0] > e[1][0];
}


#if 0
int find_projectively_close_quadruple(mpc_t * ra, mpc_t * low)
{
    int negate = 0;
    int loprec = cprec(low[0]);
    mpc_t aux[2];
    cinit (aux[0], loprec);
    cinit (aux[1], loprec);
    for(int i = 1 ; i < 4 ; i++) {
        /* Compare ra[i]/ra[0] and low[i]/low[0] */
        cmul(aux[0], ra[i], low[0]);
        cmul(aux[1], ra[0], low[i]);
        if (does_negation_bring_closer(aux[0], aux[1])) {
            negate |= 1<<i;
        }
    }
    cclear(aux[0]);
    cclear(aux[1]);
    return negate;
}
#endif

int need_negation_to_reach_sign(mpc_t y, int xsign[2])
{
    int sx, sy;
    if (fexp(creal(y)) > fexp(cimag(y))) {
        sy = fsgn(creal(y));
        sx = xsign[0];
    } else {
        sy = fsgn(cimag(y));
        sx = xsign[1];
    }
    /* Our decision strategy assumes that y is bounded away from
     * zero. In fact, calling this function with something which _is_
     * close to zero is not necessarily a problem, though, since the
     * indetermination in the sign choice is probably much less
     * likely to be a problem. So the check below is just a safety
     * guards because this escapes our primary intended use, not
     * more.
     */
    if (!sx || !sy) {
        fprintf(stderr, "Unexpected zero sign in need_negation_to_reach_sign. Please investigate\n");
        abort();
    }
    return (sx * sy == -1);
}

/*int find_projectively_close_quadruple(mpc_t * y, int xsign[4][2])
{
    int negate = 0;
    for(int i = 0 ; i < 4 ; i++) {
        negate |= need_negation_to_reach_sign(y[i], xsign[i]) << i;
    }
*/
    /* If this condition is encountered, though, then we are very
     * probably in trouble. Much of the rest of the code makes the
     * implicit assumption that the ``right'' choice has real part > 0,
     * at least for the first coordinate in a quadruple.
     */
/*    assert(!(negate & 1));
    return negate;
}
*/
