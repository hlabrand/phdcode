#ifndef MISC_H_
#define MISC_H_

/* misc.h -- headers for misc.c
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

#ifdef __cplusplus
extern "C" {
#endif

  
// Compute ceil(log_2(i))
unsigned long int ceilLog2OfInteger(unsigned long int i);
  
  
/* TODO: It seems that these two functions subtly differ, in that the
 * second one wants the most accomodating divisor (largest module), and
 * the first one wants the property to hold for both divisors.
 */

/* n such that |a-b|/|a| <= 1/2^n and |a-b|/|b| <= 1/2^n */
int creldist (mpc_t a, mpc_t b);

/* This takes two complex numbers a and b. It returns a precision p such
 * that |(a-b)/X| <= 2^-p, with X the complex number of largest module
 * between a and b.
 */
//mpfr_prec_t agreeing_bits_n (mpc_t * f, mpc_t * g, int n);


/* which of a or -a is closest to b ? */
int does_negation_bring_closer(mpc_t a, mpc_t b);

/* Find (y[0]:y[1]:y[2]:y[3]) projectively close to (x[0]:x[1]:x[2]:x[3])
 *
 * The x[] are not given. Only the sign of their real and imaginary part
 * is.  We assume that none of the y[i] is zero, and that there exists a
 * vector of signs s[0..3] such that (y[i]::) = (s[i]*x[i]::).
 *
 * y[i] is chosen so that the largest (in absolute value) of its real and
 * imaginary part has sign matching the one given for x[i].
 *
 * We further assert that the real part of y[0] is positive in the end.
 * This is, I believe, essential for the rest of the computation.
 *
 * A integer is returned, where bit i is 0 if s[i] = 1, and 1 if s[i] = -1.
 */
int find_projectively_close_quadruple(mpc_t * y, int xsign[4][2]);

//int need_negation_to_reach_sign(mpc_t y, int xsign[2]);

#ifdef __cplusplus
}
#endif

#endif	/* MISC_H_ */
