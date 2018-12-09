#ifndef CMH_MACROS_
#define CMH_MACROS_

/* macros.h -- wrapping macros around mpfr/mpc
 *
 * Copyright (C) 2006, 2010, 2011, 2012, 2013 INRIA
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

#include "mpc.h"
#if MPC_VERSION_MAJOR < 1
/* Compatibility fix */
#define mpc_mul_2ui mpc_mul_2exp
#define mpc_div_2ui mpc_div_2exp
#endif


/* {{{ macros */
#define ASSERT(x)	assert(x)

#define croak__(x,y) do {						\
        fprintf(stderr,"%s in %s at %s:%d -- %s\n",			\
                (x),__func__,__FILE__,__LINE__,(y));			\
    } while (0)

#define ASSERT_ALWAYS(x)						\
    do {								\
        if (!(x)) {							\
            croak__("code BUG() : condition " #x " failed",		\
                    "Abort");						\
            abort();							\
        }								\
    } while (0)
/* }}} */

#define max(a,b)		(((a)>(b)) ? (a) : (b))
#define min(a,b)		(((a)>(b)) ? (b) : (a))

/* To handle mpc_t */

#define cinit(a,n)  		mpc_init2((a), (n))
#define cclear(a)		mpc_clear((a))
#define cset_str(a,s)           mpc_set_str((a), (s), 10, MPC_RNDNN)
#define cadd(a,b,c) 		mpc_add((a), (b), (c), MPC_RNDNN)
#define cadd_ui(a,b,c)		mpc_add_ui((a), (b), (unsigned long int) (c), MPC_RNDNN)
#define cadd_i(a,b,c)           mpc_add_si((a), (b), (long int) (c), MPC_RNDNN)
#define cadd_fr(a,b,c)		mpc_add_fr ((a), (b), (c), MPC_RNDNN)
#define csub(a,b,c) 		mpc_sub((a), (b), (c), MPC_RNDNN)
#define csub_ui(a,b,c) 		mpc_sub_ui((a), (b), (unsigned long int) (c), MPC_RNDNN)
#define cmul(a,b,c) 		mpc_mul((a), (b), (c), MPC_RNDNN)
#define cmul_ui(a,b,c)		mpc_mul_ui((a), (b), (unsigned long int) (c), MPC_RNDNN)
#define cmul_i(a,b,c)           mpc_mul_si((a), (b), (long int) (c), MPC_RNDNN)
#define cmul_by_i(a,b)          mpc_mul_i((a), (b), 1, MPC_RNDNN)
#define cmul_2ui(a,b,c)	        mpc_mul_2ui((a), (b), (unsigned long int) (c), MPC_RNDNN)
#define cmul_fr(a,b,c)		mpc_mul_fr((a), (b), (c), MPC_RNDNN)
#define cdiv(a,b,c) 		mpc_div((a), (b), (c), MPC_RNDNN)
#define cdiv_ui(a,b,c)		mpc_div_ui((a), (b), (unsigned long int) (c), MPC_RNDNN)
#define cdiv_2ui(a,b,c)	        mpc_div_2ui((a), (b), (unsigned long int) (c), MPC_RNDNN)
#define csqr(a,b) 		mpc_sqr((a), (b), MPC_RNDNN)
#define csqrt(a,b) 		mpc_sqrt((a), (b), MPC_RNDNN)
#define cinv(a,b)   		mpc_ui_div((a), (unsigned long int) 1, (b), MPC_RNDNN)
#define cexp(a,b)		mpc_exp((a), (b), MPC_RNDNN)
#define clog(a,b)		mpc_log((a), (b), MPC_RNDNN)
#define cpow_d(a,b,c)		mpc_pow_d((a), (b), (c), MPC_RNDNN)
#define cpow_fr(a,b,c)		mpc_pow_fr((a), (b), (c), MPC_RNDNN)
#define cset(a,b)		mpc_set((a), (b), MPC_RNDNN)
#define cset_z(a,b)		mpc_set_z((a), (b), MPC_RNDNN)
#define cset_prec(a,b)		mpc_set_prec((a), (b))
#define cset_ui(a,b)		mpc_set_ui((a), (unsigned long int) (b), MPC_RNDNN)
#define cneg(a,b)		mpc_neg((a), (b), MPC_RNDNN)
#define czero(a)		mpc_set_ui((a), (unsigned long int) 0, MPC_RNDNN)
#define cone(a)			mpc_set_ui((a), (unsigned long int) 1, MPC_RNDNN)
#define cnorm(a,b)		mpc_norm((a), (b), MPC_RNDNN)
#define cabs(a,b)		mpc_abs((a), (b), GMP_RNDN)
#define cconj(a,b)		mpc_conj((a), (b), GMP_RNDN)
#define cprint1(a,n)    	mpc_out_str(stderr, 10, (n), (a), MPC_RNDNN)
#define cprint2(a,n)    	mpc_out_str(stderr, 2, (n), (a), MPC_RNDNN)
#define cswap(a,b)		mpc_swap((a), (b))
#define creal(a)                mpc_realref((a))
#define cimag(a)                mpc_imagref((a))

#define ccmp(a,b)		mpc_cmp((a), (b))
#define ccos(a,b)		mpc_cos((a), (b), MPC_RNDNN)
#define cpow(a,b,c)		mpc_pow((a), (b), (c), MPC_RNDNN)

#define MPC_RE(op)      mpc_realref(op)
#define MPC_IM(op)      mpc_imagref(op)

/* To handle mpfr_t */
#ifndef MPFR_PREC
#define MPFR_PREC(x)      ((x)->_mpfr_prec)
#define MPFR_EXP(x)       ((x)->_mpfr_exp)
#define MPFR_SIGN(x)       ((x)->_mpfr_sign)
#define MPFR_MANT(x)      ((x)->_mpfr_d)
#define MPFR_LIMB_SIZE(x) ((MPFR_PREC((x))-1)/GMP_NUMB_BITS+1)
#define MPFR_GET_ALLOC_SIZE(x) \
 ( ((mp_size_t*) MPFR_MANT(x))[-1] + 0)
#define MPFR_SET_ALLOC_SIZE(x, n) \
 ( ((mp_size_t*) MPFR_MANT(x))[-1] = n)
#endif

#define finit(a,n)  		mpfr_init2((a), (n))
#define fclear(a)		mpfr_clear((a))
#define fadd(a,b,c) 		mpfr_add((a), (b), (c), GMP_RNDN)
#define fadd_ui(a,b,c) 		mpfr_add_ui((a), (b), (unsigned long int) (c), GMP_RNDN)
#define fadd_z(a,b,c) 		mpfr_add_z((a), (b), (c), GMP_RNDN)
#define fsub(a,b,c) 		mpfr_sub((a), (b), (c), GMP_RNDN)
#define fsub_ui(a,b,c) 		mpfr_sub_ui((a), (b), (unsigned long int) (c), GMP_RNDN)
#define fsub_z(a,b,c) 		mpfr_sub_z((a), (b), (c), GMP_RNDN)
#define fmul(a,b,c) 		mpfr_mul((a), (b), (c), GMP_RNDN)
#define fmul_ui(a,b,c)		mpfr_mul_ui((a), (b), (unsigned long int) (c), GMP_RNDN)
#define fmul_2ui(a,b,c) 	mpfr_mul_2ui((a), (b), (unsigned long int) (c), GMP_RNDN)
#define fdiv(a,b,c) 		mpfr_div((a), (b), (c), GMP_RNDN)
#define fdiv_ui(a,b,c) 		mpfr_div_ui((a), (b), (unsigned long int) (c), GMP_RNDN)
#define fdiv_z(a,b,c)           mpfr_div_z((a), (b), (c), GMP_RNDN)
#define fui_div(a,b,c)		mpfr_ui_div((a), (unsigned long int) (b), (c), GMP_RNDN)
#define fdiv_2ui(a,b,c) 	mpfr_div_2ui((a), (b), (unsigned long int) (c), GMP_RNDN)
#define finv(a,b)   		mpfr_ui_div((a), (unsigned long int) 1, (b), GMP_RNDN)
#define fsqrt(a,b)		mpfr_sqrt((a), (b), GMP_RNDN)
#define fset(a,b)		mpfr_set((a), (b), GMP_RNDN)
#define fset_ui(a,b)		mpfr_set_ui((a), (unsigned long int) (b), GMP_RNDN)
#define fset_si(a,b)		mpfr_set_si((a), (long int) (b), GMP_RNDN)
#define fneg(a,b)		mpfr_neg((a), (b), GMP_RNDN)
#define fabs(a,b)		mpfr_abs((a), (b), GMP_RNDN)
#define fcmp(a,b)		mpfr_cmp((a), (b))
#define fcmp_ui(a,b)		mpfr_cmp_ui((a), (unsigned long int) (b))
#define fcmp_d(a,b)		mpfr_cmp_d((a), (double) (b))
#define fsgn(a)			mpfr_sgn((a))
#define frand(a)		mpfr_random((a))
#define fprint1(a,n)    	mpfr_out_str(stderr, 10, (n), (a), GMP_RNDN)
#define fprint2(a,n)    	mpfr_out_str(stderr, 2, (n), (a), GMP_RNDN)
#define fprec(a)		mpfr_get_prec((a))
#define cprec_round(a,n)     do {    \
    mpfr_prec_round(MPC_RE((a)), (n), GMP_RNDN);        \
    mpfr_prec_round(MPC_IM((a)), (n), GMP_RNDN);        \
} while (0)

#define fone(a)			mpfr_set_ui((a), (unsigned long int) 1, GMP_RNDN)
#define fzero(a)		mpfr_set_ui((a), (unsigned long int) 0, GMP_RNDN)
#define flog(a,b)		mpfr_log((a), (b), GMP_RNDN)
#define fget_ui(a)		mpfr_get_ui((a), GMP_RNDN)

/* Misc. functions */

#define fexp(a)                 mpfr_get_exp((a))
#define cprec(a)    		mpfr_get_prec(MPC_RE(a))
#define frint(a,b)		mpfr_rint((a), (b), GMP_RNDN);


/* mpfr/mpc output for debuggging */
#define MPC_OUT(x, prec)                                        \
do {                                                            \
  printf (#x "[%lu,%lu]=",                                      \
      (unsigned long int) mpfr_get_prec (mpc_realref (x)),      \
      (unsigned long int) mpfr_get_prec (mpc_imagref (x)));     \
  mpc_out_str (stdout, 10, prec, x, MPC_RNDNN);                 \
  printf ("\n");                                                \
} while (0)

#define MPFR_OUT(x, prec)                                       \
do {                                                            \
  printf (#x "[%lu]=", (unsigned long int) mpfr_get_prec (x));  \
  mpfr_out_str (stdout, 10, prec, x, GMP_RNDN);                 \
  printf ("\n");                                                \
} while (0)

#define MPFRX_OUT(x, prec)                                      \
do {                                                            \
  printf (#x "[%lu]=", (unsigned long int) mpfrx_get_prec (x)); \
  mpfrx_out_str (stdout, 10, prec, x);                          \
  printf ("\n");                                                \
} while (0)


#endif
