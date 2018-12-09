/* naivethetas.c -- computation of theta{00,01,10}({0,z},tau) by summing the series
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>

#include "macros.h"
#include "misc.h"
#include "naivethetas.h"







/*
 * Computes Theta(z, tau) at the given precision
 */

void NaiveTheta(mpc_t res, mpc_t z, mpc_t tau, unsigned long int p)
{
  mpc_t myz, mytau, q, w, pi, tmp, terms[4], theta, v1;
  mpfr_t bound, norme, ftmp, indice, un;
  // inits
  cinit(res, p);
  int myprec = 2*p;	// in pratice p + log B + 7 is enough
  cinit(theta, myprec);
  cinit(myz, myprec); cset(myz, z);
  cinit(mytau, myprec); cset(mytau, tau);
  cinit(tmp, myprec); czero(tmp);
  cinit(q, myprec);
  cinit(w, myprec);
  cinit(pi,myprec);
  cinit(v1, myprec);
  int i; for (i=0; i<4; i++) cinit(terms[i], myprec);

  finit(bound, myprec);
  finit(norme, myprec);
  finit(ftmp, myprec);
  finit(indice, myprec); fzero(indice);
  finit(un, myprec); fone(un);

  // pi
  cone(pi); cneg(pi, pi); clog(pi, pi); cmul_by_i(pi, pi);
  //fprintf(stderr, "pi :\n"); mpc_out_str(stderr, 10, 0, pi, MPC_RNDNN); fprintf(stderr, "\n");
  // w = e^(2 i z pi)
  cmul(w, z, pi); cmul_2ui(w, w,1); cmul_by_i(w,w); cexp(w, w);
  //fprintf(stderr, "w :\n"); mpc_out_str(stderr, 10, 0, w, MPC_RNDNN); fprintf(stderr, "\n");
  // q = e^(i pi tau)
  cmul(q, tau, pi); cmul_by_i(q,q); cexp(q,q);
  //fprintf(stderr, "q :\n"); mpc_out_str(stderr, 10, 0, q, MPC_RNDNN); fprintf(stderr, "\n");



  // If tau is fairly big, we can have theta = 1.00000000000000000000thing, and our approximation will be 1 ;
  //         then when we go to the finitedifferences, we find theta10(0)=0, and div by 0 = bye bye
  // So here, if we detect that the first term is really small, we call the function with a bigger precision so that there's at least 10 digits after the 1
  // CalculThetas won't mind, it starts with prec = whatever was sent back
  cset(tmp, w); cinv(theta, tmp); cadd(tmp, tmp, theta); cmul(tmp, tmp, q); cnorm(norme, tmp); fsqrt(norme, norme);
  fone(bound); fdiv_2ui(bound, bound, p);
  if ( fcmp(norme, bound) < 0) {
    // the second term is too small to be seen at that precision
    // we need to increase the precision so that there's at least a few digits at the end of the result
    mpfr_t log; finit(log, myprec); fset_ui(log, 2); flog(log, log); flog(norme, norme); fdiv(log, norme, log); fneg(log, log);
    //fprintf(stderr, "That tau is big!\n");
    NaiveTheta(res, z, tau, fget_ui(log) + p +1);
    return ;
  }

  // Determine the number of terms we need
  //     we use our Dupont-like bound
  fset(bound, creal(pi)); fset(ftmp, cimag(tau)); fmul(bound, bound, ftmp);      finv(bound, bound);
  // log_2(e) = log(e)/log(2) = 1/log2
  fset_ui(ftmp, 2); flog(ftmp, ftmp); fmul(bound, bound, ftmp); fmul_ui(bound, bound, p+3);
  fadd(bound, bound, un);
  fsqrt(bound, bound);
  fadd(bound, bound, un);


  // Summing
  cone(theta);
  for (i=0; i<2; i++) { cone(terms[i]); }
  int flag = 0;
  while (fcmp(bound, indice) >= 0)          // >=0 plays the role of ceil(bound)
  { // terms[0] = e^{i pi tau n}, terms[1] = e^{i pi tau n^2}, terms[2] = (e^{2 i n pi z} + e^{-2in pi z}), terms[3] = (e^{2 i (n-1) pi z} + e^{-2i (n-1) pi z})
    // n = n+1
    fadd(indice, indice, un);
    //fprintf(stderr, "indice :\n"); mpfr_out_str(stderr, 10, 0, indice, MPC_RNDNN); fprintf(stderr, "\n");

    // Compute the rank n+1

    // terms[1] = terms[1] * q^(2n+1)
    csqr(tmp, terms[0]); cmul(terms[1], terms[1], tmp); cmul(terms[1], terms[1], q);
    cmul(terms[0], terms[0], q);

    // v0 = 2, v1 = e^(2 i pi z)+e^(-2i pi z)
    if (flag == 0) {
      flag = 1;
      cinv(tmp, w); cadd(v1, w, tmp);
      cset(terms[2], v1); cset_ui(terms[3], 2);
      cset(tmp, v1); cmul(tmp, tmp, terms[1]); cadd(theta, theta, tmp);
      continue;
    }
    // (terms[2], terms[3]) = (terms[2]*v1-terms[3], terms[2])
    cmul(tmp, terms[2], v1); csub(tmp, tmp, terms[3]);
    cset(terms[3], terms[2]);
    cset(terms[2], tmp);

    // Add the term to the result
    cmul(tmp, tmp, terms[1]);
    cadd(theta, theta, tmp);
  }
  cset(res, theta);

  // clears
  cclear(myz); cclear(mytau);
  cclear(q); cclear(w);
  cclear(pi);
  cclear(tmp);
  cclear(theta);
  cclear(v1);
  for (i=0; i<4; i++) cclear(terms[i]);
  fclear(bound);
  fclear(norme);
  fclear(ftmp);
  fclear(indice);
  fclear(un);

  return;
}









/*
 * Returns (in that order) Theta00(z,tau), Theta00(0,tau), Theta01(z, tau), Theta01(0,tau) with precision p
 * Has some code to avoid returning theta00 = theta01 (which gives a division by 0 in our code)
 */


void InternalNaiveThetaFunctionAndThetaConstants_00_01(mpc_t res[4], mpc_t z, mpc_t tau, unsigned long int p)
{
  mpc_t myz, mytau, q, w, pi, tmp, terms[4], theta[4], v1;
  mpfr_t bound, norme, ftmp, indice, un;
  int sign;
  // inits
  int i; for (i=0; i<4; i++) { cinit(res[i], p); }
  int myprec = 2*p;	// in pratice p + log B + 7 is enough
  cinit(myz, myprec); cset(myz, z);
  cinit(mytau, myprec); cset(mytau, tau);
  cinit(tmp, myprec); czero(tmp);
  cinit(q, myprec);
  cinit(w, myprec);
  cinit(pi,myprec);
  cinit(v1,myprec);
  for (i=0; i<4; i++) cinit(terms[i], myprec);
  for (i=0; i<4; i++) cinit(theta[i], myprec);

  finit(bound, myprec);
  finit(norme, myprec);
  finit(ftmp, myprec);
  finit(indice, myprec); fzero(indice);
  finit(un, myprec); fone(un);

  // pi
  cone(pi); cneg(pi, pi); clog(pi, pi); cmul_by_i(pi, pi);
  //fprintf(stderr, "pi :\n"); mpc_out_str(stderr, 10, 0, pi, MPC_RNDNN); fprintf(stderr, "\n");
  // w = e^(2 i z pi)
  cmul(w, z, pi); cmul_2ui(w, w,1); cmul_by_i(w,w); cexp(w, w);
  //fprintf(stderr, "w :\n"); mpc_out_str(stderr, 10, 0, w, MPC_RNDNN); fprintf(stderr, "\n");
  // q = e^(i pi tau)
  cmul(q, tau, pi); cmul_by_i(q,q); cexp(q,q);
  //fprintf(stderr, "q :\n"); mpc_out_str(stderr, 10, 0, q, MPC_RNDNN); fprintf(stderr, "\n");


  // If tau is fairly big, we can have theta = 1.00000000000000000000thing, and our approximation will be 1 ;
  //         then when we go to the finitedifferences, we find theta10(0)=0, and div by 0 = bye bye
  // So here, if we detect that the first term is really small, we call the function with a bigger precision so that there's at least 10 digits after the 1
  // CalculThetas won't mind, it starts with prec = whatever was sent back
  cset(tmp, w); cinv(theta[0], tmp); cadd(tmp, tmp, theta[0]); cmul(tmp, tmp, q); cnorm(norme, tmp); fsqrt(norme, norme);
  fone(bound); fdiv_2ui(bound, bound, p);
  if ( fcmp(norme, bound) < 0) {
    // the second term is too small to be seen at that precision
    // we need to increase the precision so that there's at least a few digits at the end of the result
    mpfr_t log; finit(log, myprec); fset_ui(log, 2); flog(log, log); flog(norme, norme); fdiv(log, norme, log); fneg(log, log);
    //fprintf(stderr, "That tau is big!\n");
    NaiveThetaFunctionAndThetaConstants_00_01(res, z, tau, fget_ui(log) + p +1);
    return ;
  }

  // Determine the number of terms we need
  //     we use our Dupont-like bound
  fset(bound, creal(pi)); fset(ftmp, cimag(tau)); fmul(bound, bound, ftmp);      finv(bound, bound);
  // log_2(e) = log(e)/log(2) = 1/log2
  fset_ui(ftmp, 2); flog(ftmp, ftmp); fmul(bound, bound, ftmp); fmul_ui(bound, bound, p+3);
  fadd(bound, bound, un);
  fsqrt(bound, bound);
  fadd(bound, bound, un);

  // Summing
  for (i=0; i<4; i++) { cone(theta[i]); cone(terms[i]); }
  sign = 1;
  int flag = 0;
  while (fcmp(bound, indice) >= 0)          // >=0 plays the role of ceil(bound)
  { // terms[0] = e^{i pi tau n}, terms[1] = e^{i pi tau n^2}, terms[2] = (e^{2 i n pi z} + e^{-2in pi z}), terms[3] = (e^{2 i (n-1) pi z} + e^{-2i (n-1) pi z})
    // n = n+1
    fadd(indice, indice, un);
    //fprintf(stderr, "indice :\n"); mpfr_out_str(stderr, 10, 0, indice, MPC_RNDNN); fprintf(stderr, "\n");

    // Compute the rank n+1
    sign = -sign;

    csqr(tmp, terms[0]); cmul(terms[1], terms[1], tmp); cmul(terms[1], terms[1], q);
    cmul(terms[0], terms[0], q);

    // v0 = 2, v1 = e^(2 i pi z)+e^(-2i pi z)
    if (flag == 0) {
      flag = 1;
      cinv(tmp, w); cadd(v1, w, tmp);
      cset(terms[2], v1); cset_ui(terms[3], 2);
      cset(tmp, v1); cmul(tmp, tmp, terms[1]);
      cadd(theta[0], theta[0], tmp); cadd(theta[1], theta[1], terms[1]); cadd(theta[1], theta[1], terms[1]);
      csub(theta[2], theta[2], tmp); csub(theta[3], theta[3], terms[1]); csub(theta[3], theta[3], terms[1]);
      continue;
    }
    // (terms[2], terms[3]) = (terms[2]*v1-terms[3], terms[2])
    cmul(tmp, terms[2], v1); csub(tmp, tmp, terms[3]);
    cset(terms[3], terms[2]);
    cset(terms[2], tmp);

    // Add the terms to the result
    cmul(tmp, tmp, terms[1]);
    cadd(theta[0], theta[0], tmp);   cadd(theta[1], theta[1], terms[1]); cadd(theta[1], theta[1], terms[1]);
    if (sign == 1) {
      cadd(theta[2], theta[2], tmp); cadd(theta[3], theta[3], terms[1]); cadd(theta[3], theta[3], terms[1]);
    } else {
      csub(theta[2], theta[2], tmp); csub(theta[3], theta[3], terms[1]); csub(theta[3], theta[3], terms[1]);
    }
  }
  for (i=0; i<4; i++) { cset(res[i], theta[i]); }

  // clears
  cclear(myz); cclear(mytau);
  cclear(q); cclear(w);
  cclear(pi);
  cclear(tmp);
  cclear(v1);
  for (i=0; i<4; i++) { cclear(terms[i]); cclear(theta[i]); }
  fclear(bound);
  fclear(norme);
  fclear(ftmp);
  fclear(indice);
  fclear(un);

  return;
}





/*
 * Returns (in that order) Theta00(z,tau), Theta00(0,tau), Theta01(z, tau), Theta01(0,tau) with precision p
 */


void NaiveThetaFunctionAndThetaConstants_00_01(mpc_t res[4], mpc_t z, mpc_t tau, unsigned long int p)
{
  mpc_t myz, mytau, q, w, pi, tmp, tmp2, terms[4], theta[6], v1, saveq;
  mpfr_t bound, norme, ftmp, indice, un;
  int sign;
  int i;
  // inits
  unsigned long int resultprec = cprec(z);
  for (i=0; i<4; i++) { cinit(res[i], resultprec+2); }			// 2 because our result will be bounded by 4
  
  
  // p should be bigger than prec(z) to account for lost precision
  // the 2 is there because our quantities are bounded by 4, so representing them requires p+2 bits
  unsigned long int myprec;
  if (p <= resultprec) { myprec = resultprec + (ceilLog2OfInteger(p) + 10) + 2; }
  else { myprec = p + 2; }
  
  cinit(myz, myprec); cset(myz, z);
  cinit(mytau, myprec); cset(mytau, tau);
  cinit(tmp, myprec); czero(tmp);
  cinit(tmp2, myprec); czero(tmp);
  cinit(q, myprec);
  cinit(w, myprec);
  cinit(pi,myprec);
  cinit(v1,myprec);
  cinit(saveq, myprec);
  for (i=0; i<4; i++) cinit(terms[i], myprec);
  for (i=0; i<4; i++) cinit(theta[i], myprec);

  finit(bound, myprec);
  finit(norme, myprec);
  finit(ftmp, myprec);
  finit(indice, myprec); fzero(indice);
  finit(un, myprec); fone(un);

  // pi
  cone(pi); cneg(pi, pi); clog(pi, pi); cmul_by_i(pi, pi);
  //fprintf(stderr, "pi :\n"); mpc_out_str(stderr, 10, 0, pi, MPC_RNDNN); fprintf(stderr, "\n");
  // w = e^(2 i z pi)
  cmul(w, z, pi); cmul_2ui(w, w,1); cmul_by_i(w,w); cexp(w, w);
  //fprintf(stderr, "w :\n"); mpc_out_str(stderr, 10, 0, w, MPC_RNDNN); fprintf(stderr, "\n");
  // q = e^(i pi tau)
  cmul(q, tau, pi); cmul_by_i(q,q); cexp(q,q);
  //fprintf(stderr, "q :\n"); mpc_out_str(stderr, 10, 0, q, MPC_RNDNN); fprintf(stderr, "\n");


  // Determine the number of terms we need
  //     we use our Dupont-like bound
  fset(bound, creal(pi)); fset(ftmp, cimag(tau)); fmul(bound, bound, ftmp);      finv(bound, bound);
  // log_2(e) = log(e)/log(2) = 1/log2
  fset_ui(ftmp, 2); flog(ftmp, ftmp); fmul(bound, bound, ftmp); fmul_ui(bound, bound, p+3);
  fadd(bound, bound, un);
  fsqrt(bound, bound);
  fadd(bound, bound, un);

  // Summing
  for (i=0; i<4; i++) { cone(theta[i]); cone(terms[i]); }
  sign = 1;
  int flag = 0;
  while (fcmp(bound, indice) >= 0)          // >=0 plays the role of ceil(bound)
  { // terms[0] = e^{i pi tau n}, terms[1] = e^{i pi tau n^2}
    // terms[2] = q^{n^2}(w^{2n}+w^{-2n}), terms[3] = q^{(n-1)^2} (w^{2n-2}+w^{-2n+2})

    // n = n+1
    fadd(indice, indice, un);
    //fprintf(stderr, "indice :\n"); mpfr_out_str(stderr, 10, 0, indice, MPC_RNDNN); fprintf(stderr, "\n");

    // Compute the rank n+1
    sign = -sign;
    
    cset(saveq, terms[0]);
    csqr(tmp, terms[0]); cmul(terms[1], terms[1], tmp); cmul(terms[1], terms[1], q);
    cmul(terms[0], terms[0], q);


    // v0 = 2, v1 = q(w^2+w^{-2})
    if (flag == 0) {
      flag = 1;

      cmul(tmp, z, pi); cmul_2ui(tmp, tmp,1); cmul_by_i(tmp,tmp);  // 2 i pi z
      cexp(v1, tmp); cneg(tmp, tmp); cexp(tmp, tmp); cadd(v1, v1, tmp); // (w^2+w^{-2})
      cmul(v1, v1, q);
      cset(terms[2], v1); cset_ui(terms[3], 2);
      
      cadd(theta[0], theta[0], v1); cadd(theta[1], theta[1], terms[1]); cadd(theta[1], theta[1], terms[1]);
      csub(theta[2], theta[2], v1); csub(theta[3], theta[3], terms[1]); csub(theta[3], theta[3], terms[1]);

      continue;
    }

    // Lucas sequences
    //     (terms[2], terms[3]) = ( q^(2n)*v1*terms[2]-q^(4n) terms[3], terms[2])
    csqr(tmp2, saveq); csqr(tmp, tmp2); // q^2n, q^4n
    cmul(tmp2, tmp2, v1); cmul(tmp2, tmp2, terms[2]);
    cmul(tmp, tmp, terms[3]);
    csub(tmp, tmp2, tmp);
    cset(terms[3], terms[2]);
    cset(terms[2], tmp);

    // Add the terms to the result
    cadd(theta[0], theta[0], tmp);   cadd(theta[1], theta[1], terms[1]); cadd(theta[1], theta[1], terms[1]);
    if (sign == 1) {
      cadd(theta[2], theta[2], tmp); cadd(theta[3], theta[3], terms[1]); cadd(theta[3], theta[3], terms[1]);
    } else {
      csub(theta[2], theta[2], tmp); csub(theta[3], theta[3], terms[1]); csub(theta[3], theta[3], terms[1]);
    }
  }
  
  for (i=0; i<4; i++) { cset(res[i], theta[i]); }
  
  

  // clears
  cclear(myz); cclear(mytau);
  cclear(q); cclear(w);
  cclear(pi); cclear(v1);
  cclear(tmp); cclear(tmp2);
  for (i=0; i<4; i++) cclear(terms[i]);
  for (i=0; i<4; i++) cclear(theta[i]);
  fclear(bound);
  fclear(norme);
  fclear(ftmp);
  fclear(indice);
  fclear(un);


  return;
}









/*
 * Needs : z with precision P, tau with precision P, p = P + estimate number of lost bits
 *      if you don't give anything in p, we'll compute a bound of p + log p
 * Returns (in that order) Theta00(z,tau), Theta00(0,tau), Theta01(z, tau), Theta01(0,tau), Theta10(z,tau), Theta10(0,tau) with precision P
 * CAREFUL : do not EVER EVER use this to get a first approximation of the quotient:
 *              in \mathfrak{F}, you recover theta10(z,tau) by dividing an equation by theta10(0,tau): the latter is 0 if the quotient is 1
 *              the function above (_00_01) has code to make it so that the first approximation doesn't give a quotient of exactly 1, to avoid this problem
 *              this function below does NOT have this code (since it changes the asymptotic complexity); if you use it wrong, you'll get divisions by 0!
 */


void NaiveThetaFunctionAndThetaConstants_00_01_10(mpc_t res[6], mpc_t z, mpc_t tau, unsigned long int p)
{
  mpc_t myz, mytau, q, w, pi, tmp, terms[4], theta[6], v1, extraterms[4], extrav1[2];
  mpfr_t bound, norme, ftmp, indice, un;
  int sign;
  int i;
  // inits
  unsigned long int resultprec = cprec(z);
  for (i=0; i<6; i++) { cinit(res[i], resultprec+2); }			// 2 because our result will be bounded by 4
  
  
  // p should be bigger than prec(z) to account for lost precision
  // the 2 is there because our quantities are bounded by 4, so representing them requires p+2 bits
  unsigned long int myprec;
  if (p <= resultprec) { myprec = resultprec + (ceilLog2OfInteger(p) + 10) + 2; }
  else { myprec = p + 2; }
  
  cinit(myz, myprec); cset(myz, z);
  cinit(mytau, myprec); cset(mytau, tau);
  cinit(tmp, myprec); czero(tmp);
  cinit(q, myprec);
  cinit(w, myprec);
  cinit(pi,myprec);
  cinit(v1,myprec);
  for (i=0; i<4; i++) cinit(terms[i], myprec);
  for (i=0; i<4; i++) cinit(extraterms[i], myprec);
  for (i=0; i<2; i++) cinit(extrav1[i], myprec);
  for (i=0; i<6; i++) cinit(theta[i], myprec);

  finit(bound, myprec);
  finit(norme, myprec);
  finit(ftmp, myprec);
  finit(indice, myprec); fzero(indice);
  finit(un, myprec); fone(un);

  // pi
  cone(pi); cneg(pi, pi); clog(pi, pi); cmul_by_i(pi, pi);
  //fprintf(stderr, "pi :\n"); mpc_out_str(stderr, 10, 0, pi, MPC_RNDNN); fprintf(stderr, "\n");
  // w = e^(2 i z pi)
  cmul(w, z, pi); cmul_2ui(w, w,1); cmul_by_i(w,w); cexp(w, w);
  //fprintf(stderr, "w :\n"); mpc_out_str(stderr, 10, 0, w, MPC_RNDNN); fprintf(stderr, "\n");
  // q = e^(i pi tau)
  cmul(q, tau, pi); cmul_by_i(q,q); cexp(q,q);
  //fprintf(stderr, "q :\n"); mpc_out_str(stderr, 10, 0, q, MPC_RNDNN); fprintf(stderr, "\n");


  // Determine the number of terms we need
  //     we use our Dupont-like bound
  fset(bound, creal(pi)); fset(ftmp, cimag(tau)); fmul(bound, bound, ftmp);      finv(bound, bound);
  // log_2(e) = log(e)/log(2) = 1/log2
  fset_ui(ftmp, 2); flog(ftmp, ftmp); fmul(bound, bound, ftmp); fmul_ui(bound, bound, p+3);
  fadd(bound, bound, un);
  fsqrt(bound, bound);
  fadd(bound, bound, un);

  // Summing
  for (i=0; i<4; i++) { cone(theta[i]); cone(terms[i]); } cone(theta[4]); cone(theta[5]);
  sign = 1;
  int flag = 0;
  while (fcmp(bound, indice) >= 0)          // >=0 plays the role of ceil(bound)
  { // terms[0] = e^{i pi tau n}, terms[1] = e^{i pi tau n^2}, terms[2] = (e^{2 i n pi z} + e^{-2in pi z}), terms[3] = (e^{2 i (n-1) pi z} + e^{-2i (n-1) pi z})
    // extraterms[0] = (e^{2 i n pi (z+tau/2)} + e^{-2in pi (z+tau/2)}), extraterms[1] = (e^{2 i (n-1) pi (z+tau/2)} + e^{-2i (n-1) pi (z+tau/2)})
    // extraterms[2] = (e^{i n pi tau} + e^{-in pi tau}), extraterms[1] = (e^{i (n-1) pi tau} + e^{-i (n-1) pi tau})

    // n = n+1
    fadd(indice, indice, un);
    //fprintf(stderr, "indice :\n"); mpfr_out_str(stderr, 10, 0, indice, MPC_RNDNN); fprintf(stderr, "\n");

    // Compute the rank n+1
    sign = -sign;

    csqr(tmp, terms[0]); cmul(terms[1], terms[1], tmp); cmul(terms[1], terms[1], q);
    cmul(terms[0], terms[0], q);


    // v0 = 2, v1 = e^(2 i pi z)+e^(-2i pi z)
    if (flag == 0) {
      flag = 1;

      cmul(tmp, z, pi); cmul_2ui(tmp, tmp,1); cmul_by_i(tmp,tmp);
      cexp(v1, tmp); cneg(tmp, tmp); cexp(tmp, tmp); cadd(v1, v1, tmp);
      cset(terms[2], v1); cset_ui(terms[3], 2);
      cset(tmp, v1); cmul(tmp, tmp, terms[1]);
      cadd(theta[0], theta[0], tmp); cadd(theta[1], theta[1], terms[1]); cadd(theta[1], theta[1], terms[1]);
      csub(theta[2], theta[2], tmp); csub(theta[3], theta[3], terms[1]); csub(theta[3], theta[3], terms[1]);

      cset(tmp, tau); cdiv_2ui(tmp, tmp, 1); cadd(tmp, tmp, z); cmul(tmp, tmp, pi); cmul_2ui(tmp, tmp, 1); cmul_by_i(tmp, tmp);
      cexp(extrav1[0], tmp); cneg(tmp, tmp); cexp(tmp, tmp); cadd(extrav1[0], extrav1[0], tmp);
      cset(extraterms[0], extrav1[0]); cset_ui(extraterms[1], 2);
      cmul(tmp, terms[1], extraterms[0]); cadd(theta[4], theta[4], tmp);

      cset(tmp, tau); cmul(tmp, tmp, pi); cmul_by_i(tmp, tmp);
      cexp(extrav1[1], tmp); cneg(tmp, tmp); cexp(tmp, tmp); cadd(extrav1[1], extrav1[1], tmp);
      cset(extraterms[2], extrav1[1]); cset_ui(extraterms[3], 2);
      cmul(tmp, terms[1], extraterms[2]); cadd(theta[5], theta[5], tmp);

      continue;
    }

    // Lucas sequences
    //     (terms[2], terms[3]) = (terms[2]*v1-terms[3], terms[2])
    cmul(tmp, terms[2], v1); csub(tmp, tmp, terms[3]);
    cset(terms[3], terms[2]);
    cset(terms[2], tmp);
    //    (extraterms[0], extraterms[1]) = (extraterms[0]*extrav1[0]-extraterms[1], extraterms[0])
    cmul(tmp, extraterms[0], extrav1[0]); csub(tmp, tmp, extraterms[1]);
    cset(extraterms[1], extraterms[0]);
    cset(extraterms[0], tmp);
    //    (extraterms[2], extraterms[3]) = (extraterms[2]*extrav1[1]-extraterms[3], extraterms[2])
    cmul(tmp, extraterms[2], extrav1[1]); csub(tmp, tmp, extraterms[3]);
    cset(extraterms[3], extraterms[2]);
    cset(extraterms[2], tmp);

    // Add the terms to the result
    cmul(tmp, terms[2], terms[1]);
    cadd(theta[0], theta[0], tmp);   cadd(theta[1], theta[1], terms[1]); cadd(theta[1], theta[1], terms[1]);
    if (sign == 1) {
      cadd(theta[2], theta[2], tmp); cadd(theta[3], theta[3], terms[1]); cadd(theta[3], theta[3], terms[1]);
    } else {
      csub(theta[2], theta[2], tmp); csub(theta[3], theta[3], terms[1]); csub(theta[3], theta[3], terms[1]);
    }
    cmul(tmp, terms[1], extraterms[0]); cadd(theta[4], theta[4], tmp);
    cmul(tmp, terms[1], extraterms[2]); cadd(theta[5], theta[5], tmp);
  }
  // Don't forget the factors in the definition of Theta_10
  // w = e^(i z pi)
  cmul(w, z, pi); cmul_by_i(w,w); cexp(w, w);
  //fprintf(stderr, "w :\n"); mpc_out_str(stderr, 10, 0, w, MPC_RNDNN); fprintf(stderr, "\n");
  // q = e^(i pi tau / 4)
  cmul(q, tau, pi); cmul_by_i(q,q); cdiv_2ui(q, q, 2); cexp(q,q);
  //fprintf(stderr, "q :\n"); mpc_out_str(stderr, 10, 0, q, MPC_RNDNN); fprintf(stderr, "\n");
  cmul(theta[4], theta[4], q);  cmul(theta[4], theta[4], w);
  cmul(theta[5], theta[5], q);
  
  for (i=0; i<6; i++) { cset(res[i], theta[i]); }
  
  

  // clears
  cclear(myz); cclear(mytau);
  cclear(q); cclear(w);
  cclear(pi); cclear(v1);
  cclear(tmp);
  for (i=0; i<4; i++) cclear(terms[i]);
  for (i=0; i<4; i++) cclear(extraterms[i]);
  for (i=0; i<2; i++) cclear(extrav1[i]);
  for (i=0; i<6; i++) cclear(theta[i]);
  fclear(bound);
  fclear(norme);
  fclear(ftmp);
  fclear(indice);
  fclear(un);


  return;
}











/*
 * Returns (in that order) Theta00(0,tau), Theta01(0,tau) with precision p
 */


void NaiveThetaConstants_00_01(mpc_t res[2], mpc_t tau, unsigned long int p)
{
  mpc_t mytau, q, pi, tmp, terms[2], theta[2];
  mpfr_t bound, norme, ftmp, indice, un;
  int sign;
  // inits
  int i; for (i=0; i<2; i++) { cinit(res[i], p); }
  int myprec = 2*p;	// in pratice p + log B + 7 is enough
  cinit(mytau, myprec); cset(mytau, tau);
  cinit(tmp, myprec); czero(tmp);
  cinit(q, myprec);
  cinit(pi,myprec);
  for (i=0; i<2; i++) cinit(terms[i], myprec);
  for (i=0; i<2; i++) cinit(theta[i], myprec);

  finit(bound, myprec);
  finit(norme, myprec);
  finit(ftmp, myprec);
  finit(indice, myprec); fzero(indice);
  finit(un, myprec); fone(un);

  // pi
  cone(pi); cneg(pi, pi); clog(pi, pi); cmul_by_i(pi, pi);
  //fprintf(stderr, "pi :\n"); mpc_out_str(stderr, 10, 0, pi, MPC_RNDNN); fprintf(stderr, "\n");
  //fprintf(stderr, "w :\n"); mpc_out_str(stderr, 10, 0, w, MPC_RNDNN); fprintf(stderr, "\n");
  // q = e^(i pi tau)
  cmul(q, tau, pi); cmul_by_i(q,q); cexp(q,q);
  //fprintf(stderr, "q :\n"); mpc_out_str(stderr, 10, 0, q, MPC_RNDNN); fprintf(stderr, "\n");



  // If tau is fairly big, we can have theta = 1.00000000000000000000thing, and our approximation will be 1 ;
  //         then when we go to the finitedifferences, we find theta10(0)=0, and div by 0 = bye bye
  // So here, if we detect that the first term is really small, we call the function with a bigger precision so that there's at least 10 digits after the 1
  // CalculThetas won't mind, it starts with prec = whatever was sent back
  cset(tmp, q); cnorm(norme, tmp); fsqrt(norme, norme); fone(bound); fdiv_2ui(bound, bound, p);
  if ( fcmp(norme, bound) < 0) {
    // the second term is too small to be seen at that precision
    // we need to increase the precision so that there's at least a few digits at the end of the result
    mpfr_t log; finit(log, myprec); fset_ui(log, 2); flog(log, log); flog(norme, norme); fdiv(log, norme, log); fneg(log, log);
    //fprintf(stderr, "That tau is big!\n");
    NaiveThetaConstants_00_01(res, tau, fget_ui(log) + p +1);
    return ;
  }

  // Determine the number of terms we need
  //     we use our Dupont-like bound
  fset(bound, creal(pi)); fset(ftmp, cimag(tau)); fmul(bound, bound, ftmp);      finv(bound, bound);
  // log_2(e) = log(e)/log(2) = 1/log2
  fset_ui(ftmp, 2); flog(ftmp, ftmp); fmul(bound, bound, ftmp); fmul_ui(bound, bound, p+3);
  fadd(bound, bound, un);
  fsqrt(bound, bound);
  fadd(bound, bound, un);

  // Summing
  for (i=0; i<2; i++) { cone(theta[i]); cone(terms[i]); }
  sign = 1;
  while (fcmp(bound, indice) >= 0)          // >=0 plays the role of ceil(bound)
  { // terms[0] = e^{i pi tau n}, terms[1] = e^{i pi tau n^2}
    // n = n+1
    fadd(indice, indice, un);
    //fprintf(stderr, "indice :\n"); mpfr_out_str(stderr, 10, 0, indice, MPC_RNDNN); fprintf(stderr, "\n");

    // Compute the rank n+1
    sign = -sign;
    csqr(tmp, terms[0]); cmul(terms[1], terms[1], tmp); cmul(terms[1], terms[1], q);
    cmul(terms[0], terms[0], q);
    // Add the terms to the result
    cadd(theta[0], theta[0], terms[1]); cadd(theta[0], theta[0], terms[1]);
    if (sign == 1) {
      cadd(theta[1], theta[1], terms[1]); cadd(theta[1], theta[1], terms[1]);
    } else {
      csub(theta[1], theta[1], terms[1]); csub(theta[1], theta[1], terms[1]);
    }
  }
  for (i=0; i<2; i++) { cset(res[i], theta[i]); }

  // clears
  cclear(mytau);
  cclear(q);
  cclear(pi);
  cclear(tmp);
  for (i=0; i<2; i++) { cclear(terms[i]);}
  for (i=0; i<2; i++) { cclear(theta[i]); }
  fclear(bound);
  fclear(norme);
  fclear(ftmp);
  fclear(indice);
  fclear(un);

  return;
}

