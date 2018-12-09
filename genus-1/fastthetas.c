/* fastthetas.c -- fast computation of theta{00,01,10}({0,z},tau)
 *
 */




#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>

#include "macros.h"
#include "naivethetas.h"
#include "fastthetas.h"
#include "fastthetaconstants.h"
#include "misc.h"


/*
 * Given (x_n, y_n, z_n, t_n), compute the next step of the AGMPrime.
 *         (\sqrt(x_nz_n)+\sqrt(y_nt_n))/2, etc.
 * a is modified in place.
 * No attention given to sign.
 *
 */

void AGMPrime(mpc_t a[4])
{
  int prec;
  prec = cprec(a[0]);

  mpc_t tmp[4], prod[4];
  int i; for (i=0; i<4; i++) { cinit(tmp[i], prec); cinit(prod[i], prec); }
  
  csqrt(a[0], a[0]);
  csqrt(a[1], a[1]);
  csqrt(a[2], a[2]);
  csqrt(a[3], a[3]);

   // Optimizing the number of products
   //   (x+y)(z+t) = xz + yt + xt + yz
   //   (x-y)(z-t) = xz + yt - xt - yz
  //   (a+b)^2 + (a-b)^2 = 2(a^2+b^2)
  // Gain : 20% !
  cadd(tmp[0], a[0], a[1]);
  cadd(tmp[1], a[2], a[3]);
  csub(tmp[2], a[0], a[1]);
  csub(tmp[3], a[2], a[3]);

  cmul(prod[0], tmp[0], tmp[1]);
  cmul(prod[1], tmp[2], tmp[3]);
  csqr(prod[2], tmp[1]);
  csqr(prod[3], tmp[3]);

  cadd(a[0], prod[0], prod[1]); cdiv_2ui(a[0], a[0], 2);
  csub(a[1], prod[0], prod[1]); cdiv_2ui(a[1], a[1], 2);
  cadd(a[2], prod[2], prod[3]); cdiv_2ui(a[2], a[2], 2);
  csub(a[3], prod[2], prod[3]); cdiv_2ui(a[3], a[3], 2);

  //int i; for (i=0; i<4; i++)  { fprintf(stderr, "a[%d] = ", i); mpc_out_str(stderr, 10, 0, a[i], MPC_RNDNN); fprintf(stderr, "\n"); }

  // Clears
  for (i=0; i<4; i++) {cclear(prod[i]); cclear(tmp[i]); }
  return ;
}




/*
 * Helper function that computes x^2^n
 *
*/

void Pow2PowN(mpc_t r, mpc_t x, unsigned long int n)
{
  unsigned long int i = 0;
  int prec;
  prec = cprec(x);
  cinit(r, prec);
  cset(r, x);

  //fprintf(stderr, "n = %ul\n", n);
  
  while (n > i)
  {  csqr(r,r);
     //fprintf(stderr, "r = "); mpc_out_str(stderr, 10, 0, r, MPC_RNDNN); fprintf(stderr, "\n");
     i++;
  }

  return;
}


/* Check to see if |arg1-arg2| < bound
 * Returns 1 if yes, 0 if not
 */

int distanceSmallerThan(mpc_t arg1, mpc_t arg2, mpfr_t bound)
{
  unsigned long int prec = fprec(bound);
  mpc_t diff; cinit(diff, prec); csub(diff, arg1, arg2);
  mpfr_t absolute; finit(absolute, prec); cabs(absolute, diff);
  int result;
  if ( fcmp(absolute, bound) <= 0 ) { result=1; } else { result=0; }
  //if (result == 1) { fprintf(stderr, "absolute : "); mpfr_out_str(stderr, 10, 0, absolute, MPC_RNDNN); fprintf(stderr, "\n"); }
  cclear(diff); fclear(absolute); return result;
}



/*
 * Computes iterates of the AGMPrime until the right things are too small, and computes lambda and mu.
*/

void Finfty4args(mpc_t lambda, mpc_t mu, mpc_t a[4])
{
  int prec, i, n;
  prec = cprec(a[0]);  // min(cprec(a[i])) ?
  cinit(lambda, prec);
  cinit(mu, prec);
  
  // article: C = 83.64, e > 1/1808, b = 83.64*2, lambda/mu < 13.71, C' < 55
  // nb steps : log(|log(|z|)) + 1 < 3.44 + 1 (don't forget to look at big *and* small quotients)
  // Pow2PowN : (n+1) log C' +1 - log(C'-1) = 5.81 (logP + 5.44) + 1 - 5.78
  // the approximation we do on lambda loses 55 bits
  int AGMPrec = prec + (6*ceilLog2OfInteger(prec) + 27) + 55;
  // prec loss with n steps : 2n log b bits = 14.77 (\log P + 4.44)
  int fullprec = AGMPrec + (2*15*ceilLog2OfInteger(AGMPrec) + 66);
  
  
  mpfr_t bound; finit(bound, fullprec); fone(bound); fdiv_2ui(bound, bound, AGMPrec);
  mpc_t highPrecA[4]; for (i=0; i<4; i++) { cinit(highPrecA[i], fullprec); cset(highPrecA[i], a[i]); }
  mpc_t highPrecbuf; cinit(highPrecbuf, fullprec);

  // Computing Finfty
  n=0;
  while ( distanceSmallerThan(highPrecA[2], highPrecA[3], bound) == 0 )
  { n++;
    AGMPrime(highPrecA);
  }
  // one more time
  AGMPrime(highPrecA); n++;
  cset(mu, highPrecA[2]);
  // lambda, mu
  cdiv(highPrecbuf, highPrecA[0], highPrecA[2]);
  //fprintf(stderr, "quotient = "); mpc_out_str(stderr, 10, 0, highPrecbuf, MPC_RNDNN); fprintf(stderr, "\n");
  //fprintf(stderr, "power = %d\n", i);
  mpc_t res; cinit(res, fullprec);
  Pow2PowN(res, highPrecbuf, n);
  //fprintf(stderr, "quotientpow = "); mpc_out_str(stderr, 10, 0, res, MPC_RNDNN); fprintf(stderr, "\n");
  cmul(res, res, highPrecA[2]);
  cset(lambda, res);

  cclear(res); cclear(highPrecbuf);
  fclear(bound);
  for (i=0; i<4; i++) { cclear(highPrecA[i]); }
  //for (j=0; j<4; j++) {cclear(b[j]); cclear(diff[j]); fclear(norm[j]);}  fclear(bound);
  return ;
}


/*
 * Computes lambda and mu corresponding to Finfty([1, quo1, 1, quo2])
 *
*/

void Finfty2args(mpc_t lambda, mpc_t mu, mpc_t quo1, mpc_t quo2)
{
  mpc_t b[4];
  int prec, i;

  prec = cprec(quo1);  // min(cprec(quo1),cprec(quo2)) ?
  cinit(lambda, prec);
  cinit(mu, prec);

  for (i=0; i<4; i++)
    cinit(b[i], prec);
  cone(b[0]);
  cset(b[1], quo1);
  cone(b[2]);
  cset(b[3], quo2);
  //fprintf(stderr, "b = "); for (i=0; i<4; i++) { mpc_out_str(stderr, 10, 0, b[i], MPC_RNDNN); fprintf(stderr, "\n"); }
  
  Finfty4args(lambda, mu, b);

  for (i=0; i<4; i++)
    cclear(b[i]);

  return ;
}


/*
 * Given quotients theta01^2(z,tau)/theta00^2(z,tau) and theta01^2(0,tau)/theta00^2(z,tau), compute lambda and mu.
 *
*/

void WholeFunction(mpc_t lambda, mpc_t mu, mpc_t quo1, mpc_t quo2)
{
  mpc_t k, quo10_00_z, aux, b[2];
  int prec, i;
  prec = cprec(quo1);
  cinit(lambda, prec);
  cinit(mu, prec);
  cinit(k, prec);
  cinit(quo10_00_z, prec);
  cinit(aux, prec);
  for (i=0; i<2; i++) cinit(b[i], prec);

  // k^2(tau) = sqrt(1-(theta01^2/theta00^2)^2)       (Jacobi)
  cone(k); csqr(aux, quo2); csub(k, k, aux); csqrt(k, k);

  // theta10^2(z)/theta00^2(z) = [1-(theta01(z)/theta00(z))^2*k'(t)]/k(t)
  cone(quo10_00_z); cmul(aux, quo1, quo2);
  csub(quo10_00_z, quo10_00_z, aux); cdiv(quo10_00_z, quo10_00_z, k);
  
  // recover lambda/theta00^2(z), mu/theta00^2(0)
  Finfty2args(lambda, mu, quo10_00_z, k);
  // recover 1/theta00^2(z), 1/theta00^2(0);
  Finfty2args(b[0], b[1], quo1, quo2);
  // recover lambda and mu
  cdiv(lambda, lambda, b[0]); cdiv(mu, mu, b[1]);
  //fprintf(stderr, "lambda :"); mpc_out_str(stderr, 10, 0, lambda, MPC_RNDNN); fprintf(stderr, "\n");
  //fprintf(stderr, "mu :"); mpc_out_str(stderr, 10, 0, mu, MPC_RNDNN); fprintf(stderr, "\n"); 

  // clears
  cclear(k);
  cclear(quo10_00_z);
  cclear(aux);
  for (i=0; i<2; i++) cclear(b[i]);
}


/*
 * One step of Newton with finite differences, applied to the WholeFunction.
 * Give quo1(z) quo2(0) at precision p, and lambda mu at precision 2*p
 * Receive quo1(z) quo2(0) at precision 2p
*/

void NewtonStep(mpc_t quo2p[2], mpc_t quop[2], mpc_t lambda, mpc_t mu)
{
  mpc_t origquo[2], epsilon, difffinies[6], tmp, jacobian[4], chang[2];
  int prec, i;
  // inits
  prec = cprec(quop[0]);
  for (i=0; i<2; i++) {
    cinit(origquo[i], 2*prec); cset(origquo[i], quop[i]);
    cinit(quo2p[i], 2*prec);
    cinit(chang[i], 2*prec);
  }
  cinit(tmp, 2*prec);
  for (i=0; i<6; i++)
    cinit(difffinies[i], 2*prec);
  for (i=0; i<4; i++)
    cinit(jacobian[i], 2*prec);

  // Compute epsilon
  cinit(epsilon, 2*prec);
  cone(epsilon);
  cdiv_2ui(epsilon, epsilon, prec);

  // Finites Differences
  cadd(tmp, origquo[0], epsilon);
  WholeFunction(difffinies[0], difffinies[1], tmp, origquo[1]);
  cadd(tmp, origquo[1], epsilon);
  WholeFunction(difffinies[2], difffinies[3], origquo[0], tmp);
  WholeFunction(difffinies[4], difffinies[5], origquo[0], origquo[1]);

  //fprintf(stderr, "Jacobian :");
  // Inverse of the Jacobian matrix
  //    because F_y is a function of y, not x, the upper right corner of the jacobian is 0
  // J(1,1) = 1/(Df_x/Dx)
  csub(jacobian[0], difffinies[0], difffinies[4]);
  cmul_2ui(jacobian[0], jacobian[0], prec);    // times 1/epsilon
  //mpc_out_str(stderr, 10, 0, jacobian[0], MPC_RNDNN); fprintf(stderr, "\n");
  cinv(jacobian[0], jacobian[0]);
  // J(2,2) = 1/(Df_y/Dy)
  csub(jacobian[3], difffinies[3], difffinies[5]);
  cmul_2ui(jacobian[3], jacobian[3], prec);
  //mpc_out_str(stderr, 10, 0, jacobian[3], MPC_RNDNN); fprintf(stderr, "\n");
  cinv(jacobian[3], jacobian[3]);
  // J(2,1) = (Df_x/Dy) * J(1,1) * J(2,2)
  csub(jacobian[2], difffinies[2], difffinies[4]);
  cmul_2ui(jacobian[2], jacobian[2], prec);
  //mpc_out_str(stderr, 10, 0, jacobian[2], MPC_RNDNN); fprintf(stderr, "\n");
  cneg(jacobian[2], jacobian[2]);
  cmul(jacobian[2], jacobian[2], jacobian[0]);
  cmul(jacobian[2], jacobian[2], jacobian[3]);
  //fprintf(stderr, "Jacobian : "); for (k=0; k<4; k++) {mpc_out_str(stderr, 10, 0, jacobian[k], MPC_RNDNN); fprintf(stderr, "\n");}

  // Change vector
  csub(chang[0], difffinies[4], lambda);
  csub(chang[1], difffinies[5], mu);
  //fprintf(stderr, "change before newton : "); mpc_out_str(stderr, 10, 0, chang[0], MPC_RNDNN); fprintf(stderr, "\n"); mpc_out_str(stderr, 10, 0, chang[1], MPC_RNDNN); fprintf(stderr, "\n");

  // Apply the Jacobian : x = x*J(1,1) + y*J(2,1) ; y = y*J(2,2)
  cmul(tmp, chang[0], jacobian[0]);
  cmul(chang[0], chang[1], jacobian[2]);
  cadd(chang[0], chang[0], tmp);
  cmul(chang[1], chang[1], jacobian[3]);

  // x - (F(x)-l)J^-1
  csub(quo2p[0], origquo[0], chang[0]);
  csub(quo2p[1], origquo[1], chang[1]); 

  // clears
  for (i=0; i<2; i++) {
    cclear(origquo[i]);
    cclear(chang[i]);
  }
  for (i=0; i<4; i++)
    cclear(jacobian[i]);
  for (i=0; i<6; i++)
    cclear(difffinies[i]);
  cclear(tmp);
  cclear(epsilon);

  return;
}






/*
 * Fast computation of the theta functions.
 * Give z and tau at precision p
 * Receive [theta00(z,tau) theta00(0,tau) theta01(z,tau) theta01(0,tau)]
*/

void MathfrakF(mpc_t thetas[4], mpc_t z, mpc_t tau)
{
  mpc_t tmp[2], mu, lambda, pi;
  int goalprec, prec, i;

  goalprec = cprec(z);

  // inits
  cinit(mu, goalprec);
  cinit(lambda, goalprec);
  cinit(pi, goalprec);
  for (i=0; i<4; i++)
    cinit(thetas[i], goalprec);

  // Determine which lowprec you want (to avoid threshold effects caused by Newton)
  //   If Newton gave you exactly twice as many exact digits, you'd just do p/2+1
  //      but really it's "2P-delta", so we do p/2+delta+1
  // Experiments: delta=2 works mostly well, delta=4 works all the time
  int delta = 4;
  int myLowPrec = goalprec;
  while (myLowPrec > LOW_PREC_THETA)   { myLowPrec = myLowPrec/2 + delta/2 + 1; }

  // Compute thetas at low prec
  mpc_t smallprecthetas[4]; int k; for (k=0; k<4; k++) cinit(smallprecthetas[k], myLowPrec);
  InternalNaiveThetaFunctionAndThetaConstants_00_01(smallprecthetas, z, tau, myLowPrec);

  myLowPrec = cprec(smallprecthetas[0]);    // the prec may have changed if it was too small to fit the first term of the series

  if (myLowPrec >= goalprec)
  { // No need to keep going, we have the right result
    for (k=0; k<4; k++) { csqr(thetas[k], smallprecthetas[k]); }
  } else {
    // Need to Newton

    mpc_t quotients[2], musmallprec, lambdasmallprec;


    // Computing pi
        //  mpfr_const_pi(pi);      // à caster
    // currently we get it as the solution to e^x = -1
    cone(pi); cneg(pi, pi);
    clog(pi, pi); cmul_by_i(pi, pi);

    // lambda and mu
    cneg(mu, tau);
    cmul_by_i(mu, mu);
    cinv(mu, mu);

    csqr(lambda, z);
    cdiv(lambda, lambda, tau);
    cneg(lambda, lambda);   cmul_by_i(lambda, lambda);
    cmul(lambda, lambda, pi);   cmul_2ui(lambda, lambda, 1);
    cexp(lambda, lambda);
    cmul(lambda, lambda, mu);

    // theta01^2(z)/theta00^2(z) and theta01^2(0)/theta00^2(0)
    cinit(tmp[0], myLowPrec); cinit(tmp[1], myLowPrec);
    cinit(quotients[0], myLowPrec); cinit(quotients[1], myLowPrec);
    cdiv(quotients[0], smallprecthetas[2], smallprecthetas[0]);
    cdiv(quotients[1], smallprecthetas[3], smallprecthetas[1]);
    csqr(quotients[0], quotients[0]); csqr(quotients[1], quotients[1]);
    for (k=0; k<4; k++) cclear(smallprecthetas[k]);

    //fprintf(stderr, "quotients :"); int k; for (k=0; k<2; k++) {mpc_out_str(stderr, 10, 0, quotients[k], MPC_RNDNN); fprintf(stderr, "\n");}

    // Newton
    prec = myLowPrec;
    cinit(musmallprec, prec); cinit(lambdasmallprec, prec);
    while (prec < goalprec)
    { prec = 2*prec;
      // inits
      cset_prec(tmp[0], prec); cset_prec(tmp[1], prec); cset_prec(musmallprec, prec); cset_prec(lambdasmallprec, prec);
      cset(musmallprec, mu); cset(lambdasmallprec, lambda);
      // Diff finies
      NewtonStep(tmp, quotients, lambdasmallprec, musmallprec);
      // copying the results
      cset_prec(quotients[0], prec); cset_prec(quotients[1], prec);
      cset(quotients[0], tmp[0]); cset(quotients[1], tmp[1]);   
      //fprintf(stderr, "quotients :"); int l; for (l=0; l<2; l++) {mpc_out_str(stderr, 10, 0, quotients[l], MPC_RNDNN); fprintf(stderr, "\n");} 
    }
    cset_prec(quotients[0], goalprec); cset_prec(quotients[1], goalprec);
    cset(quotients[0], tmp[0]); cset(quotients[1], tmp[1]);
    cclear(tmp[0]); cclear(tmp[1]);

    //fprintf(stderr, "refined quotients : "); int j; for (j=0; j<2; j++) {mpc_out_str(stderr, 10, 0, quotients[j], MPC_RNDNN); fprintf(stderr, "\n");}

    // Extracting thetas
    // theta00
    Finfty2args(thetas[0], thetas[1], quotients[0], quotients[1]);
    cinv(thetas[0], thetas[0]);
    cinv(thetas[1], thetas[1]);
    // theta01
    cmul(thetas[2], thetas[0], quotients[0]);
    cmul(thetas[3], thetas[1], quotients[1]);
    //fprintf(stderr, "thetas : "); int j; for (j=0; j<4; j++) {mpc_out_str(stderr, 10, 0, thetas[j], MPC_RNDNN); fprintf(stderr, "\n");}

    // clears
    cclear(lambdasmallprec); cclear(musmallprec);
    cclear(quotients[0]); cclear(quotients[1]);
  }

  // clear
  cclear(mu);  cclear(lambda);
  cclear(pi);

  return ;

}






void FastThetas(mpc_t thetas[6], mpc_t z, mpc_t tau)
{
  unsigned long int prec = cprec(tau);
  unsigned long int lossOfPrec=0;
  int i;
  mpfr_t tmp, bound;
  mpc_t pi, tauprime, zprime, subroutine_results[4], tmp_thetas[6], ctmp[4];
  finit(tmp, prec);
  finit(bound, prec);
  cinit(pi, prec);
  // pi
  cone(pi); cneg(pi, pi); clog(pi, pi); cmul_by_i(pi, pi);
  
  // Assume that tau is in F and z is reduced
  fset_ui(tmp, prec); fdiv_ui(tmp, tmp, 25);
  if (fcmp(cimag(tau), tmp) >= 0 ) {       // P < 25 Im(tau)
    
    // this function loses in principle log B + 7 bits
    // Compute B
    fset(bound, creal(pi)); fset(tmp, cimag(tau)); fmul(bound, bound, tmp); finv(bound, bound);
    // log_2(e) = log(e)/log(2) = 1/log2
    fset_ui(tmp, 2); flog(tmp, tmp); fmul(bound, bound, tmp); fmul_ui(bound, bound, prec+3);
    fsqrt(bound, bound); fadd_ui(bound, bound, 1);
    flog(bound, bound); fdiv(bound, bound, tmp); frint(bound, bound);
    fadd_ui(bound, bound, 8);
    lossOfPrec = fget_ui(bound);
             //fprintf(stderr, "lossOfPrec : %lu \n", lossOfPrec);
    
    NaiveThetaFunctionAndThetaConstants_00_01_10(thetas, z, tau, prec+lossOfPrec);
    
    
  } else {

    // Compute s, set tau' = tau/2^{s+1}, z' = z/2^{s+2}
    unsigned long int s;
    cabs(bound, tau); flog(bound, bound);
    fset_ui(tmp, 2); flog(tmp, tmp);
    fdiv(bound, bound, tmp); mpfr_floor(bound, bound);
    s = fget_ui(bound); // fear not, it's smaller than p, and p fits in an unsigned long int too
             //fprintf(stderr, "s : %lu \n", s);
    
    // Compute the precision needed as per the paper
      // if you count right, we lose 30.48s + 23.8 Im(tau) + 31
      // leave it like that (ie don't bound by P), it's more efficient
	// the +36 is there because our quantities have integral part of maximal size log(C^2/epsilon^2) < 35
    fset(bound, cimag(tau)); mpfr_floor(bound, bound);
    unsigned long int mathcalPrec = prec + (31*s + 24*(fget_ui(bound)+1) + 31) + 36;
    // unsigned long int mathcalPrec = 2*prec+36
             //fprintf(stderr, "mathcalPrec: %lu\n", mathcalPrec);
    
    // Compute the precision lost during the call to Finfty (see paper)
    //lossOfPrec = 55*ceilLog2OfInteger(prec) + 3;
    lossOfPrec = 55*ceilLog2OfInteger(mathcalPrec) + 104;
             //fprintf(stderr, "lossOfPrec: %lu\n", lossOfPrec);
    
    // s+2 because we don't want to lose any digits when dividing by 2^s then multiplying back by it
    cinit(tauprime, mathcalPrec+(s+2)+lossOfPrec); cinit(zprime, mathcalPrec+(s+2)+lossOfPrec);
    cdiv_2ui(tauprime, tau, s+1); cdiv_2ui(zprime, z, s+2);
             //fprintf(stderr, "tau_2 : "); mpc_out_str(stderr, 10, 0, tauprime, MPC_RNDNN); fprintf(stderr, "\n");
             //fprintf(stderr, "z_2 : "); mpc_out_str(stderr, 10, 0, zprime, MPC_RNDNN); fprintf(stderr, "\n");
    
    // Call the subroutine
    MathfrakF(subroutine_results, zprime, tauprime);    // only the first 4 are good
    // at this point the thetas should be exact with absolute precision mathcalPrec
             //fprintf(stderr, "thetas : "); for (i=0; i<4; i++) { mpc_out_str(stderr, 10, 0, tmp_thetas[i], MPC_RNDNN); fprintf(stderr, "\n");}
    // Copy the results
    for (i=0; i<6; i++) { cinit(tmp_thetas[i], mathcalPrec); }
    for (i=0; i<4; i++) { cset(tmp_thetas[i], subroutine_results[i]); }
    
    // Step 8 : tau doubling
    for (i=0; i<4; i++) { cinit(ctmp[i], mathcalPrec); }
    // no need to check for sign here, since Re(theta_{00,01}(z,tau)) > 0 for Im(z) < Im(tau)/4 and Im(tau) > 0.345 (which is the case)
    for (i=0; i<4; i++) { csqrt(tmp_thetas[i], tmp_thetas[i]); }
    cadd(ctmp[0], tmp_thetas[0], tmp_thetas[2]);    // x+y
    csub(ctmp[1], tmp_thetas[0], tmp_thetas[2]);    // x-y
    cadd(ctmp[2], tmp_thetas[1], tmp_thetas[3]);    // z+t
    csub(ctmp[3], tmp_thetas[1], tmp_thetas[3]);    // z-t
    // theta constants
    csqr(tmp_thetas[1], ctmp[2]);	// (z+t)^2
    csqr(tmp_thetas[3], ctmp[3]);	// (z-t)^2
    cadd(tmp_thetas[1], tmp_thetas[1], tmp_thetas[3]);		// (z+t)^2+(z-t)^2
    cmul_2ui(tmp_thetas[3], tmp_thetas[3], 1); csub(tmp_thetas[3], tmp_thetas[1], tmp_thetas[3]);   // (z+t)^2-(z-t)^2
    cmul(tmp_thetas[5], ctmp[2], ctmp[3]); // (z+t)(z-t)
    cdiv_2ui(tmp_thetas[1], tmp_thetas[1], 2); cdiv_2ui(tmp_thetas[3], tmp_thetas[3], 2); cdiv_2ui(tmp_thetas[5], tmp_thetas[5], 1);
    // theta functions
    cmul(tmp_thetas[0], ctmp[0], ctmp[2]);   // (x+y)(z+t)
    cmul(tmp_thetas[2], ctmp[1], ctmp[3]);   // (x-y)(z-t)
    cadd(tmp_thetas[0], tmp_thetas[0], tmp_thetas[2]);	// (x+y)(z+t)+(x-y)(z-t)
    cmul_2ui(tmp_thetas[2], tmp_thetas[2], 1); csub(tmp_thetas[2], tmp_thetas[0], tmp_thetas[2]);   // (x+y)(z+t)-(x-y)(z-t)
    cdiv_2ui(tmp_thetas[0], tmp_thetas[0], 2); cdiv_2ui(tmp_thetas[2], tmp_thetas[2], 2);
    // theta10
    cmul(tmp_thetas[4], ctmp[0], ctmp[3]);   // (x+y)(z-t) 
    cmul(ctmp[3], ctmp[1], ctmp[2]);   // (x-y)(z+t)
    cadd(tmp_thetas[4], tmp_thetas[4], ctmp[3]); cdiv_2ui(tmp_thetas[4], tmp_thetas[4], 2);
             //fprintf(stderr, "thetas : "); for (i=0; i<5; i++) { mpc_out_str(stderr, 10, 0, tmp_thetas[i], MPC_RNDNN); fprintf(stderr, "\n");}
    
    // Step 9: square root of theta constants
    //     no need for a small precision approximation: Re(theta00(0)) > 0 and Re(theta01(0)) > 0, always
    //     we don't take the sqr of theta10 because:
    //             - if we go into the loop: we don't need it and it'll get replaced by theta10^2(0,tau)
    //             - if we don't go into the loop: we will need to use the equation of the variety
    csqrt(tmp_thetas[1], tmp_thetas[1]);
    csqrt(tmp_thetas[3], tmp_thetas[3]);
    
    // Step 10: z doubling
    csqr(ctmp[0], tmp_thetas[0]);
    csqr(ctmp[1], tmp_thetas[2]);
    csqr(ctmp[2], tmp_thetas[4]);
    cadd(tmp_thetas[0], ctmp[1], ctmp[2]);
    csub(tmp_thetas[2], ctmp[0], ctmp[2]);
    // third powers of theta constants
    csqr(ctmp[0], tmp_thetas[1]); cmul(ctmp[0], ctmp[0], tmp_thetas[1]);
    csqr(ctmp[1], tmp_thetas[3]); cmul(ctmp[1], ctmp[1], tmp_thetas[3]);
    cdiv(tmp_thetas[0], tmp_thetas[0], ctmp[0]);
    cdiv(tmp_thetas[2], tmp_thetas[2], ctmp[1]);
             //fprintf(stderr, "thetas : "); for (i=0; i<6; i++) { mpc_out_str(stderr, 10, 0, tmp_thetas[i], MPC_RNDNN); fprintf(stderr, "\n");}
    // At this point we have [t00(z1/2, tau_1), t00(0,tau_1), t01(z1/2,tau_1), t01(0,tau_1), t10^2(z_2,tau_1), t10^2(0,tau_1)]
    
    for (i=1; i<= s; i++) {
      // Step 12-15
      cadd(ctmp[0], tmp_thetas[0], tmp_thetas[2]);    // x+y
      csub(ctmp[1], tmp_thetas[0], tmp_thetas[2]);    // x-y
      cadd(ctmp[2], tmp_thetas[1], tmp_thetas[3]);    // z+t
      csub(ctmp[3], tmp_thetas[1], tmp_thetas[3]);    // z-t
	// Step 12
      csqr(tmp_thetas[1], ctmp[2]);	// (z+t)^2
      csqr(tmp_thetas[3], ctmp[3]);	// (z-t)^2
      cadd(tmp_thetas[1], tmp_thetas[1], tmp_thetas[3]);		// (z+t)^2+(z-t)^2
      cmul_2ui(tmp_thetas[3], tmp_thetas[3], 1); csub(tmp_thetas[3], tmp_thetas[1], tmp_thetas[3]);   // (z+t)^2-(z-t)^2
      cdiv_2ui(tmp_thetas[1], tmp_thetas[1], 2); cdiv_2ui(tmp_thetas[3], tmp_thetas[3], 2);
	// Step 13
      cmul(tmp_thetas[0], ctmp[0], ctmp[2]);   // (x+y)(z+t)
      cmul(tmp_thetas[2], ctmp[1], ctmp[3]);   // (x-y)(z-t)
      cadd(tmp_thetas[0], tmp_thetas[0], tmp_thetas[2]);	// (x+y)(z+t)+(x-y)(z-t)
      cmul_2ui(tmp_thetas[2], tmp_thetas[2], 1); csub(tmp_thetas[2], tmp_thetas[0], tmp_thetas[2]);   // (x+y)(z+t)-(x-y)(z-t)
      cdiv_2ui(tmp_thetas[0], tmp_thetas[0], 2); cdiv_2ui(tmp_thetas[2], tmp_thetas[2], 2);
	// Step 14
      if (i==s) {
	cmul(tmp_thetas[5], ctmp[2], ctmp[3]);   // (z+t)(z-t)
	cdiv_2ui(tmp_thetas[5], tmp_thetas[5], 1);
      }
	// Step 15
      cmul(tmp_thetas[4], ctmp[0], ctmp[3]);   // (x+y)(z-t) 
      cmul(ctmp[3], ctmp[1], ctmp[2]);   // (x-y)(z+t)
      cadd(tmp_thetas[4], tmp_thetas[4], ctmp[3]); cdiv_2ui(tmp_thetas[4], tmp_thetas[4], 2);
      
      // Step 16			// no need to check sign: Re(theta00)>0 and Re(theta01)>0
      csqrt(tmp_thetas[1], tmp_thetas[1]);
      csqrt(tmp_thetas[3], tmp_thetas[3]);
      // leave theta10(0) as it is because we need it in step 19 - no sense in sqrt to sqr later
      
      // Step 17
      csqr(ctmp[0], tmp_thetas[0]);
      csqr(ctmp[1], tmp_thetas[2]);
      csqr(ctmp[2], tmp_thetas[4]);
      cadd(tmp_thetas[0], ctmp[1], ctmp[2]);
      csub(tmp_thetas[2], ctmp[0], ctmp[2]);
      // third powers of theta constants
      csqr(ctmp[0], tmp_thetas[1]); cmul(ctmp[0], ctmp[0], tmp_thetas[1]);
      csqr(ctmp[1], tmp_thetas[3]); cmul(ctmp[1], ctmp[1], tmp_thetas[3]);
      cdiv(tmp_thetas[0], tmp_thetas[0], ctmp[0]);
      cdiv(tmp_thetas[2], tmp_thetas[2], ctmp[1]);
    }
    // we now have [theta(z/2,tau), theta(0,tau), theta01(z/2,tau), theta01(0,tau), theta10^2(z/4,tau), theta10^2(0,tau)]
             //fprintf(stderr, "thetas : "); for (i=0; i<6; i++) { mpc_out_str(stderr, 10, 0, tmp_thetas[i], MPC_RNDNN); fprintf(stderr, "\n");}
    
    // Step 19: equation of the variety
    cmul(ctmp[0], tmp_thetas[0], tmp_thetas[1]); csqr(ctmp[0], ctmp[0]);
    cmul(ctmp[1], tmp_thetas[2], tmp_thetas[3]); csqr(ctmp[1], ctmp[1]);
    csub(ctmp[0], ctmp[0], ctmp[1]);
    cdiv(tmp_thetas[4], ctmp[0], tmp_thetas[5]); 
       // Now we can take the sqrt of theta10^2(0); and no need to check the sign here: Re(theta10(0)) = Re(2q^{1/4}) = cos(piRe(tau)/4) > cos(pi/8) > 0
    csqrt(tmp_thetas[5], tmp_thetas[5]);
    // we now have [theta(z/2,tau), theta(0,tau), theta01(z/2,tau), theta01(0,tau), theta10^2(z/2,tau), theta10(0,tau)]
             //fprintf(stderr, "thetas : "); for (i=0; i<6; i++) { mpc_out_str(stderr, 10, 0, tmp_thetas[i], MPC_RNDNN); fprintf(stderr, "\n");}    
    
    // Step 20: last z-duplication
    csqr(ctmp[0], tmp_thetas[0]); csqr(ctmp[0], ctmp[0]);
    csqr(ctmp[1], tmp_thetas[2]); csqr(ctmp[1], ctmp[1]);
    csqr(ctmp[2], tmp_thetas[4]);
    cadd(tmp_thetas[0], ctmp[1], ctmp[2]);
    csub(tmp_thetas[2], ctmp[0], ctmp[2]);
    csub(tmp_thetas[4], ctmp[0], ctmp[1]);
    csqr(ctmp[0], tmp_thetas[1]); cmul(ctmp[0], ctmp[0], tmp_thetas[1]);
    csqr(ctmp[1], tmp_thetas[3]); cmul(ctmp[1], ctmp[1], tmp_thetas[3]);
    csqr(ctmp[2], tmp_thetas[5]); cmul(ctmp[2], ctmp[2], tmp_thetas[5]);
    cdiv(tmp_thetas[0], tmp_thetas[0], ctmp[0]);
    cdiv(tmp_thetas[2], tmp_thetas[2], ctmp[1]);
    cdiv(tmp_thetas[4], tmp_thetas[4], ctmp[2]);
    
    
    // The end!  		precision p+2 because theta is bounded by 4
    for (i=0; i<6; i++) { cinit(thetas[i], prec+2); cset(thetas[i], tmp_thetas[i]); }
    
    // Clears
    cclear(tauprime); cclear(zprime);
    for (i=0; i<4; i++) { cclear(ctmp[i]); }
    for (i=0; i<6; i++) { cclear(tmp_thetas[i]); }
  }
    

    // Clears
    fclear(tmp); fclear(bound); cclear(pi);
}




