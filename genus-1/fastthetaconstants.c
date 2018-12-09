/* fastthetaconstants.c -- fast computation of theta{00,01,10}(0,tau) using Dupont's algorithm
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
#include "fastthetaconstants.h"
#include "misc.h"




/*
  DupontOneStepOfAGM - computes one step of AGM, assuming the choice of root is good (ok if tau in Fk')
 */

void DupontOneStepOfAGM(mpc_t a[2]) {
  int prec;
  prec = cprec(a[0]);

  mpc_t aux;
  cinit(aux, prec);

  cadd(aux, a[0], a[1]); cdiv_2ui(aux, aux, 1);
  cmul(a[1], a[0], a[1]); csqrt(a[1], a[1]);
  cset(a[0], aux);
  // Note: apparently the choice of sign when we take that square root is always good... must be a property of F_k' ?

  cclear(aux);
  return;
}

/*
  DupontM - computes AGM(1,z)
 */
void DupontM(mpc_t res, mpc_t quo) {
  int prec;
  prec = cprec(quo);
  cinit(res, prec);

  mpc_t a[2]; cinit(a[0], prec); cinit(a[1], prec);
  cone(a[0]);
  cset(a[1], quo);

  // Computing Finfty
  while ( (creldist(a[0],a[1]) < prec+1) )                // <- this may make us lose absolute precision
  { DupontOneStepOfAGM(a);
  }

  cset(res, a[0]);

  cclear(a[0]); cclear(a[1]);
  return;
}

/*
  DupontFtau - computes iM(z) - tau M(sqrt(1-z^2))
 */
void DupontFtau(mpc_t res, mpc_t arg, mpc_t tau) {
  int prec;
  prec = cprec(arg);
  cinit(res, prec);
  mpc_t tmp, aux; cinit(tmp, prec); cinit(aux, prec);

  // iM(z)
  DupontM(res, arg);
  cmul_by_i(res, res);
  // tau M(sqrt(1-z^2))
  cone(aux); csqr(tmp, arg);
  csub(aux, aux, tmp);
  csqrt(aux, aux);
  DupontM(tmp, aux);
  cmul(tmp, tmp, tau);

  csub(res, res, tmp);

  cclear(tmp); cclear(aux);
  return;
}

/*
  DupontNewtonStep - given k' with prec p and tau with prec 2p, computes k' with prec 2p
 */
void DupontNewtonStep(mpc_t quo2p, mpc_t quop, mpc_t tau) {
  mpc_t epsilon, tmp, aux, origquo;
  int prec;
  prec = cprec(quop);

  cinit(quo2p, 2*prec);
  cinit(origquo, 2*prec); cset(origquo, quop);
  //fprintf(stderr, "origquo : \n"); mpc_out_str(stderr, 10, 0, origquo, MPC_RNDNN); fprintf(stderr, "\n");
  
  // Compute epsilon
  cinit(epsilon, 2*prec);
  cone(epsilon); cdiv_2ui(epsilon, epsilon, prec);
  //fprintf(stderr, "epsilon : \n"); mpc_out_str(stderr, 10, 0, epsilon, MPC_RNDNN); fprintf(stderr, "\n");

  // Finite differences times f : f(x) * e /(f(x+e)-f(x))
  cinit(tmp, 2*prec); cinit(aux, 2*prec);
  cadd(tmp, origquo, epsilon);
  DupontFtau(aux, tmp, tau);
  DupontFtau(tmp, origquo, tau);
  csub(aux, aux, tmp);
  cinv(aux, aux);
  cmul(aux, aux, tmp); cmul(aux, aux, epsilon);
  //fprintf(stderr, "newton corrective term : \n"); mpc_out_str(stderr, 10, 0, aux, MPC_RNDNN); fprintf(stderr, "\n");

  cset(quo2p, origquo);
  csub(quo2p, quo2p, aux);

  cclear(epsilon); cclear(tmp); cclear(aux); cclear(origquo);
  return;
}

/*
  DupontThetaConstants_Fkprime - computes theta[00,01,10](0,tau) with same precision as tau, assuming tau is in Fkprime
 */
void DupontThetaConstants_Fkprime(mpc_t thetas[3], mpc_t tau){
  int goalprec, prec;
  goalprec = cprec(tau);

  cinit(thetas[0], goalprec); cinit(thetas[1], goalprec); cinit(thetas[2], goalprec);

  mpc_t kprime, res, thet[4];

  // First approximation of kprime
  int j; for (j=0; j<2; j++) cinit(thet[j], LOW_PREC_THETACONSTANTS);
  NaiveThetaConstants_00_01(thet, tau, LOW_PREC_THETACONSTANTS);
  //fprintf(stderr, "1st approx of theta01 : \n"); mpc_out_str(stderr, 10, 0, thet[1], MPC_RNDNN); fprintf(stderr, "\n");
  //fprintf(stderr, "1st approx of theta00 : \n"); mpc_out_str(stderr, 10, 0, thet[0], MPC_RNDNN); fprintf(stderr, "\n");
  int myLowPrec = cprec(thet[1]);
  cinit(kprime, myLowPrec);                  cinit(res, myLowPrec);
  cdiv(kprime, thet[1], thet[0]);
  csqr(kprime, kprime);
  int i; for (i=0; i<2; i++) cclear(thet[i]);
  //fprintf(stderr, "1st approx of kprime : \n"); mpc_out_str(stderr, 10, 0, kprime, MPC_RNDNN); fprintf(stderr, "\n");

  // Newton
  prec = cprec(kprime);
  while (prec < goalprec)
  { prec = 2*prec;
    cset_prec(res, prec);
    DupontNewtonStep(res, kprime, tau);
    cset_prec(kprime, prec); cset(kprime, res);
    //fprintf(stderr, "kprime : \n"); mpc_out_str(stderr, 10, 0, kprime, MPC_RNDNN); fprintf(stderr, "\n");
  }
  cset_prec(kprime, goalprec); cset(kprime, res);  cclear(res);
  //fprintf(stderr, "kprime : \n"); mpc_out_str(stderr, 10, 0, kprime, MPC_RNDNN); fprintf(stderr, "\n");

  // Extract thetas
  DupontM(thetas[0], kprime);
  cinv(thetas[0], thetas[0]); // theta00^2
  cmul(thetas[1], thetas[0], kprime); // theta01^2

  cinit(res, goalprec); cset(res, thetas[1]); csqr(res, res);
  csqr(thetas[2], thetas[0]); csub(thetas[2], thetas[2], res);

  csqrt(thetas[0], thetas[0]); csqrt(thetas[1], thetas[1]);
  csqrt(thetas[2], thetas[2]); csqrt(thetas[2], thetas[2]);
  
  cclear(res); cclear(kprime);
}



/*
 * IsInFkprime - returns 1 if tau is in Fkprime, 0 if it isn't
 */

int IsInFkprime(mpc_t tau) {
  mpc_t unquart, undemi, buf;
  mpfr_t norm;

  int prec = cprec(tau);
  int result = 1;
  
  if ( fcmp_d(cimag(tau), 0.0) <= 0) {result = 0;}
  if ( (fcmp_d(creal(tau), -1.0) < 0) || (fcmp_d(creal(tau), 1.0) >= 0) ) {result = 0;}

  cinit(unquart, prec); cone(unquart); cdiv_2ui(unquart, unquart, 2);
  cinit(undemi, prec); cmul_2ui(undemi, unquart, 1);
  finit(norm, prec);
  cinit(buf, prec); cset(buf, tau);

  csub(buf, buf, undemi); csub(buf, buf, unquart);
  cabs(norm, buf);
      if ( fcmp(norm, creal(unquart)) <= 0) {result = 0;}
  cadd(buf, buf, undemi); cabs(norm, buf);
      if ( fcmp(norm, creal(unquart)) < 0) {result = 0;}
  cadd(buf, buf, undemi); cabs(norm, buf);
      if ( fcmp(norm, creal(unquart)) <= 0) {result = 0;}
  cadd(buf, buf, undemi); cabs(norm, buf);
      if ( fcmp(norm, creal(unquart)) < 0) {result = 0;}
 
  cclear(unquart); cclear(undemi); cclear(buf); fclear(norm);
  return result;
}


/*
 * ReduceInFkprime - Transform tau in tau' such that tau' is in Fkprime
 * Returns 0 if it didn't do anything ; if 1 is returned, you have matrix*tau' = tau
 */

int ReduceInFkprime(mpc_t tauprime, mpfr_t matrix[4], mpc_t tau)
{
  int prec = cprec(tau);
  mpfr_t round, tmp;

  cinit(tauprime, prec); cset(tauprime, tau);
  int i; for (i=0; i<4; i++) { finit(matrix[i], prec); }
  fone(matrix[0]);  fzero(matrix[1]);
  fzero(matrix[2]); fone(matrix[3]);
  finit(round, prec);
  finit(tmp, prec);

  int tauInFkprime=1;
  // Reduce tau so that it is in Fkprime
  //    We first reduce it so that |Re(tau)|<1, then do "-1/tau then |Re(z)|<1"
  while (IsInFkprime(tauprime) == 0) {
    if (tauInFkprime==0) {
      // -1/tau
      cinv(tauprime, tauprime); cneg(tauprime, tauprime);
      fset(tmp, matrix[0]);     fneg(matrix[0], matrix[2]);     fset(matrix[2],tmp);
      fset(tmp, matrix[1]);     fneg(matrix[1], matrix[3]);     fset(matrix[3],tmp);
    }
    tauInFkprime=0;
    // several times \tau -> \tau-1
    frint(round, creal(tauprime)); fneg(round, round);
    cadd_fr(tauprime, tauprime, round);
    fmul(tmp, matrix[2], round); fadd(matrix[0], matrix[0], tmp);
    fmul(tmp, matrix[3], round); fadd(matrix[1], matrix[1], tmp);

    // Check that atauprime+b/ctauprime+d = tau
    /*fprintf(stderr, "matrix[0] : \n"); mpfr_out_str(stderr, 10, 0, matrix[0], MPC_RNDNN); fprintf(stderr, "\n");
    fprintf(stderr, "matrix[1] : \n"); mpfr_out_str(stderr, 10, 0, matrix[1], MPC_RNDNN); fprintf(stderr, "\n");
    fprintf(stderr, "matrix[2] : \n"); mpfr_out_str(stderr, 10, 0, matrix[2], MPC_RNDNN); fprintf(stderr, "\n");
    fprintf(stderr, "matrix[3] : \n"); mpfr_out_str(stderr, 10, 0, matrix[3], MPC_RNDNN); fprintf(stderr, "\n");
    cinit(zeta, prec); cinit(temp, prec);
    cmul_fr(temp, tau, matrix[0]); cadd_fr(temp, temp, matrix[1]);
    cmul_fr(zeta, tau, matrix[2]); cadd_fr(zeta, zeta, matrix[3]); cdiv(temp, temp, zeta);
    fprintf(stderr, "temp : \n"); mpc_out_str(stderr, 10, 0, temp, MPC_RNDNN); fprintf(stderr, "\n");
    fprintf(stderr, "tauprime : \n"); mpc_out_str(stderr, 10, 0, tauprime, MPC_RNDNN); fprintf(stderr, "\n");
    cclear(zeta); cclear(temp); */
  }

  return 1-tauInFkprime;
}


/*
 *  PermutationOfThetas - find the sigma such that Theta_{\sigma(i)}(0, atau+b/ctau+d) = c times Theta(i)
 */

void PermutationOfThetas(int permut[3], mpfr_t matrix[4])
{  // page 53 of Dupont : he does it for theta constants, but it's independent of z and it works for z=0 so...
  int prec = fprec(matrix[1]);
  mpfr_t tmp; finit(tmp, prec);
  int matmod[4], i;

  for (i=0; i<4; i++)
  { fdiv_2ui(tmp, matrix[i], 1); frint(tmp, tmp); fmul_2ui(tmp, tmp, 1);
    if (fcmp(tmp, matrix[i]) == 0) { matmod[i] = 0; } else { matmod[i] = 1; }
  }
  /*fprintf(stderr, "matrix[0] : \n"); mpfr_out_str(stderr, 10, 0, matrix[0], MPC_RNDNN); fprintf(stderr, "\n");
  fprintf(stderr, "matrix[1] : \n"); mpfr_out_str(stderr, 10, 0, matrix[1], MPC_RNDNN); fprintf(stderr, "\n");
  fprintf(stderr, "matrix[2] : \n"); mpfr_out_str(stderr, 10, 0, matrix[2], MPC_RNDNN); fprintf(stderr, "\n");
  fprintf(stderr, "matrix[3] : \n"); mpfr_out_str(stderr, 10, 0, matrix[3], MPC_RNDNN); fprintf(stderr, "\n");
  fprintf(stderr, "matmod : %d %d %d %d", matmod[0], matmod[1], matmod[2], matmod[3]);*/

  // Dupont's formulae
/*  permut[0] = ((matmod[0]*matmod[2]*(matmod[3] + matmod[1]) ) %2)*2 + ((matmod[3]*matmod[1]*(matmod[0]+matmod[2])) %2);
  permut[1] = ((matmod[2] + matmod[0]*matmod[2]*(matmod[3] + matmod[1]) ) %2)*2 + ((matmod[3] + matmod[3]*matmod[1]*(matmod[0]+matmod[2])) %2);
  permut[2] = ((matmod[0] + matmod[0]*matmod[2]*(matmod[3] + matmod[1]) ) %2)*2 + ((matmod[1] + matmod[3]*matmod[1]*(matmod[0]+matmod[2])) %2);
  // this gives the permutation for theta_i = sqrt(ctau+d) theta_{\sigma(i)} - so we want the inverse of that permutation
  // we take the fifth power cause it's more compact, but we should really be better than this
  int permut2[3];  for (i=0; i<3; i++) {permut2[i]=permut[i];}
  for (i=0; i<3; i++) { permut[i] = permut2[permut2[permut2[permut2[permut2[i]]]]]; }
  fprintf(stderr, "permut : %d %d %d \n", permut[0], permut[1], permut[2]);*/

  // Our lookup table
  int matmodint = (1-matmod[0])*8 + (1-matmod[1])*4 + (1-matmod[2])*2 + (1-matmod[3]);
  switch (matmodint)
  {  case 1: permut[0]=1; permut[1]=2; permut[2]=0; break;        // inverse of {2,0,1}
     case 2: permut[0]=1; permut[1]=0; permut[2]=2; break;
     case 4: permut[0]=2; permut[1]=1; permut[2]=0; break;
     case 6: permut[0]=0; permut[1]=1; permut[2]=2; break;
     case 8: permut[0]=2; permut[1]=0; permut[2]=1; break;        // inverse of {1,2,0}
     case 9: permut[0]=0; permut[1]=2; permut[2]=1; break;
     default: fprintf(stderr, "BUG in permutation: contact the authors with your z and your tau\n"); break;
  }

  //fprintf(stderr, "permut : %d %d %d \n", permut[0], permut[1], permut[2]);

  fclear(tmp);
}





/*
  DupontThetaConstants - computes theta[00,01](0,tau) with same precision as tau
 */
void DupontThetaConstants(mpc_t theta[3], mpc_t tau){
  int prec = cprec(tau);

  mpc_t tauprime;
  mpfr_t matrix[4];
  int i;


  int tauInFkprime=ReduceInFkprime(tauprime, matrix, tau);  

  // Call the function
  DupontThetaConstants_Fkprime(theta, tauprime);

  if (tauInFkprime==1) {
    // Retrieve the right result
    mpc_t approx[6], temp, eighthroot;
    cinit(temp, prec);

    // sqrt(c tau + d)
    cmul_fr(temp, tau, matrix[2]); cadd_fr(temp, temp, matrix[3]); csqrt(temp, temp);
    cinv(temp, temp);
    cmul(theta[0], theta[0], temp);
    cmul(theta[1], theta[1], temp);
    cmul(theta[2], theta[2], temp);

    // the right permutation of theta_i (page 53 of Dupont)
    int permut[3];
    PermutationOfThetas(permut, matrix);
    //fprintf(stderr, "permut : %d %d %d \n", permut[0], permut[1], permut[2]);

    //fprintf(stderr, "before :\n"); for (i=0; i<3; i++) { mpc_out_str(stderr, 10, 0, theta[i], MPC_RNDNN); } fprintf(stderr, "\n");
    mpc_t thetabuf[3]; for (i=0; i<3; i++) { cinit(thetabuf[i], prec); cset(thetabuf[i], theta[i]); }
    for (i=0; i<3; i++) { cset(theta[i], thetabuf[permut[i]]); }
    //fprintf(stderr, "after :\n"); for (i=0; i<3; i++) { mpc_out_str(stderr, 10, 0, theta[i], MPC_RNDNN); } fprintf(stderr, "\n");
    for (i=0; i<3; i++) { cclear(thetabuf[i]); }


    cinit(eighthroot, prec);
    // determine and remove the eighth root of unity
    //     (careful ! Dupont's PhD thesis page 53 is wrong: the formula for T forgets an "i" for the case j=2,
    //                so the eighthroot is NOT the same for everyone!)
    czero(temp);
    NaiveThetaFunctionAndThetaConstants_00_01_10(approx, temp, tau, LOW_PREC_EIGHTHROOT);
    //fprintf(stderr, "approx :\n"); for (i=0; i<3; i++) { mpc_out_str(stderr, 10, 0, approx[2*i+1], MPC_RNDNN); } fprintf(stderr, "\n");
    GetEighthRoot(eighthroot, theta[0], approx[1]);
    cdiv(theta[0], theta[0], eighthroot);
    GetEighthRoot(eighthroot, theta[1], approx[3]);
    cdiv(theta[1], theta[1], eighthroot);
    GetEighthRoot(eighthroot, theta[2], approx[5]);
    cdiv(theta[2], theta[2], eighthroot);


    for (i=0; i<6; i++) { cclear(approx[i]); }
    cclear(temp); cclear(eighthroot);
  }

  // clears
  cclear(tauprime);
  for (i=0; i<4; i++) { fclear(matrix[i]); }
}


/*
  GetEighthRoot - given z up to an eighth root of unity and an approx of the final result, recognize the eighth root
*/

void GetEighthRoot(mpc_t eighth, mpc_t z, mpc_t approx)
{
    //fprintf(stderr, "theta : \n"); mpc_out_str(stderr, 10, 0, z, MPC_RNDNN); fprintf(stderr, "\n");
    //fprintf(stderr, "approx : \n"); mpc_out_str(stderr, 10, 0, approx, MPC_RNDNN); fprintf(stderr, "\n");
    int prec = cprec(z);
    cinit(eighth, prec);

    mpc_t zeta;
    mpfr_t tmp;
    cinit(zeta, LOW_PREC_EIGHTHROOT);
    finit(tmp, LOW_PREC_EIGHTHROOT);

    cdiv(zeta, z, approx);
    //fprintf(stderr, "zeta : \n"); mpc_out_str(stderr, 10, 0, zeta, MPC_RNDNN); fprintf(stderr, "\n");

    cinit(eighth, prec); cone(eighth); cmul_by_i(eighth, eighth); csqrt(eighth, eighth);    // e^{i pi / 4}

    int zeta_pow;
    // Look at the first digit of Re() and Im()
    fabs(tmp, creal(zeta)); fmul_ui(tmp, tmp, 10); frint(tmp, tmp);
    if ( fcmp_ui(tmp, 7) == 0) {  // 0.7... : it's 1,3,5 or 7
      if (fcmp_ui(creal(zeta), 0) < 0) {
        if (fcmp_ui(cimag(zeta), 0) < 0) { zeta_pow = 5; } else { zeta_pow = 3; }
      } else {
        if (fcmp_ui(cimag(zeta), 0) < 0) { zeta_pow = 7; } else { zeta_pow = 1; }
      }
    } else { // it's 0,2,4 or 6
      // turn everything by 45° to have better tests
      cmul(zeta, zeta, eighth);
      if (fcmp_ui(creal(zeta), 0) < 0) {
        if (fcmp_ui(cimag(zeta), 0) < 0) { zeta_pow = 5; } else { zeta_pow = 3; }
      } else {
        if (fcmp_ui(cimag(zeta), 0) < 0) { zeta_pow = 7; } else { zeta_pow = 1; }
      }
      zeta_pow = zeta_pow-1;
    }
    //fprintf(stderr, "zeta_pow : %d\n", zeta_pow);
    cpow_d(eighth, eighth, (double) zeta_pow);

    cclear(zeta); fclear(tmp);
}

