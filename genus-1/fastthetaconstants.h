#ifndef FAST_THETA_CONSTANTS_H_
#define FAST_THETA_CONSTANTS_H_

/* fastthetaconstants.h -- headers for fastthetaconstants.c
 *
 */


#ifdef __cplusplus
extern "C" {
#endif


// Dupont's PhD says it's faster than naive for P>2500
#ifndef LOW_PREC_THETACONSTANTS
#define LOW_PREC_THETACONSTANTS 2500
#endif

// Prec needed to recognize the eighthroot (very small is ok)
#ifndef LOW_PREC_EIGHTHROOT
#define LOW_PREC_EIGHTHROOT 50
#endif

int IsInFkprime(mpc_t tau);
int ReduceInFkprime(mpc_t tauprime, mpfr_t matrix[4], mpc_t tau);
void PermutationOfThetas(int permut[3], mpfr_t matrix[4]);
void GetEighthRoot(mpc_t eighth, mpc_t z, mpc_t approx);

void DupontOneStepOfAGM(mpc_t a[2]);
void DupontM(mpc_t res, mpc_t quo);
void DupontFtau(mpc_t res, mpc_t arg, mpc_t tau);
void DupontNewtonStep(mpc_t quo2p, mpc_t quop, mpc_t tau);
void DupontThetaConstants_Fkprime(mpc_t thetas[3], mpc_t tau);
void DupontThetaConstants(mpc_t thetas[3], mpc_t tau);


#ifdef __cplusplus
}
#endif


#endif
